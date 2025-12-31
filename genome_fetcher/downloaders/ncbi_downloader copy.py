"""
NCBI downloader for resolving proteins to genomes
"""
import time
import logging
import requests
from pathlib import Path
from typing import Optional, Tuple
from genome_fetcher.config import NCBI_BASE_URL, NCBI_EMAIL, NCBI_TOOL, NCBI_RATE_LIMIT
from genome_fetcher.utils import DownloadLogger, ResponseCache


logger = logging.getLogger(__name__)


class NCBIDownloader:
    """Download and resolve genome information from NCBI"""
    
    # def __init__(self, email: str = NCBI_EMAIL):
    def __init__(self, email: str = NCBI_EMAIL, use_cache: bool = True):
        """
        Initialize NCBI downloader
        
        Args:
            email: Email for NCBI E-utilities (required by NCBI)
        """
        from threading import Lock
        self.request_lock = Lock()
        self.email = email
        self.tool = NCBI_TOOL
        self.last_request_time = 0
        self.rate_limit = NCBI_RATE_LIMIT
        self.session = requests.Session()
        self.cache = ResponseCache(enabled=use_cache)  # ← ADD THIS

    
    def _sanitize_json_response(self, text: str) -> str:
        """
        Sanitize JSON response by removing invalid control characters
        
        NCBI sometimes returns responses with control characters that break JSON parsing.
        This function removes them while preserving valid content.
        """
        # Remove control characters except newline, carriage return, and tab
        # Keep characters: \n (10), \r (13), \t (9)
        sanitized = ''.join(
            char for char in text 
            if ord(char) >= 32 or ord(char) in [9, 10, 13]
        )
        return sanitized
    
    # def _rate_limited_request(self, url: str, params: dict = None) -> Optional[requests.Response]:
    #     """Make a rate-limited request to NCBI with robust error handling"""
    #     # Rate limiting
    #     current_time = time.time()
    #     time_since_last = current_time - self.last_request_time
        
    #     if time_since_last < self.rate_limit:
    #         time.sleep(self.rate_limit - time_since_last)
        
    #     self.last_request_time = time.time()
        
    #     # Add email and tool to params
    #     if params is None:
    #         params = {}
    #     params['email'] = self.email
    #     params['tool'] = self.tool
        
    #     try:
    #         # response = requests.get(url, params=params, timeout=150)
    #         response = self.session.get(url, params=params, timeout=150)

    #         response.raise_for_status()
    #         return response
    #     except requests.exceptions.Timeout:
    #         logger.error(f"NCBI request timeout for {url}")
    #         return None
    #     except requests.exceptions.HTTPError as e:
    #         logger.error(f"NCBI HTTP error: {e}")
    #         return None
    #     except Exception as e:
    #         logger.error(f"NCBI request failed: {e}")
    #         return None




    # def _rate_limited_request(self, url: str, params: dict = None, max_retries: int = 3) -> Optional[requests.Response]:
    #     """Make a rate-limited request to NCBI with retry logic"""
        
    #     if params is None:
    #         params = {}
    #     params['email'] = self.email
    #     params['tool'] = self.tool
        
    #     for attempt in range(max_retries):
    #         # Rate limiting
    #         current_time = time.time()
    #         time_since_last = current_time - self.last_request_time
            
    #         if time_since_last < self.rate_limit:
    #             time.sleep(self.rate_limit - time_since_last)
            
    #         self.last_request_time = time.time()
            
    #         try:
    #             response = self.session.get(url, params=params, timeout=150)
    #             response.raise_for_status()
    #             return response
                
    #         except requests.exceptions.Timeout:
    #             logger.warning(f"NCBI timeout (attempt {attempt + 1}/{max_retries})")
    #             if attempt < max_retries - 1:
    #                 time.sleep(2 ** attempt)  # 1s, 2s, 4s
    #                 continue
    #             logger.error(f"NCBI request timeout after {max_retries} attempts")
    #             return None
                
    #         except requests.exceptions.HTTPError as e:
    #             if e.response.status_code == 429:  # Rate limit
    #                 logger.warning(f"NCBI rate limit hit (attempt {attempt + 1}/{max_retries})")
    #                 time.sleep(5 * (attempt + 1))  # 5s, 10s, 15s
    #                 continue
    #             elif e.response.status_code >= 500:  # Server error
    #                 logger.warning(f"NCBI server error {e.response.status_code} (attempt {attempt + 1}/{max_retries})")
    #                 if attempt < max_retries - 1:
    #                     time.sleep(2 ** attempt)
    #                     continue
    #             logger.error(f"NCBI HTTP error: {e}")
    #             return None
                
    #         except Exception as e:
    #             logger.error(f"NCBI request failed: {e}")
    #             return None
        
    #     return None


    def _rate_limited_request(self, url: str, params: dict = None, max_retries: int = 3) -> Optional[requests.Response]:
        """Make a rate-limited request to NCBI with caching and retry logic"""
        
        # ADDED FOR SAFE MULTIPROCESSING
        with self.request_lock:  # ← ADD THIS
            # Rate limiting
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            
            if time_since_last < self.rate_limit:
                time.sleep(self.rate_limit - time_since_last)
            
            self.last_request_time = time.time()

        if params is None:
            params = {}
        params['email'] = self.email
        params['tool'] = self.tool
        
        # Build full URL for caching
        from urllib.parse import urlencode
        cache_url = f"{url}?{urlencode(params)}"
        
        # Check cache first
        if self.cache.enabled:
            cached_text = self.cache.get(cache_url)
            if cached_text:
                from requests import Response
                response = Response()
                response._content = cached_text.encode('utf-8')
                response.status_code = 200
                return response
        
        # Make request with retry logic (rest of your existing retry code)
        for attempt in range(max_retries):
            # Rate limiting
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            
            if time_since_last < self.rate_limit:
                time.sleep(self.rate_limit - time_since_last)
            
            self.last_request_time = time.time()
            
            try:
                response = self.session.get(url, params=params, timeout=150)
                response.raise_for_status()
                
                # Cache successful response
                if self.cache.enabled:
                    self.cache.set(cache_url, response.text)
                
                return response
            except requests.exceptions.Timeout:
                logger.warning(f"NCBI timeout (attempt {attempt + 1}/{max_retries})")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                logger.error(f"NCBI request timeout after {max_retries} attempts")
                return None
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 429:
                    logger.warning(f"NCBI rate limit hit (attempt {attempt + 1}/{max_retries})")
                    time.sleep(5 * (attempt + 1))
                    continue
                elif e.response.status_code >= 500:
                    logger.warning(f"NCBI server error {e.response.status_code} (attempt {attempt + 1}/{max_retries})")
                    if attempt < max_retries - 1:
                        time.sleep(2 ** attempt)
                        continue
                logger.error(f"NCBI HTTP error: {e}")
                return None
            except Exception as e:
                logger.error(f"NCBI request failed: {e}")
                return None
        
        return None

    
    def _parse_json_response(self, response: requests.Response) -> Optional[dict]:
        """
        Safely parse JSON response with sanitization
        
        Args:
            response: requests.Response object
            
        Returns:
            Parsed JSON dict or None if parsing fails
        """
        try:
            # First try standard JSON parsing
            return response.json()
        except json.JSONDecodeError as e:
            logger.warning(f"Standard JSON parsing failed: {e}")
            logger.info("Attempting to sanitize response and retry...")
            
            try:
                # Sanitize and retry
                sanitized_text = self._sanitize_json_response(response.text)
                return json.loads(sanitized_text)
            except json.JSONDecodeError as e2:
                logger.error(f"JSON parsing failed even after sanitization: {e2}")
                logger.error(f"Response preview: {response.text[:200]}...")
                return None
        except Exception as e:
            logger.error(f"Unexpected error parsing JSON: {e}")
            return None


    # NOT BEING USED, leave for later
    # def _download_large_file_streaming(self, url: str, output_file: Path, params: dict = None) -> bool:
    #     """
    #     Download large file with streaming to avoid memory issues
        
    #     Args:
    #         url: URL to download
    #         output_file: Path to save file
    #         params: Optional query parameters
            
    #     Returns:
    #         True if successful
    #     """
    #     try:
    #         if params is None:
    #             params = {}
    #         params['email'] = self.email
    #         params['tool'] = self.tool
            
    #         # Rate limiting
    #         current_time = time.time()
    #         time_since_last = current_time - self.last_request_time
    #         if time_since_last < self.rate_limit:
    #             time.sleep(self.rate_limit - time_since_last)
    #         self.last_request_time = time.time()
            
    #         # Stream the download
    #         response = self.session.get(url, params=params, timeout=300, stream=True)
    #         response.raise_for_status()
            
    #         with open(output_file, 'wb') as f:
    #             for chunk in response.iter_content(chunk_size=8192):
    #                 if chunk:
    #                     f.write(chunk)
            
    #         return True
    #     except Exception as e:
    #         logger.error(f"Streaming download failed: {e}")
    #         return False
    
    def resolve_protein_to_genome(
        self,
        protein_accession: str,
        terminal_log: DownloadLogger
    ) -> Optional[Tuple[str, str]]:
        """
        Resolve a protein accession to its source genome
        
        Args:
            protein_accession: NCBI protein accession (e.g., "WP_000111473.1")
            terminal_log: Logger for this terminal
        
        Returns:
            Tuple of (genome_accession, genome_taxid) or None if failed
        """
        terminal_log.info(f"Resolving protein {protein_accession} to genome via NCBI")
        logger.info(f"Resolving protein {protein_accession} to genome")
        
        # Step 1: Get protein record using esummary
        esummary_url = f"{NCBI_BASE_URL}/esummary.fcgi"
        params = {
            'db': 'protein',
            'id': protein_accession,
            'retmode': 'json'
        }
        
        response = self._rate_limited_request(esummary_url, params)
        if not response:
            terminal_log.error("Failed to get protein summary from NCBI")
            return None
        
        # Parse JSON with sanitization
        data = self._parse_json_response(response)
        if not data:
            terminal_log.error("Failed to parse NCBI response")
            return None
        
        try:
            # Extract information from the response
            result = data.get('result', {})
            uids = result.get('uids', [])
            
            if not uids:
                terminal_log.error(f"No protein found for accession {protein_accession}")
                return None
            
            protein_info = result.get(uids[0], {})
            
            # Try to get genome/assembly information
            # Method 1: Check for genome accession in protein record
            genome_acc = protein_info.get('genome', '')
            taxid = str(protein_info.get('taxid', ''))
            
            if genome_acc and taxid:
                terminal_log.success(f"Resolved to genome {genome_acc}, taxid {taxid}")
                logger.info(f"Resolved {protein_accession} to genome {genome_acc}")
                return (genome_acc, taxid)
            
            # Method 2: Use elink to find linked genomes
            terminal_log.info("Attempting to find linked genome via elink...")
            elink_url = f"{NCBI_BASE_URL}/elink.fcgi"
            link_params = {
                'dbfrom': 'protein',
                'db': 'nuccore',
                'id': protein_accession,
                'retmode': 'json'
            }
            
            link_response = self._rate_limited_request(elink_url, link_params)
            if link_response:
                link_data = self._parse_json_response(link_response)
                if link_data:
                    linksets = link_data.get('linksets', [])
                    
                    if linksets and 'linksetdbs' in linksets[0]:
                        for linkdb in linksets[0]['linksetdbs']:
                            if linkdb.get('linkname') == 'protein_nuccore':
                                links = linkdb.get('links', [])
                                if links:
                                    genome_id = str(links[0])
                                    terminal_log.success(f"Found linked genome ID: {genome_id}")
                                    logger.info(f"Found genome {genome_id} for protein {protein_accession}")
                                    
                                    # Get taxid if we don't have it
                                    if not taxid:
                                        taxid = self.get_genome_taxid(genome_id) or ''
                                    
                                    return (genome_id, taxid)
            
            terminal_log.warning(f"Could not resolve protein {protein_accession} to a genome")
            logger.warning(f"Failed to resolve protein {protein_accession} to genome")
            return None
            
        except Exception as e:
            terminal_log.error(f"Error parsing NCBI response: {e}")
            logger.error(f"Error resolving protein {protein_accession}: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def get_genome_taxid(self, genome_accession: str) -> Optional[str]:
        """
        Get taxonomy ID for a genome accession
        
        Args:
            genome_accession: NCBI genome/nuccore accession
        
        Returns:
            Taxonomy ID or None
        """
        esummary_url = f"{NCBI_BASE_URL}/esummary.fcgi"
        params = {
            'db': 'nuccore',
            'id': genome_accession,
            'retmode': 'json'
        }
        
        response = self._rate_limited_request(esummary_url, params)
        if not response:
            return None
        
        data = self._parse_json_response(response)
        if not data:
            return None
        
        try:
            result = data.get('result', {})
            uids = result.get('uids', [])
            
            if uids:
                info = result.get(uids[0], {})
                return str(info.get('taxid', ''))
        except Exception as e:
            logger.error(f"Error getting taxid for {genome_accession}: {e}")
        
        return None
    


    def download_genome_files(
        self,
        genome_accession: str,
        output_dir: Path,
        terminal_log: DownloadLogger
    ) -> bool:


        """
        Download genome files from NCBI
        
        Args:
            genome_accession: NCBI genome accession
            output_dir: Output directory
            terminal_log: Logger for this terminal
            
        Returns:
            True if successful
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        terminal_log.info(f"Downloading genome {genome_accession} from NCBI")
        logger.info(f"Downloading genome {genome_accession}")
        
        genome_fasta = output_dir / "genome.fasta"
        
        # Check if already exists
        if genome_fasta.exists() and genome_fasta.stat().st_size > 100:
            terminal_log.info("Genome FASTA already exists")
            return True
        
        # Download genome sequence
        # Download genome sequence
        efetch_url = f"{NCBI_BASE_URL}/efetch.fcgi"
        params = {
            'db': 'nuccore',
            'id': genome_accession,
            'rettype': 'fasta',
            'retmode': 'text'
        }
        
        response = self._rate_limited_request(efetch_url, params)  # ← Current approach
        if not response:
            terminal_log.error("Failed to download genome from NCBI")
            return False
        
        try:
            content = response.text  # ← Loads entire response into memory
            if content and len(content) > 100:
                # Convert to uppercase and write
                self._write_fasta_uppercase(content, genome_fasta)
                terminal_log.success(f"Downloaded genome.fasta ({len(content)} bytes)")
                return True
            else:
                terminal_log.error("Genome sequence response empty or too small")
                return False
        except Exception as e:
            terminal_log.error(f"Error saving genome: {e}")
            return False
    
    def _write_fasta_uppercase(self, fasta_content: str, output_file: Path):
        """
        Write FASTA content with sequences in uppercase
        
        Args:
            fasta_content: FASTA formatted text
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            for line in fasta_content.split('\n'):
                if line.startswith('>'):
                    # Keep header as-is
                    f.write(line + '\n')
                else:
                    # Convert sequence to uppercase
                    f.write(line.upper() + '\n')
    
    def download_protein_sequences(
        self,
        protein_accession: str,
        output_dir: Path,
        terminal_log: DownloadLogger,
        sequence_id: str = None
    ) -> bool:
        """
        Download protein amino acid and nucleotide sequences
        
        Args:
            protein_accession: Protein accession
            output_dir: Output directory
            terminal_log: Logger for this terminal
            sequence_id: Optional custom sequence ID for output
            
        Returns:
            True if successful
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if sequence_id is None:
            sequence_id = protein_accession
        
        terminal_log.info(f"Downloading protein sequences for {protein_accession}")
        
        protein_aa = output_dir / "protein_aminoacid.fasta"
        protein_nt = output_dir / "protein_nucleotide.fasta"
        
        success = True
        
        # Download amino acid sequence
        if not protein_aa.exists() or protein_aa.stat().st_size < 50:
            efetch_url = f"{NCBI_BASE_URL}/efetch.fcgi"
            params = {
                'db': 'protein',
                'id': protein_accession,
                'rettype': 'fasta',
                'retmode': 'text'
            }
            
            response = self._rate_limited_request(efetch_url, params)
            if response and response.text:
                content = response.text
                # Replace header with custom sequence_id
                lines = content.split('\n')
                with open(protein_aa, 'w') as f:
                    f.write(f">{sequence_id}\n")
                    for line in lines[1:]:
                        if line.strip():
                            # Convert to uppercase
                            f.write(line.strip().upper() + '\n')
                terminal_log.success(f"Downloaded protein amino acid sequence")
            else:
                terminal_log.warning("Failed to download amino acid sequence")
                success = False
        
        # Download nucleotide sequence (CDS)
        if not protein_nt.exists() or protein_nt.stat().st_size < 50:
            efetch_url = f"{NCBI_BASE_URL}/efetch.fcgi"
            params = {
                'db': 'protein',
                'id': protein_accession,
                'rettype': 'fasta_cds_na',
                'retmode': 'text'
            }
            
            response = self._rate_limited_request(efetch_url, params)
            if response and response.text:
                content = response.text
                # Replace header with custom sequence_id
                lines = content.split('\n')
                with open(protein_nt, 'w') as f:
                    f.write(f">{sequence_id}\n")
                    for line in lines[1:]:
                        if line.strip():
                            # Convert to uppercase
                            f.write(line.strip().upper() + '\n')
                terminal_log.success(f"Downloaded protein nucleotide sequence")
            else:
                terminal_log.warning("Failed to download nucleotide sequence")
                success = False
        
        return success