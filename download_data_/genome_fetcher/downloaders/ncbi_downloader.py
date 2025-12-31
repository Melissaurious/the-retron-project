# """
# NCBI downloader for resolving proteins to genomes
# """
# import time
# import logging
# import requests
# from pathlib import Path
# from typing import Optional, Tuple
# from genome_fetcher.config import NCBI_BASE_URL, NCBI_EMAIL, NCBI_TOOL, NCBI_RATE_LIMIT
# from genome_fetcher.utils import DownloadLogger, ResponseCache


# logger = logging.getLogger(__name__)

"""
Enhanced NCBI downloader with BaseDownloader interface.

This preserves all original functionality (E-utilities, protein resolution, caching)
while adding the BaseDownloader interface for consistency with other downloaders.
"""

import logging
import time
import json
import gzip
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from threading import Lock

# Assuming these exist in your project
try:
    from .base_downloader import BaseDownloader
except ImportError:
    # Fallback if running standalone
    import sys
    sys.path.append(str(Path(__file__).parent))
    from base_downloader import BaseDownloader

logger = logging.getLogger(__name__)

# These should come from your config
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_RATE_LIMIT = 0.34  # ~3 requests per second


class NCBIDownloader(BaseDownloader):
    """
    Enhanced NCBI downloader with full E-utilities support.
    
    Features:
    - Protein → genome resolution via E-utilities
    - Genome and protein sequence downloads
    - Response caching to reduce API calls
    - Rate limiting for API compliance
    - Streaming downloads for large files
    - Inherits BaseDownloader interface for consistency
    """
    
    GENBANK_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
    REFSEQ_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
    
    def __init__(self, email: str = None, use_cache: bool = True, database_name: str = "ncbi"):
        """
        Initialize NCBI downloader.
        
        Args:
            email: Email for NCBI E-utilities (required by NCBI)
            use_cache: Enable response caching
            database_name: Database name for BaseDownloader
        """
        super().__init__(database_name)
        
        self.email = email or "your_email@example.com"
        self.tool = "GenomeFetcher"
        self.last_request_time = 0
        self.rate_limit = NCBI_RATE_LIMIT
        self.session = requests.Session()
        self.request_lock = Lock()
        
        # Simple in-memory cache
        self.cache_enabled = use_cache
        self.cache = {}
        
        logger.info(f"Initialized NCBI downloader (email: {self.email}, cache: {use_cache})")
    
    # =========================================================================
    # BaseDownloader Interface Methods
    # =========================================================================
    
    def download_genome(
        self,
        genome_record: Dict,
        output_dir: Path,
        required_files: List[str] = None,
        compress: bool = True,
        terminal_log: 'DownloadLogger' = None
    ) -> bool:
        """
        Download genome from NCBI (BaseDownloader interface).
        
        Args:
            genome_record: Genome metadata dictionary
            output_dir: Output directory
            required_files: Files to download (default: ['genome.fasta'])
            compress: Whether to compress output
            terminal_log: Optional logger
        
        Returns:
            True if successful
        """
        if required_files is None:
            required_files = ['genome.fasta']
        
        self.update_stats('attempted')
        
        # Extract genome ID (assembly accession)
        genome_id = self.get_genome_id(genome_record)
        
        if terminal_log:
            terminal_log.info(f"Processing NCBI genome: {genome_id}")
        
        # Validate
        if not self.validate_genome_id(genome_id):
            if terminal_log:
                terminal_log.error(f"Invalid NCBI accession: {genome_id}")
            self.update_stats('failed')
            return False
        
        # Create genome directory
        genome_dir = output_dir / genome_id
        genome_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if complete
        if self.is_genome_complete(genome_dir, required_files):
            if terminal_log:
                terminal_log.info("Genome already complete, skipping")
            self.update_stats('skipped')
            return True
        
        # Download genome
        success = self.download_genome_files(
            genome_id,
            genome_dir,
            terminal_log,
            compress=compress
        )
        
        if success:
            self.update_stats('successful')
            if terminal_log:
                terminal_log.success(f"Downloaded {genome_id}")
        else:
            self.update_stats('failed')
            if terminal_log:
                terminal_log.error(f"Failed to download {genome_id}")
        
        return success
    
    def get_metadata(self, genome_record: Dict) -> Dict:
        """Extract metadata from NCBI record."""
        genome_id = self.get_genome_id(genome_record)
        
        return {
            'genome_id': genome_id,
            'database': self.database_name,
            'assembly_accession': genome_record.get('assembly_accession', genome_id),
            'bioproject': genome_record.get('bioproject'),
            'biosample': genome_record.get('biosample'),
            'organism_name': genome_record.get('organism_name'),
            'assembly_level': genome_record.get('assembly_level'),
            'genome_rep': genome_record.get('genome_rep'),
            'ftp_path': genome_record.get('ftp_path')
        }
    
    def validate_genome_id(self, genome_id: str) -> bool:
        """Validate NCBI assembly accession format."""
        import re
        pattern = r'^(GCA|GCF)_\d{9}\.\d+$'
        return bool(re.match(pattern, genome_id))
    
    # =========================================================================
    # Original Methods (Preserved)
    # =========================================================================
    
    def _sanitize_json_response(self, text: str) -> str:
        """Remove invalid control characters from JSON response."""
        sanitized = ''.join(
            char for char in text 
            if ord(char) >= 32 or ord(char) in [9, 10, 13]
        )
        return sanitized
    
    def _rate_limited_request(self, url: str, params: dict = None, max_retries: int = 3) -> Optional[requests.Response]:
        """Make rate-limited request with retry logic."""
        
        if params is None:
            params = {}
        params['email'] = self.email
        params['tool'] = self.tool
        
        # Build cache key
        from urllib.parse import urlencode
        cache_key = f"{url}?{urlencode(sorted(params.items()))}"
        
        # Check cache
        if self.cache_enabled and cache_key in self.cache:
            logger.debug(f"Cache hit for: {url}")
            response = requests.Response()
            response._content = self.cache[cache_key].encode('utf-8')
            response.status_code = 200
            return response
        
        # Rate limiting
        with self.request_lock:
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            
            if time_since_last < self.rate_limit:
                time.sleep(self.rate_limit - time_since_last)
            
            self.last_request_time = time.time()
        
        # Make request with retries
        for attempt in range(max_retries):
            try:
                response = self.session.get(url, params=params, timeout=150, stream=True)
                response.raise_for_status()
                
                # Cache successful response (but not file downloads)
                if self.cache_enabled and 'rettype=fasta' not in url:
                    self.cache[cache_key] = response.text
                
                return response
                
            except requests.exceptions.Timeout:
                logger.warning(f"NCBI timeout (attempt {attempt + 1}/{max_retries})")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                return None
                
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 429:
                    logger.warning(f"Rate limit hit (attempt {attempt + 1}/{max_retries})")
                    time.sleep(5 * (attempt + 1))
                    continue
                elif e.response.status_code >= 500:
                    logger.warning(f"Server error {e.response.status_code}")
                    if attempt < max_retries - 1:
                        time.sleep(2 ** attempt)
                        continue
                logger.error(f"HTTP error: {e}")
                return None
                
            except Exception as e:
                logger.error(f"Request failed: {e}")
                return None
        
        return None
    
    def _parse_json_response(self, response: requests.Response) -> Optional[dict]:
        """Safely parse JSON response with sanitization."""
        try:
            return response.json()
        except json.JSONDecodeError as e:
            logger.warning(f"JSON parsing failed: {e}, attempting sanitization...")
            try:
                sanitized_text = self._sanitize_json_response(response.text)
                return json.loads(sanitized_text)
            except json.JSONDecodeError as e2:
                logger.error(f"JSON parsing failed after sanitization: {e2}")
                return None
        except Exception as e:
            logger.error(f"Unexpected error parsing JSON: {e}")
            return None
    
    def resolve_protein_to_genome(
        self,
        protein_accession: str,
        terminal_log: 'DownloadLogger'
    ) -> Optional[Tuple[str, str]]:
        """
        Resolve protein accession to source genome.
        
        Args:
            protein_accession: NCBI protein accession
            terminal_log: Logger
        
        Returns:
            Tuple of (genome_accession, taxid) or None
        """
        if terminal_log:
            terminal_log.info(f"Resolving protein {protein_accession} to genome")
        
        logger.info(f"Resolving protein {protein_accession}")
        
        # Get protein summary
        esummary_url = f"{NCBI_BASE_URL}/esummary.fcgi"
        params = {
            'db': 'protein',
            'id': protein_accession,
            'retmode': 'json'
        }
        
        response = self._rate_limited_request(esummary_url, params)
        if not response:
            if terminal_log:
                terminal_log.error("Failed to get protein summary")
            return None
        
        data = self._parse_json_response(response)
        if not data:
            if terminal_log:
                terminal_log.error("Failed to parse NCBI response")
            return None
        
        try:
            result = data.get('result', {})
            uids = result.get('uids', [])
            
            if not uids:
                if terminal_log:
                    terminal_log.error(f"No protein found: {protein_accession}")
                return None
            
            protein_info = result.get(uids[0], {})
            genome_acc = protein_info.get('genome', '')
            taxid = str(protein_info.get('taxid', ''))
            
            if genome_acc and taxid:
                if terminal_log:
                    terminal_log.success(f"Resolved to genome {genome_acc}, taxid {taxid}")
                logger.info(f"Resolved {protein_accession} → {genome_acc}")
                return (genome_acc, taxid)
            
            # Try elink to find linked genomes
            if terminal_log:
                terminal_log.info("Attempting elink to find genome...")
            
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
                                    if terminal_log:
                                        terminal_log.success(f"Found linked genome: {genome_id}")
                                    
                                    if not taxid:
                                        taxid = self.get_genome_taxid(genome_id) or ''
                                    
                                    return (genome_id, taxid)
            
            if terminal_log:
                terminal_log.warning(f"Could not resolve {protein_accession} to genome")
            return None
            
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"Error parsing response: {e}")
            logger.error(f"Error resolving protein {protein_accession}: {e}")
            return None
    
    def get_genome_taxid(self, genome_accession: str) -> Optional[str]:
        """Get taxonomy ID for genome accession."""
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
            logger.error(f"Error getting taxid: {e}")
        
        return None
    
    def download_genome_files(
        self,
        genome_accession: str,
        output_dir: Path,
        terminal_log: 'DownloadLogger',
        compress: bool = True
    ) -> bool:
        """
        Download genome FASTA with streaming to avoid memory issues.
        
        Args:
            genome_accession: NCBI genome accession
            output_dir: Output directory
            terminal_log: Logger
            compress: Whether to compress output
        
        Returns:
            True if successful
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if terminal_log:
            terminal_log.info(f"Downloading genome {genome_accession}")
        logger.info(f"Downloading genome {genome_accession}")
        
        # Determine output file
        if compress:
            genome_file = output_dir / "genome.fasta.gz"
        else:
            genome_file = output_dir / "genome.fasta"
        
        # Check if exists
        if genome_file.exists() and genome_file.stat().st_size > 100:
            if terminal_log:
                terminal_log.info("Genome FASTA already exists")
            return True
        
        # Download via E-fetch
        efetch_url = f"{NCBI_BASE_URL}/efetch.fcgi"
        params = {
            'db': 'nuccore',
            'id': genome_accession,
            'rettype': 'fasta',
            'retmode': 'text'
        }
        
        response = self._rate_limited_request(efetch_url, params)
        if not response:
            if terminal_log:
                terminal_log.error("Failed to download genome")
            return False
        
        try:
            # Stream to disk with uppercase conversion
            if compress:
                self._stream_fasta_uppercase_compressed(response, genome_file)
            else:
                self._stream_fasta_uppercase(response, genome_file)
            
            file_size = genome_file.stat().st_size
            if file_size < 100:
                if terminal_log:
                    terminal_log.error("Downloaded file too small")
                return False
            
            if terminal_log:
                terminal_log.success(f"Downloaded genome.fasta ({file_size:,} bytes)")
            
            return True
            
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"Error saving genome: {e}")
            logger.error(f"Error downloading genome: {e}")
            return False
    
    def _stream_fasta_uppercase(self, response, output_file: Path):
        """Stream FASTA to disk with uppercase conversion."""
        with open(output_file, 'w') as f:
            for line in response.iter_lines(decode_unicode=True):
                if line.startswith('>'):
                    f.write(line + '\n')
                else:
                    f.write(line.upper() + '\n')
    
    def _stream_fasta_uppercase_compressed(self, response, output_file: Path):
        """Stream FASTA to compressed file with uppercase conversion."""
        with gzip.open(output_file, 'wt') as f:
            for line in response.iter_lines(decode_unicode=True):
                if line.startswith('>'):
                    f.write(line + '\n')
                else:
                    f.write(line.upper() + '\n')
    
    def download_protein_sequences(
        self,
        protein_accession: str,
        output_dir: Path,
        terminal_log: 'DownloadLogger',
        sequence_id: str = None
    ) -> bool:
        """
        Download protein amino acid and nucleotide sequences.
        
        Args:
            protein_accession: Protein accession
            output_dir: Output directory
            terminal_log: Logger
            sequence_id: Custom sequence ID
        
        Returns:
            True if successful
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if sequence_id is None:
            sequence_id = protein_accession
        
        if terminal_log:
            terminal_log.info(f"Downloading protein sequences for {protein_accession}")
        
        protein_aa = output_dir / "protein_aminoacid.fasta"
        protein_nt = output_dir / "protein_nucleotide.fasta"
        
        success = True
        
        # Download amino acid
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
                lines = content.split('\n')
                with open(protein_aa, 'w') as f:
                    f.write(f">{sequence_id}\n")
                    for line in lines[1:]:
                        if line.strip():
                            f.write(line.strip().upper() + '\n')
                if terminal_log:
                    terminal_log.success("Downloaded amino acid sequence")
            else:
                if terminal_log:
                    terminal_log.warning("Failed to download amino acid sequence")
                success = False
        
        # Download nucleotide
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
                lines = content.split('\n')
                with open(protein_nt, 'w') as f:
                    f.write(f">{sequence_id}\n")
                    for line in lines[1:]:
                        if line.strip():
                            f.write(line.strip().upper() + '\n')
                if terminal_log:
                    terminal_log.success("Downloaded nucleotide sequence")
            else:
                if terminal_log:
                    terminal_log.warning("Failed to download nucleotide sequence")
                success = False
        
        return success