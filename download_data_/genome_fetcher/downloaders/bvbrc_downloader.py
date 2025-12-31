"""
BV-BRC (PATRIC) downloader for genome files via FTP
"""
import time
import logging
import urllib.request
import socket
import ftplib 
import requests
import hashlib
import pickle

from pathlib import Path
from typing import Optional, Tuple
from genome_fetcher.config import BVBRC_FTP_BASE
from genome_fetcher.utils import (
    DownloadLogger,
    convert_fasta_to_uppercase,
    is_file_complete, ResponseCache
)

logger = logging.getLogger(__name__)

class BVBRCDownloader:
    """Download genome files from BV-BRC Web API"""
    
    # def __init__(self, rate_limit: float = 0.1):
    def __init__(self, rate_limit: float = 0.1, use_cache: bool = True):
        
        """
        Initialize BV-BRC downloader
        
        Args:
            rate_limit: Minimum seconds between API requests
        """

        from threading import Lock
        self.request_lock = Lock()  # ← ADD THIS



        self.rate_limit = rate_limit
        self.last_request_time = 0
        self.base_url = "https://www.bv-brc.org"
        self.session = requests.Session()
        self.cache = ResponseCache(enabled=use_cache)  # ← ADD THIS


    
    def _get_cache_key(self, url: str) -> str:
        """Generate cache key from URL"""
        return hashlib.md5(url.encode()).hexdigest()
    
    def _get_cached_response(self, url: str) -> Optional[str]:
        """Get cached response if available"""
        if not self.use_cache:
            return None
        
        cache_file = self.cache_dir / f"{self._get_cache_key(url)}.cache"
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            except:
                return None
        return None
    
    def _cache_response(self, url: str, response_text: str):
        """Cache successful response"""
        if not self.use_cache:
            return
        
        cache_file = self.cache_dir / f"{self._get_cache_key(url)}.cache"
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(response_text, f)
        except:
            pass
    
    # def _rate_limited_request(
    #     self, 
    #     url: str, 
    #     headers: dict = None, 
    #     **kwargs
    # ) -> Optional[requests.Response]:

    # def _rate_limited_request(self, url: str, headers: dict = None, **kwargs) -> Optional[requests.Response]:
    #     """
    #     Make a rate-limited HTTP request
        
    #     Args:
    #         url: URL to request
    #         headers: Optional headers dictionary
    #         **kwargs: Additional arguments for requests.get()
    #     """

    #     # Rate limiting
    #     current_time = time.time()
    #     time_since_last = current_time - self.last_request_time
        
    #     if time_since_last < self.rate_limit:
    #         time.sleep(self.rate_limit - time_since_last)
        
    #     self.last_request_time = time.time()
        
    #     try:
    #         # Merge headers if provided
    #         if headers:
    #             if 'headers' in kwargs:
    #                 kwargs['headers'].update(headers)
    #             else:
    #                 kwargs['headers'] = headers
            
    #         # response = requests.get(url, timeout=150, **kwargs)
    #         response = self.session.get(url, timeout=120, **kwargs)  
    #         return response
    #     except Exception as e:
    #         logger.warning(f"Request failed for {url}: {e}")
    #         return None




    # def _rate_limited_request(self, url: str, headers: dict = None, **kwargs) -> Optional[requests.Response]:
    #     """
    #     Make a rate-limited HTTP request
        
    #     Args:
    #         url: URL to request
    #         headers: Optional headers dictionary
    #         **kwargs: Additional arguments for requests.get()
    #     """
        
    #     # Check cache first
    #     if self.use_cache:
    #         cached = self._get_cached_response(url)
    #         if cached:
    #             logger.debug(f"Using cached response for {url}")
    #             # Create mock response
    #             from requests import Response
    #             response = Response()
    #             response._content = cached.encode()
    #             response.status_code = 200
    #             return response
        
    #     # Rate limiting
    #     current_time = time.time()
    #     time_since_last = current_time - self.last_request_time
        
    #     if time_since_last < self.rate_limit:
    #         time.sleep(self.rate_limit - time_since_last)
        
    #     self.last_request_time = time.time()
        
    #     try:
    #         # Merge headers if provided
    #         if headers:
    #             if 'headers' in kwargs:
    #                 kwargs['headers'].update(headers)
    #             else:
    #                 kwargs['headers'] = headers
            
    #         response = self.session.get(url, timeout=120, **kwargs)
            
    #         # Cache on success
    #         if response and response.status_code == 200 and self.use_cache:
    #             self._cache_response(url, response.text)
            
    #         return response
    #     except Exception as e:
    #         logger.warning(f"Request failed for {url}: {e}")
    #         return None


    # def _rate_limited_request(self, url: str, headers: dict = None, **kwargs) -> Optional[requests.Response]:
    #     """
    #     Make a rate-limited HTTP request with caching support
        
    #     Args:
    #         url: URL to request
    #         headers: Optional headers dictionary
    #         **kwargs: Additional arguments for requests.get()
    #     """
    #     with self.request_lock:  # ← ADD THIS
    #         # Rate limiting
    #         current_time = time.time()
    #         time_since_last = current_time - self.last_request_time
            
    #         if time_since_last < self.rate_limit:
    #             time.sleep(self.rate_limit - time_since_last)
            
    #         self.last_request_time = time.time()


    #     # Check cache first (only for GET requests without special headers that modify response)
    #     if self.cache.enabled and headers is None:
    #         cached_text = self.cache.get(url)
    #         if cached_text:
    #             # Create a mock Response object
    #             from requests import Response
    #             response = Response()
    #             response._content = cached_text.encode('utf-8')
    #             response.status_code = 200
    #             response.headers['Content-Type'] = 'application/json'
    #             return response
        
    #     # Rate limiting
    #     current_time = time.time()
    #     time_since_last = current_time - self.last_request_time
        
    #     if time_since_last < self.rate_limit:
    #         time.sleep(self.rate_limit - time_since_last)
        
    #     self.last_request_time = time.time()
        
    #     try:
    #         # Merge headers if provided
    #         if headers:
    #             if 'headers' in kwargs:
    #                 kwargs['headers'].update(headers)
    #             else:
    #                 kwargs['headers'] = headers
            
    #         response = self.session.get(url, timeout=150, **kwargs)
            
    #         # Cache successful responses
    #         if response.status_code == 200 and self.cache.enabled:
    #             self.cache.set(url, response.text)
            
    #         return response
    #     except Exception as e:
    #         logger.warning(f"Request failed for {url}: {e}")
    #         return None



    def _rate_limited_request(self, url: str, headers: dict = None, max_retries: int = 2, **kwargs) -> Optional[requests.Response]:
        """
        Make a rate-limited HTTP request with caching and retry support
        
        Args:
            url: URL to request
            headers: Optional headers dictionary
            max_retries: Number of retry attempts on timeout/failure
            **kwargs: Additional arguments for requests.get()
        """
        # Check cache first (only for GET requests without special headers that modify response)
        if self.cache.enabled and headers is None:
            cached_text = self.cache.get(url)
            if cached_text:
                # Create a mock Response object
                from requests import Response
                response = Response()
                response._content = cached_text.encode('utf-8')
                response.status_code = 200
                response.headers['Content-Type'] = 'application/json'
                return response
        
        # Try with retries
        for attempt in range(max_retries):
            # Rate limiting (thread-safe)
            with self.request_lock:
                current_time = time.time()
                time_since_last = current_time - self.last_request_time
                
                if time_since_last < self.rate_limit:
                    time.sleep(self.rate_limit - time_since_last)
                
                self.last_request_time = time.time()
            
            try:
                # Merge headers if provided
                if headers:
                    if 'headers' in kwargs:
                        kwargs['headers'].update(headers)
                    else:
                        kwargs['headers'] = headers
                
                # Increased timeout from 150s to 300s (5 minutes)
                response = self.session.get(url, timeout=300, **kwargs)
                
                # Cache successful responses
                if response.status_code == 200 and self.cache.enabled:
                    self.cache.set(url, response.text)
                
                return response
                
            except requests.exceptions.Timeout:
                logger.warning(f"BV-BRC timeout on {url} (attempt {attempt + 1}/{max_retries})")
                if attempt < max_retries - 1:
                    wait_time = 5 * (attempt + 1)  # 5s, 10s
                    logger.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                    continue
                logger.error(f"BV-BRC request timeout after {max_retries} attempts: {url}")
                return None
                
            except requests.exceptions.RequestException as e:
                logger.warning(f"BV-BRC request error on {url}: {e} (attempt {attempt + 1}/{max_retries})")
                if attempt < max_retries - 1:
                    time.sleep(5)
                    continue
                logger.error(f"BV-BRC request failed after {max_retries} attempts: {url}")
                return None
                
            except Exception as e:
                logger.error(f"Unexpected error for {url}: {e}")
                return None
        
        return None
    
    def _write_fasta_uppercase(self, content: str, output_file: Path):
        """
        Write FASTA content with sequences in uppercase
        
        Args:
            content: FASTA formatted text
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            for line in content.split('\n'):
                if line.startswith('>'):
                    # Keep header as-is
                    f.write(line + '\n')
                else:
                    # Convert sequence to uppercase
                    f.write(line.upper() + '\n')
    
    def download_genome(
        self,
        genome_id: str,
        output_dir: Path,
        terminal_log: DownloadLogger
    ) -> bool:
        """
        Download genome files from BV-BRC Web API
        
        Args:
            genome_id: BV-BRC genome ID (e.g., "670897.3")
            output_dir: Output directory for files
            terminal_log: Logger for this terminal
        
        Returns:
            True if successful
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        terminal_log.info(f"Downloading genome {genome_id} from BV-BRC Web API")
        logger.info(f"Downloading genome {genome_id} via Web API")
        
        # Define output files
        genome_fasta = output_dir / "genome.fasta"
        genome_faa = output_dir / "genome.faa"
        genome_ffn = output_dir / "genome.ffn"
        genome_features = output_dir / "genome.features.tab"
        
        # Check if files already exist
        if all(is_file_complete(f) for f in [genome_fasta, genome_faa, genome_ffn]):
            terminal_log.info("All genome files already exist and are complete")
            logger.info(f"Genome {genome_id} already downloaded")
            return True
        
        success_count = 0
        
        # Download genome FASTA (contigs) - JSON format, convert to FASTA
        if not is_file_complete(genome_fasta):
            terminal_log.info("Downloading genome.fasta (contigs)...")
            url = f"{self.base_url}/api/genome_sequence/?eq(genome_id,{genome_id})&limit(25000)"
            
            response = self._rate_limited_request(url)
            if response and response.status_code == 200:
                try:
                    data = response.json()  # Parse as JSON
                    if data:
                        with open(genome_fasta, 'w') as f:
                            for contig in data:
                                sequence_id = contig.get('sequence_id', contig.get('accession', 'unknown'))
                                sequence = contig.get('sequence', '')
                                description = contig.get('description', '')
                                if sequence:
                                    # Write in FASTA format with uppercase sequence
                                    f.write(f">{sequence_id} {description}\n")
                                    # Convert to uppercase and write in 60-character lines
                                    sequence_upper = sequence.upper()
                                    for i in range(0, len(sequence_upper), 60):
                                        f.write(sequence_upper[i:i+60] + "\n")
                        terminal_log.success(f"Downloaded genome.fasta ({len(data)} contigs)")
                        success_count += 1
                    else:
                        terminal_log.warning("No contig data returned")
                except Exception as e:
                    terminal_log.error(f"Error processing contig data: {e}")
            else:
                terminal_log.warning(f"Failed to download genome sequence (status: {response.status_code if response else 'timeout'})")
        else:
            success_count += 1
            terminal_log.info("genome.fasta already exists")
        
        # Download protein sequences (FAA) - Use FASTA endpoint with Accept header
        if not is_file_complete(genome_faa):
            terminal_log.info("Downloading genome.faa (proteins)...")
            url = f"{self.base_url}/api/genome_feature/?and(eq(genome_id,{genome_id}),eq(annotation,PATRIC),eq(feature_type,CDS))&limit(25000)"
            
            # Use Accept header for protein FASTA format
            headers = {'Accept': 'application/protein+fasta'}
            
            response = self._rate_limited_request(url, headers=headers)
            if response and response.status_code == 200:
                content = response.text
                if content and len(content) > 100:
                    # Write with uppercase conversion
                    self._write_fasta_uppercase(content, genome_faa)
                    
                    # Count sequences
                    seq_count = content.count('>')
                    terminal_log.success(f"Downloaded genome.faa ({seq_count} proteins)")
                    success_count += 1
                else:
                    terminal_log.warning("Protein FASTA response empty or too small")
            else:
                terminal_log.warning(f"Failed to download proteins (status: {response.status_code if response else 'timeout'})")
        else:
            success_count += 1
            terminal_log.info("genome.faa already exists")
        
        # Download CDS nucleotide sequences (FFN) - Use DNA FASTA endpoint
        if not is_file_complete(genome_ffn):
            terminal_log.info("Downloading genome.ffn (CDS nucleotides)...")
            url = f"{self.base_url}/api/genome_feature/?and(eq(genome_id,{genome_id}),eq(annotation,PATRIC),eq(feature_type,CDS))&limit(25000)"
            
            # Use Accept header for DNA FASTA format
            headers = {'Accept': 'application/dna+fasta'}
            
            response = self._rate_limited_request(url, headers=headers)
            if response and response.status_code == 200:
                content = response.text
                if content and len(content) > 100:
                    # Write with uppercase conversion
                    self._write_fasta_uppercase(content, genome_ffn)
                    
                    # Count sequences
                    seq_count = content.count('>')
                    terminal_log.success(f"Downloaded genome.ffn ({seq_count} CDS)")
                    success_count += 1
                else:
                    terminal_log.warning("CDS FASTA response empty or too small")
            else:
                terminal_log.warning(f"Failed to download CDS (status: {response.status_code if response else 'timeout'})")
        else:
            success_count += 1
            terminal_log.info("genome.ffn already exists")
        
        # Download features table (optional but useful)
        if not is_file_complete(genome_features):
            terminal_log.info("Downloading genome.features.tab (annotations)...")
            url = f"{self.base_url}/api/genome_feature/?and(eq(genome_id,{genome_id}),eq(annotation,PATRIC))&limit(25000)"
            
            response = self._rate_limited_request(url)
            if response and response.status_code == 200:
                try:
                    data = response.json()
                    if data:
                        # Write TSV header
                        with open(genome_features, 'w') as f:
                            # Get all keys from first feature
                            if data:
                                keys = list(data[0].keys())
                                f.write('\t'.join(keys) + '\n')
                                
                                # Write data
                                for feature in data:
                                    values = [str(feature.get(k, '')) for k in keys]
                                    f.write('\t'.join(values) + '\n')
                        terminal_log.success(f"Downloaded genome.features.tab ({len(data)} features)")
                    else:
                        terminal_log.warning("No feature data returned")
                except Exception as e:
                    terminal_log.error(f"Error processing features: {e}")
            else:
                terminal_log.warning(f"Failed to download features (status: {response.status_code if response else 'timeout'})")
        
        # Success if we got at least the core files
        if success_count >= 3:
            terminal_log.success(f"Successfully downloaded genome {genome_id} ({success_count} files)")
            logger.info(f"Successfully downloaded genome {genome_id}")
            return True
        else:
            terminal_log.error(f"Failed to download genome {genome_id} (only {success_count}/3 core files)")
            logger.error(f"Incomplete download for genome {genome_id}")
            return False


# class BVBRCDownloader:
#     """Download genome files from BV-BRC Web API"""
    
#     def __init__(self, rate_limit: float = 0.5):
#         """
#         Initialize BV-BRC downloader
        
#         Args:
#             rate_limit: Minimum seconds between API requests
#         """
#         self.rate_limit = rate_limit
#         self.last_request_time = 0
#         self.base_url = "https://www.bv-brc.org"


#     def _rate_limited_request(self, url: str, headers: dict = None, **kwargs) -> Optional[requests.Response]:
#         """
#         Make a rate-limited HTTP request
        
#         Args:
#             url: URL to request
#             headers: Optional headers dictionary
#             **kwargs: Additional arguments for requests.get()
#         """
#         # Rate limiting
#         current_time = time.time()
#         time_since_last = current_time - self.last_request_time
        
#         if time_since_last < self.rate_limit:
#             time.sleep(self.rate_limit - time_since_last)
        
#         self.last_request_time = time.time()
        
#         try:
#             # Merge headers if provided
#             if headers:
#                 if 'headers' in kwargs:
#                     kwargs['headers'].update(headers)
#                 else:
#                     kwargs['headers'] = headers
            
#             response = requests.get(url, timeout=30, **kwargs)
#             return response
#         except Exception as e:
#             logger.warning(f"Request failed for {url}: {e}")
#             return None


#     def download_genome(self, genome_id: str, output_dir: Path, terminal_log: DownloadLogger) -> bool:

#         """
#         Download genome files from BV-BRC Web API
        
#         Args:
#             genome_id: BV-BRC genome ID (e.g., "670897.3")
#             output_dir: Output directory for files
#             terminal_log: Logger for this terminal
        
#         Returns:
#             True if successful
#         """
#         output_dir = Path(output_dir)
#         output_dir.mkdir(parents=True, exist_ok=True)
        
#         terminal_log.info(f"Downloading genome {genome_id} from BV-BRC Web API")
#         logger.info(f"Downloading genome {genome_id} via Web API")
        
#         # Define output files
#         genome_fasta = output_dir / "genome.fasta"
#         genome_faa = output_dir / "genome.faa"
#         genome_ffn = output_dir / "genome.ffn"
#         genome_features = output_dir / "genome.features.tab"
        
#         # Check if files already exist
#         if all(is_file_complete(f) for f in [genome_fasta, genome_faa, genome_ffn]):
#             terminal_log.info("All genome files already exist and are complete")
#             logger.info(f"Genome {genome_id} already downloaded")
#             return True
        
#         success_count = 0
        
#         # Download genome FASTA (contigs)
#         if not is_file_complete(genome_fasta):
#             terminal_log.info("Downloading genome.fasta (contigs)...")
#             url = f"{self.base_url}/api/genome_sequence/?eq(genome_id,{genome_id})&limit(25000)"
            
#             response = self._rate_limited_request(url)
#             if response and response.status_code == 200:
#                 try:
#                     data = response.json()  # Parse as JSON instead of text
#                     if data:
#                         with open(genome_fasta, 'w') as f:
#                             for contig in data:
#                                 sequence_id = contig.get('sequence_id', contig.get('accession', 'unknown'))
#                                 sequence = contig.get('sequence', '')
#                                 description = contig.get('description', '')
#                                 if sequence:
#                                     # Write in FASTA format
#                                     f.write(f">{sequence_id} {description}\n")
#                                     # Write sequence in 60-character lines
#                                     for i in range(0, len(sequence), 60):
#                                         f.write(sequence[i:i+60] + "\n")
#                         terminal_log.success(f"Downloaded genome.fasta ({len(data)} contigs)")
#                         success_count += 1
#                     else:
#                         terminal_log.warning("No contig data returned")
#                 except Exception as e:
#                     terminal_log.error(f"Error processing contig data: {e}")
#             else:
#                 terminal_log.warning(f"Failed to download genome sequence (status: {response.status_code if response else 'timeout'})")
#         else:
#             success_count += 1
#             terminal_log.info("genome.fasta already exists")
                
#         # Download protein sequences (FAA)
#         if not is_file_complete(genome_faa):
#             terminal_log.info("Downloading genome.faa (proteins)...")
#             # Use the feature endpoint with FASTA accept header
#             url = f"{self.base_url}/api/genome_feature/?and(eq(genome_id,{genome_id}),eq(annotation,PATRIC),eq(feature_type,CDS))&limit(25000)"
            
#             # Add Accept header for protein FASTA format
#             headers = {'Accept': 'application/protein+fasta'}
            
#             response = self._rate_limited_request(url, headers=headers)
#             if response and response.status_code == 200:
#                 content = response.text
#                 if content and len(content) > 100:
#                     with open(genome_faa, 'w') as f:
#                         f.write(content)
                    
#                     # Count sequences
#                     seq_count = content.count('>')
#                     terminal_log.success(f"Downloaded genome.faa ({seq_count} proteins)")
#                     success_count += 1
#                 else:
#                     terminal_log.warning("Protein FASTA response empty or too small")
#             else:
#                 terminal_log.warning(f"Failed to download proteins")

        
#         # Download CDS nucleotide sequences (FFN)
#         if not is_file_complete(genome_ffn):
#             terminal_log.info("Downloading genome.ffn (CDS nucleotides)...")
#             url = f"{self.base_url}/api/genome_feature/?and(eq(genome_id,{genome_id}),eq(annotation,PATRIC),eq(feature_type,CDS))&limit(25000)"
            
#             # Add Accept header for DNA FASTA format
#             headers = {'Accept': 'application/dna+fasta'}
            
#             response = self._rate_limited_request(url, headers=headers)
#             if response and response.status_code == 200:
#                 content = response.text
#                 if content and len(content) > 100:
#                     with open(genome_ffn, 'w') as f:
#                         f.write(content)
                    
#                     seq_count = content.count('>')
#                     terminal_log.success(f"Downloaded genome.ffn ({seq_count} CDS)")
#                     success_count += 1
#                 else:
#                     terminal_log.warning("CDS FASTA response empty or too small")
#             else:
#                 terminal_log.warning(f"Failed to download CDS")
        
#         # Download features table (optional but useful)
#         if not is_file_complete(genome_features):
#             terminal_log.info("Downloading genome.features.tab (annotations)...")
#             url = f"{self.base_url}/api/genome_feature/?and(eq(genome_id,{genome_id}),eq(annotation,PATRIC))&limit(25000)"
            
#             response = self._rate_limited_request(url)
#             if response and response.status_code == 200:
#                 try:
#                     data = response.json()
#                     if data:
#                         # Write TSV header
#                         with open(genome_features, 'w') as f:
#                             # Get all keys from first feature
#                             if data:
#                                 keys = list(data[0].keys())
#                                 f.write('\t'.join(keys) + '\n')
                                
#                                 # Write data
#                                 for feature in data:
#                                     values = [str(feature.get(k, '')) for k in keys]
#                                     f.write('\t'.join(values) + '\n')
#                         terminal_log.success(f"Downloaded genome.features.tab ({len(data)} features)")
#                     else:
#                         terminal_log.warning("No feature data returned")
#                 except Exception as e:
#                     terminal_log.error(f"Error processing features: {e}")
#             else:
#                 terminal_log.warning(f"Failed to download features (status: {response.status_code if response else 'timeout'})")
        
#         # Success if we got at least the core files
#         if success_count >= 3:
#             terminal_log.success(f"Successfully downloaded genome {genome_id} ({success_count} files)")
#             logger.info(f"Successfully downloaded genome {genome_id}")
#             return True
#         else:
#             terminal_log.error(f"Failed to download genome {genome_id} (only {success_count}/3 core files)")
#             logger.error(f"Incomplete download for genome {genome_id}")
#             return False