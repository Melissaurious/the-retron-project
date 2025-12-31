"""
NCBI downloader for bacteria and archaea genomes.

Downloads directly from NCBI FTP using assembly accessions (GCA_/GCF_).
This is a simpler, cleaner version that follows the BaseDownloader interface.
"""

import logging
import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin
import requests
import shutil

from .base_downloader import BaseDownloader

logger = logging.getLogger(__name__)


class NCBIDownloaderSimple(BaseDownloader):
    """
    Simple NCBI downloader for direct FTP downloads.
    
    Uses assembly accessions (GCA_/GCF_) to download from NCBI FTP.
    For more complex operations (E-utilities, protein downloads),
    use the full NCBIDownloader class.
    """
    
    # NCBI FTP base URLs
    GENBANK_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
    REFSEQ_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
    
    def __init__(self, database_name: str = "ncbi"):
        """
        Initialize NCBI downloader.
        
        Args:
            database_name: Database name (e.g., 'ncbi_bacteria', 'ncbi_archaea')
        """
        super().__init__(database_name)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'GenomeFetcher/1.0'
        })
    
    def download_genome(
        self,
        genome_record: Dict,
        output_dir: Path,
        required_files: List[str] = None,
        compress: bool = True,
        terminal_log: 'DownloadLogger' = None
    ) -> bool:
        """
        Download genome from NCBI FTP.
        
        Args:
            genome_record: Dictionary with NCBI metadata
            output_dir: Base output directory
            required_files: List of files to download (default: ['genome.fasta'])
            compress: Whether to compress output (default: True)
            terminal_log: Logger for progress
        
        Returns:
            True if successful
        """
        if required_files is None:
            required_files = ['genome.fasta']
        
        self.update_stats('attempted')
        
        # Extract assembly accession
        assembly_acc = self.get_genome_id(genome_record)
        
        if terminal_log:
            terminal_log.info(f"Processing NCBI genome: {assembly_acc}")
        
        # Validate accession
        if not self.validate_genome_id(assembly_acc):
            if terminal_log:
                terminal_log.error(f"Invalid NCBI accession format: {assembly_acc}")
            self.update_stats('failed')
            return False
        
        # Create genome directory
        genome_dir = output_dir / assembly_acc
        genome_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if already complete
        if self.is_genome_complete(genome_dir, required_files):
            if terminal_log:
                terminal_log.info(f"Genome already complete, skipping")
            self.update_stats('skipped')
            return True
        
        # Get FTP path from metadata (if available)
        ftp_path = genome_record.get('ftp_path')
        
        if ftp_path:
            # Use provided FTP path
            success = self._download_from_ftp_path(
                ftp_path,
                genome_dir,
                assembly_acc,
                compress,
                terminal_log
            )
        else:
            # Construct FTP path from accession
            success = self._download_from_accession(
                assembly_acc,
                genome_dir,
                compress,
                terminal_log
            )
        
        if success:
            self.update_stats('successful')
            if terminal_log:
                terminal_log.success(f"Successfully downloaded {assembly_acc}")
        else:
            self.update_stats('failed')
            if terminal_log:
                terminal_log.error(f"Failed to download {assembly_acc}")
        
        return success
    
    def _download_from_ftp_path(
        self,
        ftp_path: str,
        genome_dir: Path,
        assembly_acc: str,
        compress: bool,
        terminal_log: 'DownloadLogger'
    ) -> bool:
        """
        Download genome using provided FTP path.
        
        Args:
            ftp_path: Full FTP path to assembly directory
            genome_dir: Output directory
            assembly_acc: Assembly accession
            compress: Whether to compress
            terminal_log: Logger
        
        Returns:
            True if successful
        """
        try:
            # Extract assembly name from path
            # Example: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/600/855/GCA_036600855.1_ASM3660085v1
            assembly_name = ftp_path.rstrip('/').split('/')[-1]
            
            if terminal_log:
                terminal_log.info(f"Assembly: {assembly_name}")
            
            # Build FASTA URL
            fasta_filename = f"{assembly_name}_genomic.fna.gz"
            fasta_url = urljoin(ftp_path + '/', fasta_filename)
            
            if terminal_log:
                terminal_log.info(f"Downloading: {fasta_url}")
            
            return self._download_fasta(fasta_url, genome_dir, compress, terminal_log)
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"Download error: {e}")
            logger.error(f"FTP path download error for {assembly_acc}: {e}")
            return False
    
    def _download_from_accession(
        self,
        assembly_acc: str,
        genome_dir: Path,
        compress: bool,
        terminal_log: 'DownloadLogger'
    ) -> bool:
        """
        Download genome by constructing FTP URL from accession.
        
        Args:
            assembly_acc: Assembly accession
            genome_dir: Output directory
            compress: Whether to compress
            terminal_log: Logger
        
        Returns:
            True if successful
        """
        try:
            # Build FTP directory URL
            ftp_dir_url = self._build_ftp_url(assembly_acc)
            if not ftp_dir_url:
                return False
            
            if terminal_log:
                terminal_log.info(f"FTP directory: {ftp_dir_url}")
            
            # Get assembly name from directory listing
            assembly_name = self._get_assembly_name(ftp_dir_url, assembly_acc, terminal_log)
            if not assembly_name:
                return False
            
            # Build FASTA URL
            fasta_filename = f"{assembly_name}_genomic.fna.gz"
            fasta_url = urljoin(ftp_dir_url + assembly_name + '/', fasta_filename)
            
            if terminal_log:
                terminal_log.info(f"Downloading: {fasta_url}")
            
            return self._download_fasta(fasta_url, genome_dir, compress, terminal_log)
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"Download error: {e}")
            logger.error(f"NCBI download error for {assembly_acc}: {e}")
            return False
    
    def _build_ftp_url(self, assembly_acc: str) -> Optional[str]:
        """
        Build NCBI FTP URL from accession.
        
        Args:
            assembly_acc: Assembly accession (e.g., 'GCA_036600855.1')
        
        Returns:
            FTP base URL or None
        """
        match = re.match(r'(GCA|GCF)_(\d{3})(\d{3})(\d{3})\.(\d+)', assembly_acc)
        if not match:
            logger.error(f"Invalid NCBI accession format: {assembly_acc}")
            return None
        
        prefix, part1, part2, part3, version = match.groups()
        base_url = self.GENBANK_FTP if prefix == 'GCA' else self.REFSEQ_FTP
        path = f"{part1}/{part2}/{part3}/"
        
        return base_url + path
    
    def _get_assembly_name(self, ftp_dir_url: str, assembly_acc: str, terminal_log) -> Optional[str]:
        """
        Get full assembly directory name from FTP listing.
        
        Args:
            ftp_dir_url: Base FTP directory URL
            assembly_acc: Assembly accession
            terminal_log: Logger
        
        Returns:
            Full assembly name or None
        """
        try:
            response = self.session.get(ftp_dir_url, timeout=30)
            if response.status_code != 200:
                if terminal_log:
                    terminal_log.error(f"Cannot access FTP: HTTP {response.status_code}")
                return None
            
            # Parse directory listing
            pattern = f'{assembly_acc}_[^/"]+/'
            matches = re.findall(pattern, response.text)
            
            if matches:
                assembly_name = matches[0].rstrip('/')
                if terminal_log:
                    terminal_log.info(f"Found assembly: {assembly_name}")
                return assembly_name
            
            if terminal_log:
                terminal_log.error("Assembly not found in FTP listing")
            return None
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"FTP listing error: {e}")
            return None
    
    def _download_fasta(
        self,
        fasta_url: str,
        genome_dir: Path,
        compress: bool,
        terminal_log: 'DownloadLogger'
    ) -> bool:
        """
        Download FASTA file with appropriate compression handling.
        
        Args:
            fasta_url: URL to FASTA file
            genome_dir: Output directory
            compress: Whether to compress output
            terminal_log: Logger
        
        Returns:
            True if successful
        """
        try:
            response = self.session.get(fasta_url, timeout=300, stream=True)
            if response.status_code != 200:
                if terminal_log:
                    terminal_log.error(f"Download failed: HTTP {response.status_code}")
                return False
            
            # Determine output file
            if compress:
                output_file = genome_dir / "genome.fasta.gz"
            else:
                output_file = genome_dir / "genome.fasta"
            
            # Download with appropriate compression
            if compress and fasta_url.endswith('.gz'):
                self._download_direct(response, output_file)
            elif not compress and fasta_url.endswith('.gz'):
                self._download_and_decompress(response, output_file)
            else:
                if compress:
                    self._download_and_compress(response, output_file)
                else:
                    self._download_direct(response, output_file)
            
            file_size = output_file.stat().st_size
            if terminal_log:
                terminal_log.success(f"Downloaded genome.fasta ({file_size:,} bytes)")
            
            return True
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"FASTA download error: {e}")
            return False
    
    def _download_direct(self, response, output_file: Path):
        """Download file directly."""
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    
    def _download_and_compress(self, response, output_file: Path):
        """Download and compress on the fly."""
        with gzip.open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    
    def _download_and_decompress(self, response, output_file: Path):
        """Download and decompress on the fly."""
        with gzip.GzipFile(fileobj=response.raw) as gz_file:
            with open(output_file, 'wb') as f:
                shutil.copyfileobj(gz_file, f)
    
    def get_metadata(self, genome_record: Dict) -> Dict:
        """
        Extract metadata from NCBI record.
        
        Args:
            genome_record: NCBI genome record
        
        Returns:
            Standardized metadata dictionary
        """
        assembly_acc = self.get_genome_id(genome_record)
        
        metadata = {
            'genome_id': assembly_acc,
            'database': self.database_name,
            'assembly_accession': assembly_acc,
            'bioproject': genome_record.get('bioproject'),
            'biosample': genome_record.get('biosample'),
            'organism_name': genome_record.get('organism_name'),
            'assembly_level': genome_record.get('assembly_level'),
            'genome_rep': genome_record.get('genome_rep'),
            'ftp_path': genome_record.get('ftp_path')
        }
        
        return metadata
    
    def validate_genome_id(self, genome_id: str) -> bool:
        """
        Validate NCBI assembly accession format.
        
        Args:
            genome_id: Assembly accession
        
        Returns:
            True if valid format
        
        NCBI format: GCA_XXXXXXXXX.X or GCF_XXXXXXXXX.X
        Example: GCA_036600855.1, GCF_000005845.2
        """
        pattern = r'^(GCA|GCF)_\d{9}\.\d+$'
        return bool(re.match(pattern, genome_id))