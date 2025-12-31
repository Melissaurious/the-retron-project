"""
MGnify downloader for human gut, marine, and soil catalogues.

MGnify provides direct FTP download URLs in their metadata, making this
one of the simplest downloaders to implement.
"""

import logging
import gzip
import shutil
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urlparse
import requests

from .base_downloader import BaseDownloader

logger = logging.getLogger(__name__)


class MGnifyDownloader(BaseDownloader):
    """
    Downloader for MGnify genome catalogues.
    
    Supports:
    - Human Gut v2.0
    - Marine v2.0
    - Soil v1.0
    - UHGG (Unified Human Gastrointestinal Genome) v2.0
    
    These catalogues provide direct FTP URLs in the metadata.
    """
    
    def __init__(self, database_name: str = "mgnify"):
        """
        Initialize MGnify downloader.
        
        Args:
            database_name: Database name (e.g., 'mgnify_human_gut', 'mgnify_marine')
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
        Download genome from MGnify FTP.
        
        Args:
            genome_record: Dictionary with MGnify metadata (must have 'FTP_download')
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
        
        # Extract genome ID
        genome_id = self.get_genome_id(genome_record)
        
        if terminal_log:
            terminal_log.info(f"Processing MGnify genome: {genome_id}")
        
        # Validate genome ID
        if not self.validate_genome_id(genome_id):
            if terminal_log:
                terminal_log.error(f"Invalid MGnify genome ID format: {genome_id}")
            self.update_stats('failed')
            return False
        
        # Create genome directory
        genome_dir = output_dir / genome_id
        genome_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if already complete
        if self.is_genome_complete(genome_dir, required_files):
            if terminal_log:
                terminal_log.info(f"Genome already complete, skipping")
            self.update_stats('skipped')
            return True
        
        # Get FTP URL
        ftp_url = genome_record.get('FTP_download')
        if not ftp_url:
            if terminal_log:
                terminal_log.error("No FTP_download URL in metadata")
            self.update_stats('failed')
            return False
        
        if terminal_log:
            terminal_log.info(f"FTP URL: {ftp_url}")
        
        # Download genome.fasta from FTP
        # MGnify provides GFF files by default, we need to get the FASTA
        success = self._download_from_ftp(
            ftp_url,
            genome_dir,
            genome_id,
            required_files,
            compress,
            terminal_log
        )
        
        if success:
            self.update_stats('successful')
            if terminal_log:
                terminal_log.success(f"Successfully downloaded {genome_id}")
        else:
            self.update_stats('failed')
            if terminal_log:
                terminal_log.error(f"Failed to download {genome_id}")
        
        return success
    
    def _download_from_ftp(
        self,
        ftp_url: str,
        genome_dir: Path,
        genome_id: str,
        required_files: List[str],
        compress: bool,
        terminal_log: 'DownloadLogger'
    ) -> bool:
        """
        Download files from MGnify FTP.
        
        MGnify FTP structure:
        - MGYG000000001.gff.gz (annotation)
        - We need to derive the FASTA URL from GFF URL
        
        Args:
            ftp_url: FTP URL from metadata
            genome_dir: Output directory
            genome_id: Genome identifier
            required_files: Files to download
            compress: Whether to compress
            terminal_log: Logger
        
        Returns:
            True if successful
        """
        try:
            # MGnify typically provides .gff.gz URLs
            # We need to derive the FASTA URL
            # Pattern: Replace .gff.gz with .fna.gz or .fasta.gz
            
            base_url = ftp_url.replace('.gff.gz', '')
            
            # Try different FASTA extensions
            fasta_extensions = ['.fna.gz', '.fasta.gz', '.fa.gz']
            
            for ext in fasta_extensions:
                fasta_url = base_url + ext
                
                if terminal_log:
                    terminal_log.info(f"Trying: {fasta_url}")
                
                try:
                    response = self.session.get(fasta_url, timeout=120, stream=True)
                    
                    if response.status_code == 200:
                        # Determine output filename
                        if compress:
                            output_file = genome_dir / "genome.fasta.gz"
                        else:
                            output_file = genome_dir / "genome.fasta"
                        
                        # Download with streaming
                        if compress and not fasta_url.endswith('.gz'):
                            # Need to compress
                            self._download_and_compress(response, output_file, terminal_log)
                        elif not compress and fasta_url.endswith('.gz'):
                            # Need to decompress
                            self._download_and_decompress(response, output_file, terminal_log)
                        else:
                            # Direct download
                            self._download_direct(response, output_file, terminal_log)
                        
                        file_size = output_file.stat().st_size
                        if terminal_log:
                            terminal_log.success(f"Downloaded genome.fasta ({file_size:,} bytes)")
                        
                        return True
                
                except Exception as e:
                    if terminal_log:
                        terminal_log.warning(f"Failed with {ext}: {e}")
                    continue
            
            # If we get here, none of the extensions worked
            if terminal_log:
                terminal_log.error("Could not find FASTA file at FTP location")
            return False
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"FTP download error: {e}")
            logger.error(f"MGnify FTP error for {genome_id}: {e}")
            return False
    
    def _download_direct(self, response, output_file: Path, terminal_log):
        """Download file directly with uppercase conversion."""
        with open(output_file, 'w', encoding='utf-8') as f:
            for line in response.iter_lines():
                # Handle both bytes and str
                if isinstance(line, bytes):
                    line = line.decode('utf-8', errors='ignore')
                if not line:
                    continue
                if line.startswith('>'):
                    f.write(line + '\n')
                else:
                    f.write(line.upper() + '\n')
    
    def _download_and_compress(self, response, output_file: Path, terminal_log):
        """Download and compress with uppercase conversion."""
        with gzip.open(output_file, 'wt', encoding='utf-8') as f:
            for line in response.iter_lines():
                # Handle both bytes and str
                if isinstance(line, bytes):
                    line = line.decode('utf-8', errors='ignore')
                if not line:
                    continue
                if line.startswith('>'):
                    f.write(line + '\n')
                else:
                    f.write(line.upper() + '\n')
    
    def _download_and_decompress(self, response, output_file: Path, terminal_log):
        """Download, decompress, and convert to uppercase."""
        import io
        
        # Decompress and convert
        with gzip.GzipFile(fileobj=response.raw) as gz_file:
            content = gz_file.read().decode('utf-8')
        
        # Write with uppercase
        with open(output_file, 'w', encoding='utf-8') as f:
            for line in content.split('\n'):
                if not line.strip():
                    continue
                if line.startswith('>'):
                    f.write(line + '\n')
                else:
                    f.write(line.upper() + '\n')
    
    def get_metadata(self, genome_record: Dict) -> Dict:
        """
        Extract metadata from MGnify record.
        
        Args:
            genome_record: MGnify genome record
        
        Returns:
            Standardized metadata dictionary
        """
        genome_id = self.get_genome_id(genome_record)
        
        metadata = {
            'genome_id': genome_id,
            'database': self.database_name,
            'catalogue': genome_record.get('catalogue', 'unknown'),
            'genome_type': genome_record.get('Genome_type'),
            'length': genome_record.get('Length'),
            'n_contigs': genome_record.get('N_contigs'),
            'n50': genome_record.get('N50'),
            'gc_content': genome_record.get('GC_content'),
            'completeness': genome_record.get('Completeness'),
            'contamination': genome_record.get('Contamination'),
            'taxonomy': genome_record.get('Lineage'),
            'country': genome_record.get('Country'),
            'continent': genome_record.get('Continent'),
            'ftp_url': genome_record.get('FTP_download')
        }
        
        return metadata
    
    def validate_genome_id(self, genome_id: str) -> bool:
        """
        Validate MGnify genome ID format.
        
        Args:
            genome_id: Genome identifier
        
        Returns:
            True if valid format
        
        MGnify format: MGYG followed by 9 digits
        Example: MGYG000000001, MGYG000296000
        """
        if not genome_id.startswith('MGYG'):
            return False
        
        # Should be MGYG + 9 digits
        if len(genome_id) != 13:
            return False
        
        # Check digits
        try:
            int(genome_id[4:])
            return True
        except ValueError:
            return False