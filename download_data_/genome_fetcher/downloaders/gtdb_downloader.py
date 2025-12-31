"""
GTDB downloader for bacteria and archaea genomes.

GTDB uses genome accessions with prefixes:
- GB_ = GenBank assembly (RefSeq not available)
- RS_ = RefSeq assembly (preferred)

These map to NCBI accessions (GCA_/GCF_) which we download from NCBI FTP.
"""

import logging
import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin
import requests

from .base_downloader import BaseDownloader

logger = logging.getLogger(__name__)


class GTDBDownloader(BaseDownloader):
    """
    Downloader for GTDB (Genome Taxonomy Database) genomes.
    
    GTDB provides curated bacterial and archaeal genomes with
    standardized taxonomy. Genomes are sourced from NCBI but
    use GTDB-specific identifiers.
    """
    
    # NCBI FTP base URLs
    GENBANK_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
    REFSEQ_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
    
    def __init__(self, database_name: str = "gtdb"):
        """
        Initialize GTDB downloader.
        
        Args:
            database_name: Database name (e.g., 'gtdb_bacteria', 'gtdb_archaea')
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
        Download genome from NCBI via GTDB accession.
        
        Args:
            genome_record: Dictionary with GTDB metadata
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
        gtdb_id = self.get_genome_id(genome_record)
        
        if terminal_log:
            terminal_log.info(f"Processing GTDB genome: {gtdb_id}")
        
        # Validate and convert GTDB ID to NCBI accession
        ncbi_accession = self._gtdb_to_ncbi(gtdb_id)
        if not ncbi_accession:
            if terminal_log:
                terminal_log.error(f"Invalid GTDB ID format: {gtdb_id}")
            self.update_stats('failed')
            return False
        
        if terminal_log:
            terminal_log.info(f"NCBI accession: {ncbi_accession}")
        
        # Create genome directory
        genome_dir = output_dir / gtdb_id
        genome_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if already complete
        if self.is_genome_complete(genome_dir, required_files):
            if terminal_log:
                terminal_log.info(f"Genome already complete, skipping")
            self.update_stats('skipped')
            return True
        
        # Download from NCBI FTP
        success = self._download_from_ncbi_ftp(
            ncbi_accession,
            genome_dir,
            gtdb_id,
            required_files,
            compress,
            terminal_log
        )
        
        if success:
            self.update_stats('successful')
            if terminal_log:
                terminal_log.success(f"Successfully downloaded {gtdb_id}")
        else:
            self.update_stats('failed')
            if terminal_log:
                terminal_log.error(f"Failed to download {gtdb_id}")
        
        return success
    
    def _gtdb_to_ncbi(self, gtdb_id: str) -> Optional[str]:
        """
        Convert GTDB ID to NCBI accession.
        
        Args:
            gtdb_id: GTDB identifier (e.g., 'GB_GCA_041569275.1', 'RS_GCF_034719275.1')
        
        Returns:
            NCBI accession (e.g., 'GCA_041569275.1') or None if invalid
        """
        # Remove GB_ or RS_ prefix
        if gtdb_id.startswith('GB_'):
            return gtdb_id[3:]  # Remove 'GB_'
        elif gtdb_id.startswith('RS_'):
            return gtdb_id[3:]  # Remove 'RS_'
        else:
            logger.warning(f"GTDB ID missing prefix: {gtdb_id}")
            # Might already be an NCBI accession
            if gtdb_id.startswith('GCA_') or gtdb_id.startswith('GCF_'):
                return gtdb_id
            return None
    
    def _build_ncbi_ftp_url(self, ncbi_accession: str) -> Optional[str]:
        """
        Build NCBI FTP URL from accession.
        
        NCBI FTP structure:
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/041/569/275/GCA_041569275.1_ASM4156927v1/
        
        Args:
            ncbi_accession: NCBI accession (e.g., 'GCA_041569275.1')
        
        Returns:
            FTP base URL or None if invalid format
        """
        # Parse accession: GCA_041569275.1
        match = re.match(r'(GCA|GCF)_(\d{3})(\d{3})(\d{3})\.(\d+)', ncbi_accession)
        if not match:
            logger.error(f"Invalid NCBI accession format: {ncbi_accession}")
            return None
        
        prefix, part1, part2, part3, version = match.groups()
        
        # Choose base URL
        base_url = self.GENBANK_FTP if prefix == 'GCA' else self.REFSEQ_FTP
        
        # Build path: GCA/041/569/275/
        path = f"{part1}/{part2}/{part3}/"
        
        return base_url + path
    
    def _get_assembly_name(self, ftp_dir_url: str, ncbi_accession: str, terminal_log) -> Optional[str]:
        """
        Get the full assembly directory name from NCBI FTP.
        
        NCBI names directories like: GCA_041569275.1_ASM4156927v1
        We need to find this exact name.
        
        Args:
            ftp_dir_url: Base FTP directory URL
            ncbi_accession: NCBI accession
            terminal_log: Logger
        
        Returns:
            Full assembly directory name or None
        """
        try:
            # List directory
            response = self.session.get(ftp_dir_url, timeout=30)
            if response.status_code != 200:
                if terminal_log:
                    terminal_log.error(f"Cannot access FTP directory: {response.status_code}")
                return None
            
            # Parse HTML directory listing (handle both text and bytes)
            if isinstance(response.content, bytes):
                html = response.content.decode('utf-8', errors='ignore')
            else:
                html = response.text
            
            pattern = f'{ncbi_accession}_[^/"]+/'
            matches = re.findall(pattern, html)
            
            if matches:
                # Remove trailing slash
                assembly_name = matches[0].rstrip('/')
                if terminal_log:
                    terminal_log.info(f"Found assembly: {assembly_name}")
                return assembly_name
            
            if terminal_log:
                terminal_log.error(f"Assembly directory not found in FTP listing")
            return None
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"Error accessing FTP: {e}")
            logger.error(f"FTP listing error for {ncbi_accession}: {e}")
            return None
    
    def _download_from_ncbi_ftp(
        self,
        ncbi_accession: str,
        genome_dir: Path,
        gtdb_id: str,
        required_files: List[str],
        compress: bool,
        terminal_log: 'DownloadLogger'
    ) -> bool:
        """
        Download genome files from NCBI FTP.
        
        Args:
            ncbi_accession: NCBI accession
            genome_dir: Output directory
            gtdb_id: GTDB identifier (for logging)
            required_files: Files to download
            compress: Whether to compress
            terminal_log: Logger
        
        Returns:
            True if successful
        """
        try:
            # Build FTP directory URL
            ftp_dir_url = self._build_ncbi_ftp_url(ncbi_accession)
            if not ftp_dir_url:
                return False
            
            if terminal_log:
                terminal_log.info(f"FTP directory: {ftp_dir_url}")
            
            # Get full assembly name
            assembly_name = self._get_assembly_name(ftp_dir_url, ncbi_accession, terminal_log)
            if not assembly_name:
                return False
            
            # Build full FTP URL for genomic FASTA
            # File format: GCA_041569275.1_ASM4156927v1_genomic.fna.gz
            fasta_filename = f"{assembly_name}_genomic.fna.gz"
            fasta_url = urljoin(ftp_dir_url + assembly_name + '/', fasta_filename)
            
            if terminal_log:
                terminal_log.info(f"Downloading: {fasta_url}")
            
            # Download FASTA
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
            
            # Download with appropriate compression handling
            if compress and fasta_url.endswith('.gz'):
                # Already compressed, download directly
                self._download_direct(response, output_file, terminal_log)
            elif not compress and fasta_url.endswith('.gz'):
                # Need to decompress
                self._download_and_decompress(response, output_file, terminal_log)
            else:
                # Rare case: uncompressed source
                if compress:
                    self._download_and_compress(response, output_file, terminal_log)
                else:
                    self._download_direct(response, output_file, terminal_log)
            
            file_size = output_file.stat().st_size
            if terminal_log:
                terminal_log.success(f"Downloaded genome.fasta ({file_size:,} bytes)")
            
            return True
        
        except Exception as e:
            if terminal_log:
                terminal_log.error(f"Download error: {e}")
            logger.error(f"NCBI FTP download error for {gtdb_id}: {e}")
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
        import shutil
        import io
        
        # Decompress to memory buffer first
        with gzip.GzipFile(fileobj=response.raw) as gz_file:
            content = gz_file.read().decode('utf-8')
        
        # Write with uppercase conversion
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
        Extract metadata from GTDB record.
        
        Args:
            genome_record: GTDB genome record
        
        Returns:
            Standardized metadata dictionary
        """
        gtdb_id = self.get_genome_id(genome_record)
        
        metadata = {
            'genome_id': gtdb_id,
            'database': self.database_name,
            'ncbi_accession': self._gtdb_to_ncbi(gtdb_id),
            'taxonomy': genome_record.get('gtdb_taxonomy'),
            'completeness': genome_record.get('checkm_completeness'),
            'contamination': genome_record.get('checkm_contamination'),
            'genome_size': genome_record.get('genome_size'),
            'protein_count': genome_record.get('protein_count'),
            'contig_count': genome_record.get('contig_count'),
            'n50_contigs': genome_record.get('n50_contigs'),
            'scaffold_count': genome_record.get('scaffold_count')
        }
        
        return metadata
    
    def validate_genome_id(self, genome_id: str) -> bool:
        """
        Validate GTDB genome ID format.
        
        Args:
            genome_id: Genome identifier
        
        Returns:
            True if valid format
        
        GTDB format: GB_GCA_XXXXXXXXX.X or RS_GCF_XXXXXXXXX.X
        Example: GB_GCA_041569275.1, RS_GCF_034719275.1
        """
        pattern = r'^(GB|RS)_(GCA|GCF)_\d{9}\.\d+$'
        return bool(re.match(pattern, genome_id))