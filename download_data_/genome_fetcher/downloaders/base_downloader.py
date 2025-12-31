"""
Abstract base class for genome downloaders.

All database-specific downloaders (GTDB, MGnify, NCBI, etc.) inherit from this class
and implement the required methods. This ensures a consistent interface across all
download sources.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Optional, List, Tuple
import logging

logger = logging.getLogger(__name__)


class BaseDownloader(ABC):
    """
    Abstract base class for genome downloaders.
    
    Each database-specific downloader must implement:
    1. download_genome() - Download genome files
    2. get_metadata() - Extract relevant metadata
    3. validate_genome_id() - Check if genome ID is valid for this database
    """
    
    def __init__(self, database_name: str):
        """
        Initialize base downloader.
        
        Args:
            database_name: Name of the database (e.g., 'gtdb_bacteria', 'mgnify_human_gut')
        """
        self.database_name = database_name
        self.stats = {
            'attempted': 0,
            'successful': 0,
            'failed': 0,
            'skipped': 0
        }
    
    @abstractmethod
    def download_genome(
        self,
        genome_record: Dict,
        output_dir: Path,
        required_files: List[str] = None,
        compress: bool = True,
        terminal_log: 'DownloadLogger' = None
    ) -> bool:
        """
        Download genome files for a single genome.
        
        Args:
            genome_record: Dictionary containing genome metadata from input JSON
            output_dir: Base output directory (database subfolder will be created inside)
            required_files: List of required files (e.g., ['genome.fasta'])
            compress: Whether to compress output files (.gz)
            terminal_log: Logger for this specific download
        
        Returns:
            True if download successful, False otherwise
        
        Implementation should:
        1. Extract genome_id from genome_record
        2. Create genome-specific directory: output_dir / genome_id /
        3. Download required files
        4. Compress if requested
        5. Update terminal_log with progress
        """
        pass
    
    @abstractmethod
    def get_metadata(self, genome_record: Dict) -> Dict:
        """
        Extract and format metadata from genome record.
        
        Args:
            genome_record: Raw genome record from input JSON
        
        Returns:
            Dictionary with standardized metadata fields for manifest.json
        
        Should include at minimum:
        - genome_id: Unique identifier
        - database: Database name
        - taxonomy: Taxonomic lineage (if available)
        - quality_metrics: Completeness, contamination, etc. (if available)
        - Any other database-specific fields
        """
        pass
    
    @abstractmethod
    def validate_genome_id(self, genome_id: str) -> bool:
        """
        Validate if a genome ID is properly formatted for this database.
        
        Args:
            genome_id: Genome identifier to validate
        
        Returns:
            True if valid format, False otherwise
        
        Examples:
        - GTDB: 'GB_GCA_041569275.1' or 'RS_GCF_034719275.1'
        - MGnify: 'MGYG000000001'
        - NCBI: 'GCA_036600855.1'
        """
        pass
    
    def is_genome_complete(
        self,
        genome_dir: Path,
        required_files: List[str]
    ) -> bool:
        """
        Check if all required files exist and are valid.
        
        Args:
            genome_dir: Directory containing genome files
            required_files: List of required filenames
        
        Returns:
            True if all files present and valid, False otherwise
        """
        if not genome_dir.exists():
            return False
        
        for filename in required_files:
            # Check both compressed and uncompressed versions
            filepath = genome_dir / filename
            filepath_gz = genome_dir / f"{filename}.gz"
            
            if not (filepath.exists() or filepath_gz.exists()):
                return False
            
            # Check file size (sanity check for corruption)
            actual_file = filepath if filepath.exists() else filepath_gz
            if actual_file.stat().st_size < 100:  # Too small, likely corrupted
                return False
        
        return True
    
    def get_genome_id(self, genome_record: Dict) -> str:
        """
        Extract genome ID from genome record.
        
        Different databases use different field names:
        - GTDB: 'accession'
        - MGnify: 'genome_id'
        - NCBI: 'assembly_accession'
        
        Override this method in subclasses if needed.
        """
        # Try common field names
        for field in ['genome_id', 'accession', 'assembly_accession', 'Genome']:
            if field in genome_record:
                return genome_record[field]
        
        raise ValueError(f"Could not find genome ID in record: {genome_record.keys()}")
    
    def update_stats(self, status: str):
        """
        Update download statistics.
        
        Args:
            status: One of 'attempted', 'successful', 'failed', 'skipped'
        """
        if status in self.stats:
            self.stats[status] += 1
    
    def get_stats(self) -> Dict:
        """Return current statistics."""
        return self.stats.copy()
    
    def log_summary(self):
        """Log summary statistics."""
        logger.info(f"Download Summary for {self.database_name}:")
        logger.info(f"  Attempted: {self.stats['attempted']}")
        logger.info(f"  Successful: {self.stats['successful']}")
        logger.info(f"  Failed: {self.stats['failed']}")
        logger.info(f"  Skipped: {self.stats['skipped']}")
        
        if self.stats['attempted'] > 0:
            success_rate = (self.stats['successful'] / self.stats['attempted']) * 100
            logger.info(f"  Success Rate: {success_rate:.1f}%")