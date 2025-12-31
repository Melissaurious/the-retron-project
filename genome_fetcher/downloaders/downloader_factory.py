"""
Downloader factory for automatically selecting the correct downloader.

This module provides a simple interface to get the appropriate downloader
based on the database name specified in the genome metadata.
"""

import logging
from typing import Optional

from .base_downloader import BaseDownloader
from .gtdb_downloader import GTDBDownloader
from .mgnify_downloader import MGnifyDownloader
from .ncbi_downloader_simple import NCBIDownloaderSimple

logger = logging.getLogger(__name__)


class DownloaderFactory:
    """
    Factory for creating database-specific downloaders.
    
    Usage:
        downloader = DownloaderFactory.get_downloader('gtdb_bacteria')
        success = downloader.download_genome(genome_record, output_dir)
    """
    
    # Mapping of database names to downloader classes
    DOWNLOADER_MAP = {
        # GTDB databases
        'gtdb_bacteria': GTDBDownloader,
        'gtdb_archaea': GTDBDownloader,
        'GTDB': GTDBDownloader,
        
        # MGnify databases
        'mgnify_human_gut': MGnifyDownloader,
        'mgnify_marine': MGnifyDownloader,
        'mgnify_soil': MGnifyDownloader,
        'MGnify': MGnifyDownloader,
        
        # UHGG (uses same downloader as MGnify)
        'uhgg': MGnifyDownloader,
        'UHGG': MGnifyDownloader,
        
        # NCBI databases
        'ncbi_bacteria': NCBIDownloaderSimple,
        'ncbi_archaea': NCBIDownloaderSimple,
        'NCBI': NCBIDownloaderSimple,
    }
    
    @classmethod
    def get_downloader(cls, database_name: str) -> Optional[BaseDownloader]:
        """
        Get the appropriate downloader for a database.
        
        Args:
            database_name: Name of the database (e.g., 'gtdb_bacteria', 'mgnify_human_gut')
        
        Returns:
            Instantiated downloader or None if database not supported
        
        Example:
            >>> downloader = DownloaderFactory.get_downloader('gtdb_bacteria')
            >>> downloader.download_genome(genome_record, output_dir)
        """
        downloader_class = cls.DOWNLOADER_MAP.get(database_name)
        
        if downloader_class is None:
            logger.error(f"No downloader available for database: {database_name}")
            logger.info(f"Supported databases: {list(cls.DOWNLOADER_MAP.keys())}")
            return None
        
        # Instantiate and return
        return downloader_class(database_name=database_name)
    
    @classmethod
    def get_supported_databases(cls) -> list:
        """
        Get list of supported database names.
        
        Returns:
            List of supported database names
        """
        return list(cls.DOWNLOADER_MAP.keys())
    
    @classmethod
    def is_supported(cls, database_name: str) -> bool:
        """
        Check if a database is supported.
        
        Args:
            database_name: Database name to check
        
        Returns:
            True if supported
        """
        return database_name in cls.DOWNLOADER_MAP
    
    @classmethod
    def register_downloader(cls, database_name: str, downloader_class: type):
        """
        Register a custom downloader for a database.
        
        This allows users to add custom downloaders without modifying the factory.
        
        Args:
            database_name: Name to register under
            downloader_class: Downloader class (must inherit from BaseDownloader)
        
        Example:
            >>> class MyCustomDownloader(BaseDownloader):
            ...     pass
            >>> DownloaderFactory.register_downloader('my_database', MyCustomDownloader)
        """
        if not issubclass(downloader_class, BaseDownloader):
            raise ValueError(f"{downloader_class} must inherit from BaseDownloader")
        
        cls.DOWNLOADER_MAP[database_name] = downloader_class
        logger.info(f"Registered downloader for database: {database_name}")


def get_downloader(database_name: str) -> Optional[BaseDownloader]:
    """
    Convenience function to get a downloader.
    
    Args:
        database_name: Name of the database
    
    Returns:
        Instantiated downloader or None
    
    Example:
        >>> from genome_fetcher.downloaders import get_downloader
        >>> downloader = get_downloader('gtdb_bacteria')
    """
    return DownloaderFactory.get_downloader(database_name)