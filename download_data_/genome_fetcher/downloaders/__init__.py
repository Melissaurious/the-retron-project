
"""
Genome downloaders for various databases.

This module provides a unified interface for downloading genomes from
different databases (GTDB, MGnify, NCBI, UHGG, etc.).

Usage:
    from genome_fetcher.downloaders import get_downloader
    
    downloader = get_downloader('gtdb_bacteria')
    success = downloader.download_genome(genome_record, output_dir)

Available downloaders:
- GTDBDownloader: For GTDB bacteria and archaea
- MGnifyDownloader: For MGnify catalogues (human gut, marine, soil)
- NCBIDownloaderSimple: For NCBI assemblies (FTP downloads)
- NCBIDownloader: Full-featured NCBI downloader (E-utilities, proteins, etc.)

The factory function get_downloader() automatically selects the correct
downloader based on the database name.
"""

from .base_downloader import BaseDownloader
from .gtdb_downloader import GTDBDownloader
from .mgnify_downloader import MGnifyDownloader
from .ncbi_downloader import NCBIDownloader
from .bvbrc_downloader import BVBRCDownloader
from .downloader_factory import DownloaderFactory, get_downloader

# try:
#     from .ncbi_downloader_simple import NCBIDownloaderSimple
# except ImportError:
#     NCBIDownloaderSimple = None  # Optional, not critical

# Import the full-featured NCBI downloader if it exists
# try:
#     from .ncbi_downloader import NCBIDownloader
# except ImportError:
#     NCBIDownloader = None

# # Import BVBRC downloader if it exists
# try:
#     from .bvbrc_downloader import BVBRCDownloader
# except ImportError:
#     BVBRCDownloader = None

__all__ = [
    # Base class
    'BaseDownloader',
    
    # Specific downloaders
    'GTDBDownloader',
    'MGnifyDownloader',
    'NCBIDownloaderSimple',
    'NCBIDownloader',
    'BVBRCDownloader',
    
    # Factory
    'DownloaderFactory',
    'get_downloader',
]

__version__ = '1.0.0'


# """
# Downloaders package
# """
# from .bvbrc_downloader import BVBRCDownloader
# from .ncbi_downloader import NCBIDownloader

# __all__ = ['BVBRCDownloader', 'NCBIDownloader']
