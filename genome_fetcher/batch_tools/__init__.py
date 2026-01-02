"""
Batch processing tools for large-scale genome downloads.

This module provides utilities for splitting large genome metadata files
into manageable batches and processing them with progress tracking and
resume capability.

Main components:
- batch_splitter: Split large JSON files into batches with optional SLURM job creation
- batch_processor: Download genomes in batches with automatic resume

Example usage:
    # Split into batches
    from batch_tools.batch_splitter import main as split
    
    # Process batches
    from batch_tools.batch_processor import BatchProcessor
    
    processor = BatchProcessor(
        database='gtdb_bacteria',
        batch_dir='batches/gtdb_bacteria',
        output_dir='downloads'
    )
    processor.run()
"""

__version__ = '1.0.0'
__author__ = 'Genome Fetcher Project'

# Export main classes for easy import
try:
    from .batch_processor import BatchProcessor
except ImportError:
    # Handle case where dependencies aren't available yet
    BatchProcessor = None

__all__ = ['BatchProcessor']
