"""
Configuration file for genome fetcher
"""
import logging
from pathlib import Path

# NCBI Configuration
NCBI_EMAIL = "melissa.rioszertuche@kaust.edu.sa"  # Users should update this
NCBI_TOOL = "genome_fetcher"
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_RATE_LIMIT = 0.34  # ~3 requests per second (NCBI limit)

# BV-BRC Configuration
BVBRC_FTP_BASE = "ftp://ftp.bvbrc.org/genomes"
BVBRC_API_BASE = "https://www.bv-brc.org/api"

# File naming constants
REQUIRED_FILES = [
    "genome.fasta",
    "genome.faa",
    "genome.ffn",
    "protein_aminoacid.fasta",
    "protein_nucleotide.fasta"
]

OPTIONAL_FILES = [
    "genome.features.tab",
]

# Logging configuration
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'

def setup_logging(log_file=None, level=logging.INFO):
    """Setup logging configuration"""
    handlers = [logging.StreamHandler()]
    
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format=LOG_FORMAT,
        datefmt=LOG_DATE_FORMAT,
        handlers=handlers,
        force=True
    )
    
    return logging.getLogger(__name__)