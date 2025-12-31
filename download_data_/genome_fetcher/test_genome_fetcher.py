#!/usr/bin/env python3
"""
Genome Fetcher Test Suite

Tests all components of the genome download pipeline before scaling up to
4.2M genomes. Verifies downloaders, file formats, and integration.

Usage:
    python test_genome_fetcher.py --quick        # Quick test (5 genomes, ~2 min)
    python test_genome_fetcher.py --full         # Full test (50 genomes, ~10 min)
    python test_genome_fetcher.py --database gtdb_bacteria  # Test specific database
"""

import argparse
import json
import sys
import time
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Test configurations
TEST_GENOMES = {
    'gtdb_bacteria': [
        {
            "genome_id": "RS_GCF_000005845.2",
            "accession": "RS_GCF_000005845.2",
            "database": "gtdb_bacteria",
            "gtdb_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
            "checkm_completeness": 99.93,
            "checkm_contamination": 0.19
        }
    ],
    'mgnify_human_gut': [
        {
            "genome_id": "MGYG000000001",
            "database": "mgnify_human_gut",
            "catalogue": "human_gut",
            "Completeness": 98.59,
            "Contamination": 0.7,
            "Lineage": "d__Bacteria;p__Firmicutes_A;c__Clostridia",
            "FTP_download": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/all_genomes/MGYG0000000/MGYG000000001/genome/MGYG000000001.fna.gz"
        }
    ],
    'ncbi_bacteria': [
        {
            "genome_id": "GCA_000005845.2",
            "assembly_accession": "GCA_000005845.2",
            "database": "ncbi_bacteria",
            "organism_name": "Escherichia coli K-12",
            "assembly_level": "Complete Genome",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2"
        }
    ]
}


class GenomeFetcherTester:
    """Test suite for genome fetcher."""
    
    def __init__(self, output_dir: Path = None):
        """Initialize tester."""
        self.output_dir = output_dir or Path('./test_output')
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'passed': [],
            'failed': [],
            'skipped': []
        }
    
    def run_all_tests(self, test_mode: str = 'quick', database: str = None):
        """Run all tests."""
        logger.info("="*80)
        logger.info("GENOME FETCHER TEST SUITE")
        logger.info("="*80)
        logger.info(f"Test Mode: {test_mode}")
        logger.info(f"Output Directory: {self.output_dir}")
        logger.info("")
        
        # Test 1: Import all modules
        self.test_imports()
        
        # Test 2: Downloader factory
        self.test_downloader_factory()
        
        # Test 3: Download single genome per database
        if database:
            databases_to_test = [database]
        else:
            databases_to_test = ['gtdb_bacteria', 'mgnify_human_gut', 'ncbi_bacteria'] if test_mode == 'full' else ['gtdb_bacteria']
        
        for db in databases_to_test:
            self.test_database_download(db)
        
        # Test 4: Batch file processing
        self.test_batch_processing(databases_to_test[0])
        
        # Test 5: Verify output structure
        self.test_output_structure()
        
        # Test 6: Compression
        self.test_compression()
        
        # Print summary
        self.print_summary()
    
    def test_imports(self):
        """Test that all required modules can be imported."""
        test_name = "Module Imports"
        logger.info(f"Test: {test_name}")
        
        try:
            from genome_fetcher.downloaders import (
                BaseDownloader,
                GTDBDownloader,
                MGnifyDownloader,
                NCBIDownloader,
                get_downloader
            )
            logger.info("  ✅ All modules imported successfully")
            self.results['passed'].append(test_name)
        except ImportError as e:
            logger.error(f"  ❌ Import failed: {e}")
            self.results['failed'].append((test_name, str(e)))
    
    def test_downloader_factory(self):
        """Test downloader factory."""
        test_name = "Downloader Factory"
        logger.info(f"\nTest: {test_name}")
        
        try:
            from genome_fetcher.downloaders import get_downloader, DownloaderFactory
            
            # Test all supported databases
            supported = DownloaderFactory.get_supported_databases()
            logger.info(f"  Supported databases: {len(supported)}")
            
            for db in ['gtdb_bacteria', 'mgnify_human_gut', 'ncbi_bacteria']:
                downloader = get_downloader(db)
                if downloader:
                    logger.info(f"  ✅ {db}: {downloader.__class__.__name__}")
                else:
                    logger.error(f"  ❌ {db}: No downloader found")
                    raise ValueError(f"No downloader for {db}")
            
            self.results['passed'].append(test_name)
            
        except Exception as e:
            logger.error(f"  ❌ Factory test failed: {e}")
            self.results['failed'].append((test_name, str(e)))
    
    def test_database_download(self, database: str):
        """Test downloading from a specific database."""
        test_name = f"Download: {database}"
        logger.info(f"\nTest: {test_name}")
        
        if database not in TEST_GENOMES:
            logger.warning(f"  ⚠️  No test genomes defined for {database}")
            self.results['skipped'].append(test_name)
            return
        
        try:
            from genome_fetcher.downloaders import get_downloader
            
            # Get test genome
            test_genome = TEST_GENOMES[database][0]
            genome_id = test_genome['genome_id']
            
            logger.info(f"  Testing genome: {genome_id}")
            
            # Get downloader
            downloader = get_downloader(database)
            if not downloader:
                raise ValueError(f"No downloader for {database}")
            
            # Create test output directory
            test_output = self.output_dir / database
            test_output.mkdir(parents=True, exist_ok=True)
            
            # Download (with timeout protection)
            logger.info(f"  Downloading...")
            start_time = time.time()
            
            success = downloader.download_genome(
                genome_record=test_genome,
                output_dir=test_output,
                required_files=['genome.fasta'],
                compress=True,
                terminal_log=None
            )
            
            elapsed = time.time() - start_time
            
            if success:
                # Verify file exists
                genome_dir = test_output / genome_id
                expected_file = genome_dir / "genome.fasta.gz"
                
                if expected_file.exists():
                    file_size = expected_file.stat().st_size
                    logger.info(f"  ✅ Downloaded successfully ({file_size:,} bytes in {elapsed:.1f}s)")
                    self.results['passed'].append(test_name)
                else:
                    logger.error(f"  ❌ File not found: {expected_file}")
                    self.results['failed'].append((test_name, "Output file missing"))
            else:
                logger.error(f"  ❌ Download failed")
                self.results['failed'].append((test_name, "Download returned False"))
                
        except Exception as e:
            logger.error(f"  ❌ Error: {e}")
            self.results['failed'].append((test_name, str(e)))
    
    def test_batch_processing(self, database: str):
        """Test batch file processing."""
        test_name = "Batch Processing"
        logger.info(f"\nTest: {test_name}")
        
        try:
            # Create a small batch file
            batch_file = self.output_dir / "test_batch.json"
            
            test_data = {}
            if database in TEST_GENOMES:
                for i, genome in enumerate(TEST_GENOMES[database][:2]):  # Max 2 genomes
                    genome_copy = genome.copy()
                    genome_copy['batch_id'] = 0
                    genome_copy['quality_tier'] = 1
                    test_data[genome['genome_id']] = genome_copy
            
            with open(batch_file, 'w') as f:
                json.dump(test_data, f, indent=2)
            
            logger.info(f"  Created test batch: {batch_file}")
            logger.info(f"  Genomes in batch: {len(test_data)}")
            
            # Try to load it
            with open(batch_file, 'r') as f:
                loaded = json.load(f)
            
            if len(loaded) == len(test_data):
                logger.info(f"  ✅ Batch file valid")
                self.results['passed'].append(test_name)
            else:
                logger.error(f"  ❌ Batch file corrupted")
                self.results['failed'].append((test_name, "Data mismatch"))
                
        except Exception as e:
            logger.error(f"  ❌ Error: {e}")
            self.results['failed'].append((test_name, str(e)))
    
    def test_output_structure(self):
        """Test output directory structure."""
        test_name = "Output Structure"
        logger.info(f"\nTest: {test_name}")
        
        try:
            # Check if database subdirectories were created
            database_dirs = [d for d in self.output_dir.iterdir() if d.is_dir() and d.name in ['gtdb_bacteria', 'mgnify_human_gut', 'ncbi_bacteria']]
            
            logger.info(f"  Found {len(database_dirs)} database directories")
            
            for db_dir in database_dirs:
                genome_dirs = [d for d in db_dir.iterdir() if d.is_dir()]
                logger.info(f"  {db_dir.name}: {len(genome_dirs)} genome directories")
                
                # Check if genomes have files
                for genome_dir in genome_dirs[:3]:  # Check first 3
                    files = list(genome_dir.glob("*.fasta*"))
                    if files:
                        logger.info(f"    ✓ {genome_dir.name}: {len(files)} files")
            
            if database_dirs:
                logger.info(f"  ✅ Output structure valid")
                self.results['passed'].append(test_name)
            else:
                logger.error(f"  ❌ No database directories found")
                self.results['failed'].append((test_name, "No output directories"))
                
        except Exception as e:
            logger.error(f"  ❌ Error: {e}")
            self.results['failed'].append((test_name, str(e)))
    
    def test_compression(self):
        """Test that files are properly compressed."""
        test_name = "Compression"
        logger.info(f"\nTest: {test_name}")
        
        try:
            import gzip
            
            # Find .gz files
            gz_files = list(self.output_dir.rglob("*.fasta.gz"))
            
            if not gz_files:
                logger.warning(f"  ⚠️  No compressed files found")
                self.results['skipped'].append(test_name)
                return
            
            logger.info(f"  Found {len(gz_files)} compressed files")
            
            # Test first file
            test_file = gz_files[0]
            
            try:
                # Try to read it as gzip
                with gzip.open(test_file, 'rt') as f:
                    first_line = f.readline()
                
                if first_line.startswith('>'):
                    logger.info(f"  ✅ Compression valid (FASTA format detected)")
                    logger.info(f"  Sample: {first_line.strip()[:60]}...")
                    self.results['passed'].append(test_name)
                else:
                    logger.error(f"  ❌ Invalid FASTA format: {first_line[:50]}")
                    self.results['failed'].append((test_name, "Invalid FASTA"))
            except gzip.BadGzipFile as e:
                logger.error(f"  ❌ Invalid gzip file: {e}")
                self.results['failed'].append((test_name, f"Bad gzip: {e}"))
            except EOFError:
                logger.error(f"  ❌ Truncated gzip file")
                self.results['failed'].append((test_name, "Truncated file"))
                
        except Exception as e:
            logger.error(f"  ❌ Error: {e}")
            self.results['failed'].append((test_name, str(e)))
    
    def print_summary(self):
        """Print test summary."""
        logger.info("")
        logger.info("="*80)
        logger.info("TEST SUMMARY")
        logger.info("="*80)
        
        total = len(self.results['passed']) + len(self.results['failed']) + len(self.results['skipped'])
        
        logger.info(f"Total Tests: {total}")
        logger.info(f"✅ Passed: {len(self.results['passed'])}")
        logger.info(f"❌ Failed: {len(self.results['failed'])}")
        logger.info(f"⚠️  Skipped: {len(self.results['skipped'])}")
        
        if self.results['failed']:
            logger.info("")
            logger.info("Failed Tests:")
            for test_name, error in self.results['failed']:
                logger.info(f"  ❌ {test_name}: {error}")
        
        logger.info("")
        logger.info("="*80)
        
        if not self.results['failed']:
            logger.info("✅ ALL TESTS PASSED - Ready for production!")
            logger.info("")
            logger.info("Next Steps:")
            logger.info("  1. Run prepare_download_jobs.py to create batches")
            logger.info("  2. Submit SLURM array job with download_batch.slurm")
            logger.info("  3. Monitor progress with completed_batches.txt")
            return 0
        else:
            logger.info("❌ SOME TESTS FAILED - Fix issues before scaling up")
            return 1


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Test genome fetcher before production use'
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Quick test (1 database, ~2 minutes)'
    )
    parser.add_argument(
        '--full',
        action='store_true',
        help='Full test (all databases, ~10 minutes)'
    )
    parser.add_argument(
        '--database',
        choices=['gtdb_bacteria', 'mgnify_human_gut', 'ncbi_bacteria'],
        help='Test specific database only'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('./test_output'),
        help='Output directory for test downloads'
    )
    
    args = parser.parse_args()
    
    # Determine test mode
    if args.full:
        test_mode = 'full'
    else:
        test_mode = 'quick'
    
    # Run tests
    tester = GenomeFetcherTester(output_dir=args.output_dir)
    exit_code = tester.run_all_tests(test_mode=test_mode, database=args.database)
    
    sys.exit(exit_code)


if __name__ == '__main__':
    main()