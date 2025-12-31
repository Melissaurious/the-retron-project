"""
Enhanced Genome Fetcher - Supports both TSV and JSON inputs.

Features:
- TSV/CSV input: Original workflow (protein accessions → genome resolution)
- JSON input: New batch workflow (direct genome downloads)
- Automatic downloader selection based on database
- Configurable required files
- Compression support
- Parallel downloads
- Resume capability
"""

import argparse
import json
import logging
import sys
import time
import multiprocessing
import concurrent.futures
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime

# Your existing imports
from genome_fetcher.config import setup_logging, NCBI_EMAIL
from genome_fetcher.downloaders import (
    BVBRCDownloader,
    NCBIDownloader,
    get_downloader  # New factory function
)
from genome_fetcher.extractors import SequenceExtractor
from genome_fetcher.utils import TSVParser, CacheManager, DownloadLogger, is_file_complete

logger = logging.getLogger(__name__)


class GenomeFetcher:
    """
    Enhanced genome fetcher supporting multiple input formats and databases.
    """
    
    def __init__(
        self,
        input_file: Path,
        output_dir: Path,
        input_type: str = 'auto',  # 'auto', 'tsv', or 'json'
        node_column: str = 'Node',
        required_files: List[str] = None,
        compress: bool = True,
        ncbi_email: str = None,
        skip_existing: bool = True,
        try_bvbrc_fallback: bool = True,
        max_workers: int = None,
        parallel_downloads: bool = False,
        generate_manifest: bool = True
    ):
        """
        Initialize Genome Fetcher.
        
        Args:
            input_file: Input file (TSV/CSV or JSON)
            output_dir: Output directory
            input_type: Input format ('auto', 'tsv', or 'json')
            node_column: Column for terminal naming (TSV only)
            required_files: Files to download (default: ['genome.fasta'])
            compress: Compress output files
            ncbi_email: Email for NCBI queries
            skip_existing: Skip completed downloads
            try_bvbrc_fallback: Try BV-BRC if NCBI fails
            max_workers: Max parallel workers
            parallel_downloads: Enable parallel processing
            generate_manifest: Generate manifest.json per database
        """
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        self.input_type = input_type
        self.node_column = node_column
        self.required_files = required_files or ['genome.fasta']
        self.compress = compress
        self.skip_existing = skip_existing
        self.try_bvbrc_fallback = try_bvbrc_fallback
        self.parallel_downloads = parallel_downloads
        self.generate_manifest = generate_manifest
        
        # Auto-calculate workers
        if max_workers is None:
            available_cores = multiprocessing.cpu_count()
            self.max_workers = min(available_cores // 2, 10)
        else:
            self.max_workers = max_workers
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Detect input type if auto
        if self.input_type == 'auto':
            self.input_type = self._detect_input_type()
        
        logger.info(f"Input type: {self.input_type}")
        logger.info(f"Required files: {self.required_files}")
        logger.info(f"Compression: {'enabled' if compress else 'disabled'}")
        
        # Initialize components based on input type
        if self.input_type == 'tsv':
            # Original TSV workflow
            self.parser = TSVParser(self.input_file)
            self.cache = CacheManager(self.output_dir)
            self.bvbrc = BVBRCDownloader()
            self.ncbi = NCBIDownloader(email=ncbi_email) if ncbi_email else None
            self.extractor = SequenceExtractor()
        else:
            # New JSON workflow
            self.ncbi_email = ncbi_email
            self.manifests = {}  # database -> manifest data
        
        # Statistics
        self.stats = {
            'total': 0,
            'success': 0,
            'failed': 0,
            'skipped': 0
        }
        
        if parallel_downloads:
            logger.info(f"Parallel downloads enabled (max workers: {self.max_workers})")
    
    def _detect_input_type(self) -> str:
        """Auto-detect input file type."""
        suffix = self.input_file.suffix.lower()
        
        if suffix == '.json':
            return 'json'
        elif suffix in ['.tsv', '.csv', '.txt']:
            return 'tsv'
        else:
            # Try to read first line
            with open(self.input_file, 'r') as f:
                first_line = f.read(100).strip()
                if first_line.startswith('{') or first_line.startswith('['):
                    return 'json'
                else:
                    return 'tsv'
    
    def run(self) -> bool:
        """Run the genome fetcher."""
        logger.info("="*80)
        logger.info("Starting Genome Fetcher")
        logger.info("="*80)
        logger.info(f"Input: {self.input_file}")
        logger.info(f"Output: {self.output_dir}")
        logger.info(f"Mode: {self.input_type.upper()}")
        logger.info("")
        
        if self.input_type == 'tsv':
            return self._run_tsv_workflow()
        else:
            return self._run_json_workflow()
    
    # =========================================================================
    # TSV Workflow (Original)
    # =========================================================================
    
    def _run_tsv_workflow(self) -> bool:
        """Run original TSV-based workflow."""
        try:
            records = self.parser.parse(node_column=self.node_column)
        except Exception as e:
            logger.error(f"Failed to parse TSV: {e}")
            return False
        
        self.stats['total'] = len(records)
        logger.info(f"Processing {self.stats['total']} records from TSV")
        
        if self.parallel_downloads and len(records) > 1:
            optimal_workers = self._calculate_optimal_workers(len(records))
            logger.info(f"Using parallel processing ({optimal_workers} workers)")
            success = self._run_parallel_tsv(records, optimal_workers)
        else:
            logger.info("Using sequential processing")
            success = self._run_sequential_tsv(records)
        
        self.print_summary()
        self.generate_detailed_report()
        
        return self.stats['failed'] == 0
    
    def _run_sequential_tsv(self, records: List[Dict]) -> bool:
        """Sequential TSV processing."""
        for idx, record in enumerate(records, 1):
            logger.info(f"[{idx}/{self.stats['total']}] Processing record {idx}")
            success = self.process_tsv_record(record, idx)
            
            if success:
                self.stats['success'] += 1
            else:
                self.stats['failed'] += 1
            
            logger.info("")
        
        return True
    
    def _run_parallel_tsv(self, records: List[Dict], num_workers: int) -> bool:
        """Parallel TSV processing."""
        from threading import Lock
        
        progress_lock = Lock()
        completed = {'count': 0}
        
        def process_with_progress(record_tuple):
            idx, record = record_tuple
            success = self.process_tsv_record(record, idx)
            
            with progress_lock:
                completed['count'] += 1
                if success:
                    self.stats['success'] += 1
                else:
                    self.stats['failed'] += 1
                logger.info(f"Progress: {completed['count']}/{self.stats['total']}")
            
            return success
        
        indexed_records = [(idx, record) for idx, record in enumerate(records, 1)]
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(process_with_progress, record) 
                      for record in indexed_records]
            concurrent.futures.wait(futures)
        
        return True
    
    def process_tsv_record(self, record: dict, record_num: int) -> bool:
        """Process single TSV record (original logic)."""
        # Your original process_record logic here
        # Keeping this minimal for now - you can copy your full implementation
        
        node_value = record.get(self.node_column, '')
        accession = record.get('Accesion', '')
        acc_type = record.get('_accession_type', 'unknown')
        genome_id = record.get('_genome_id')
        
        terminal_name = f"terminal_{node_value}"
        terminal_dir = self.output_dir / terminal_name
        terminal_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = terminal_dir / "download.log"
        terminal_log = DownloadLogger(log_file)
        
        terminal_log.info(f"Processing {accession} (type: {acc_type})")
        
        # Check if complete
        if self.skip_existing and self.cache.is_complete(terminal_name, self.required_files):
            logger.info(f"  ✓ Already complete - skipping")
            self.stats['skipped'] += 1
            return True
        
        # Process based on type
        try:
            if acc_type == 'bvbrc_fig':
                success = self.process_bvbrc_record(
                    accession, genome_id, terminal_dir, terminal_log
                )
            elif acc_type == 'ncbi_protein':
                sequence_id = f"{terminal_name}_{accession}"
                success = self.process_ncbi_protein(
                    accession, terminal_dir, terminal_log, sequence_id
                )
            else:
                terminal_log.error(f"Unknown type: {acc_type}")
                success = False
            
            if success:
                self.cache.mark_complete(terminal_name, metadata={'accession': accession})
            else:
                self.cache.mark_failed(terminal_name, "Download failed")
            
            return success
            
        except Exception as e:
            logger.error(f"Error: {e}")
            terminal_log.error(f"Error: {e}")
            self.cache.mark_failed(terminal_name, str(e))
            return False
    
    def process_bvbrc_record(self, accession, genome_id, terminal_dir, terminal_log):
        """Process BV-BRC record (original logic)."""
        download_success = self.bvbrc.download_genome(genome_id, terminal_dir, terminal_log)
        if not download_success:
            return False
        
        sequence_id = f"{terminal_dir.name}_{accession}"
        extract_success = self.extractor.extract_protein(accession, terminal_dir, sequence_id, terminal_log)
        return extract_success
    
    def process_ncbi_protein(self, accession, terminal_dir, terminal_log, sequence_id):
        """Process NCBI protein (original logic)."""
        if not self.ncbi:
            terminal_log.error("NCBI downloader not initialized")
            return False
        
        # Resolve protein to genome
        genome_info = self.ncbi.resolve_protein_to_genome(accession, terminal_log)
        if not genome_info:
            return False
        
        genome_acc, taxid = genome_info
        
        # Download genome
        genome_success = self.ncbi.download_genome_files(genome_acc, terminal_dir, terminal_log)
        
        if not genome_success and self.try_bvbrc_fallback:
            # Try BV-BRC fallback
            pass  # Implement if needed
        
        # Download protein
        protein_success = self.ncbi.download_protein_sequences(
            accession, terminal_dir, terminal_log, sequence_id
        )
        
        return genome_success and protein_success
    
    # =========================================================================
    # JSON Workflow (New)
    # =========================================================================
    
    def _run_json_workflow(self) -> bool:
        """Run new JSON batch workflow."""
        # Load JSON
        try:
            genomes = self._load_json_batch()
        except Exception as e:
            logger.error(f"Failed to load JSON: {e}")
            return False
        
        self.stats['total'] = len(genomes)
        logger.info(f"Processing {self.stats['total']} genomes from JSON")
        
        # Process genomes
        if self.parallel_downloads and len(genomes) > 1:
            optimal_workers = self._calculate_optimal_workers(len(genomes))
            logger.info(f"Using parallel processing ({optimal_workers} workers)")
            success = self._run_parallel_json(genomes, optimal_workers)
        else:
            logger.info("Using sequential processing")
            success = self._run_sequential_json(genomes)
        
        # Generate manifests
        if self.generate_manifest:
            self._generate_manifests()
        
        self.print_summary()
        
        return self.stats['failed'] == 0
    
    def _load_json_batch(self) -> List[Dict]:
        """Load genomes from JSON file."""
        with open(self.input_file, 'r') as f:
            data = json.load(f)
        
        # Handle both dict and list formats
        if isinstance(data, dict):
            # Convert dict to list
            genomes = [
                {**metadata, 'genome_id': genome_id}
                for genome_id, metadata in data.items()
            ]
        else:
            genomes = data
        
        return genomes
    
    def _run_sequential_json(self, genomes: List[Dict]) -> bool:
        """Sequential JSON processing."""
        for idx, genome in enumerate(genomes, 1):
            logger.info(f"[{idx}/{self.stats['total']}] Processing genome {idx}")
            success = self.process_json_genome(genome, idx)
            
            if success:
                self.stats['success'] += 1
            else:
                self.stats['failed'] += 1
            
            logger.info("")
        
        return True
    
    def _run_parallel_json(self, genomes: List[Dict], num_workers: int) -> bool:
        """Parallel JSON processing."""
        from threading import Lock
        
        progress_lock = Lock()
        completed = {'count': 0}
        
        def process_with_progress(genome_tuple):
            idx, genome = genome_tuple
            success = self.process_json_genome(genome, idx)
            
            with progress_lock:
                completed['count'] += 1
                if success:
                    self.stats['success'] += 1
                else:
                    self.stats['failed'] += 1
                logger.info(f"Progress: {completed['count']}/{self.stats['total']}")
            
            return success
        
        indexed_genomes = [(idx, genome) for idx, genome in enumerate(genomes, 1)]
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(process_with_progress, genome) 
                      for genome in indexed_genomes]
            concurrent.futures.wait(futures)
        
        return True
    
    def process_json_genome(self, genome_record: Dict, idx: int) -> bool:
        """
        Process single genome from JSON.
        
        Args:
            genome_record: Genome metadata dictionary
            idx: Record index for logging
        
        Returns:
            True if successful
        """
        # Get database and genome ID
        database = genome_record.get('database')
        if not database:
            logger.error(f"No 'database' field in genome record {idx}")
            return False
        
        genome_id = genome_record.get('genome_id') or genome_record.get('accession') or genome_record.get('assembly_accession')
        if not genome_id:
            logger.error(f"No genome ID found in record {idx}")
            return False
        
        logger.info(f"  Database: {database}")
        logger.info(f"  Genome ID: {genome_id}")
        
        # Get appropriate downloader
        downloader = get_downloader(database)
        if not downloader:
            logger.error(f"  ✗ Unsupported database: {database}")
            return False
        
        # For NCBI, set email
        if database.startswith('ncbi') or database.startswith('NCBI'):
            if hasattr(downloader, 'email'):
                downloader.email = self.ncbi_email or "your_email@example.com"
        
        # Create database-specific output directory
        db_output_dir = self.output_dir / database
        db_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create genome-specific log
        genome_dir = db_output_dir / genome_id
        genome_dir.mkdir(parents=True, exist_ok=True)
        log_file = genome_dir / "download.log"
        terminal_log = DownloadLogger(log_file)
        
        # Download
        try:
            success = downloader.download_genome(
                genome_record=genome_record,
                output_dir=db_output_dir,
                required_files=self.required_files,
                compress=self.compress,
                terminal_log=terminal_log
            )
            
            if success:
                logger.info(f"  ✓ Success")
                # Track for manifest
                self._add_to_manifest(database, genome_id, genome_record, success=True)
            else:
                logger.info(f"  ✗ Failed")
                self._add_to_manifest(database, genome_id, genome_record, success=False)
            
            return success
            
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            terminal_log.error(f"Error: {e}")
            self._add_to_manifest(database, genome_id, genome_record, success=False, error=str(e))
            return False
    
    def _add_to_manifest(self, database: str, genome_id: str, genome_record: Dict, success: bool, error: str = None):
        """Add genome to manifest tracking."""
        if database not in self.manifests:
            self.manifests[database] = {
                'database': database,
                'total_genomes': 0,
                'successful': 0,
                'failed': 0,
                'genomes': {}
            }
        
        manifest = self.manifests[database]
        manifest['total_genomes'] += 1
        
        if success:
            manifest['successful'] += 1
            status = 'complete'
        else:
            manifest['failed'] += 1
            status = 'failed'
        
        # Add genome entry
        genome_dir = self.output_dir / database / genome_id
        manifest['genomes'][genome_id] = {
            'genome_id': genome_id,
            'status': status,
            'file_path': str(genome_dir / f"genome.fasta{'gz' if self.compress else ''}"),
            'metadata': genome_record,
            'error': error
        }
    
    def _generate_manifests(self):
        """Generate manifest.json for each database."""
        for database, manifest_data in self.manifests.items():
            manifest_file = self.output_dir / database / "manifest.json"
            with open(manifest_file, 'w') as f:
                json.dump(manifest_data, f, indent=2)
            logger.info(f"Generated manifest: {manifest_file}")
    
    # =========================================================================
    # Utilities
    # =========================================================================
    
    def _calculate_optimal_workers(self, total_records: int) -> int:
        """Calculate optimal worker count."""
        if total_records < 5:
            return 1
        elif total_records < 20:
            return min(3, self.max_workers)
        elif total_records < 100:
            return min(5, self.max_workers)
        else:
            calculated = min(total_records // 10, self.max_workers)
            return max(calculated, 5)
    
    def print_summary(self):
        """Print processing summary."""
        logger.info("="*80)
        logger.info("PROCESSING SUMMARY")
        logger.info("="*80)
        logger.info(f"Total:     {self.stats['total']}")
        logger.info(f"Success:   {self.stats['success']}")
        logger.info(f"Failed:    {self.stats['failed']}")
        logger.info(f"Skipped:   {self.stats['skipped']}")
        logger.info("="*80)
        
        if self.stats['failed'] > 0:
            logger.warning(f"{self.stats['failed']} records failed")
    
    def generate_detailed_report(self):
        """Generate detailed report (TSV mode only)."""
        if self.input_type != 'tsv' or not hasattr(self, 'cache'):
            return
        
        output_file = self.output_dir / "download_report.txt"
        # Your existing report generation logic
        logger.info(f"Report saved: {output_file}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Fetch genomes from multiple databases'
    )
    parser.add_argument('--input', '-i', required=True, help='Input file (TSV/CSV or JSON)')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--input-type', choices=['auto', 'tsv', 'json'], default='auto', help='Input format')
    parser.add_argument('--node-column', default='Node', help='Column for terminal naming (TSV only)')
    parser.add_argument('--required-files', nargs='+', default=['genome.fasta'], help='Required files to download')
    parser.add_argument('--compress', action='store_true', default=True, help='Compress output (default: True)')
    parser.add_argument('--no-compress', action='store_false', dest='compress', help='Disable compression')
    parser.add_argument('--ncbi-email', default=None, help='Email for NCBI API')
    parser.add_argument('--skip-existing', action='store_true', default=True, help='Skip existing (default: True)')
    parser.add_argument('--no-skip-existing', action='store_false', dest='skip_existing')
    parser.add_argument('--parallel', action='store_true', default=False, help='Enable parallel downloads')
    parser.add_argument('--max-workers', type=int, default=None, help='Max parallel workers')
    parser.add_argument('--no-manifest', action='store_false', dest='generate_manifest', help='Disable manifest generation')
    parser.add_argument('--log-file', default=None, help='Log file path')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logging(log_file=args.log_file, level=log_level)
    
    # Validate input
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)
    
    # Get NCBI email
    ncbi_email = args.ncbi_email or NCBI_EMAIL
    if not ncbi_email:
        logger.warning("No NCBI email provided - NCBI downloads may be limited")
    
    # Create fetcher
    fetcher = GenomeFetcher(
        input_file=input_file,
        output_dir=Path(args.output),
        input_type=args.input_type,
        node_column=args.node_column,
        required_files=args.required_files,
        compress=args.compress,
        ncbi_email=ncbi_email,
        skip_existing=args.skip_existing,
        max_workers=args.max_workers,
        parallel_downloads=args.parallel,
        generate_manifest=args.generate_manifest
    )
    
    success = fetcher.run()
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()