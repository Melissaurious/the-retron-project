#!/usr/bin/env python3
"""
Batch Processor - Download genomes in batches with progress tracking.

Features:
- Processes batches created by batch_splitter.py
- Integrates with existing downloaders (gtdb, mgnify, ncbi, etc.)
- Automatic resume capability (reads completed batches)
- Simple progress tracking (text files + JSON stats)
- Parallel downloads within batch (optional)
- Retry failed downloads

Usage:
    # Process all batches
    python batch_processor.py \\
        --database gtdb_bacteria \\
        --batch-dir batches/gtdb_bacteria \\
        --output-dir downloads

    # Process single batch (for SLURM array jobs)
    python batch_processor.py \\
        --database gtdb_bacteria \\
        --batch-dir batches/gtdb_bacteria \\
        --output-dir downloads \\
        --batch-id 5 \\
        --single-batch

    # Retry failed downloads
    python batch_processor.py \\
        --database gtdb_bacteria \\
        --batch-dir batches/gtdb_bacteria \\
        --output-dir downloads \\
        --retry-failed
"""

import json
import argparse
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Optional
import sys
import logging

# Add parent directory to path for imports
script_dir = Path(__file__).resolve().parent
if script_dir.name == 'genome_fetcher':
    sys.path.insert(0, str(script_dir.parent))
else:
    sys.path.insert(0, str(script_dir))

from genome_fetcher.downloaders import get_downloader
from genome_fetcher.utils import DownloadLogger

logger = logging.getLogger(__name__)


def parse_taxonomy(genome_data: Dict) -> tuple:
    """
    Parse taxonomy from genome metadata.
    
    Extracts genus and species from various taxonomy formats:
    - GTDB: "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
    - MGnify: Similar to GTDB (Lineage field)
    - NCBI: "Escherichia coli" or "Escherichia coli str. K-12" (organism_name field)
    
    Args:
        genome_data: Genome metadata dictionary
    
    Returns:
        Tuple of (genus, species) or (None, None) if not found
    """
    # Try GTDB/MGnify style (Lineage field with GTDB format)
    lineage = genome_data.get('Lineage') or genome_data.get('lineage') or genome_data.get('taxonomy')
    
    if lineage and isinstance(lineage, str):
        # Parse GTDB format: "...;g__Escherichia;s__Escherichia_coli"
        parts = lineage.split(';')
        
        genus = None
        species = None
        
        for part in parts:
            if part.startswith('g__'):
                genus = part[3:].strip().replace(' ', '_')
            elif part.startswith('s__'):
                species = part[3:].strip().replace(' ', '_')
        
        if genus and species and species != genus:
            return (genus, species)
        elif genus:
            # Has genus but no species
            return (genus, 'unclassified')
    
    # Try NCBI style (organism_name field)
    organism_name = genome_data.get('organism_name') or genome_data.get('organism')
    
    if organism_name and isinstance(organism_name, str):
        # Parse "Genus species" or "Genus species strain ..."
        words = organism_name.split()
        
        # Filter out common non-taxonomic words
        skip_words = {'sp.', 'str.', 'strain', 'subsp.', 'subspecies', 'serovar', 'uncultured', 'bacterium', 'archaeon'}
        
        # Get first two valid words
        valid_words = [w for w in words if w not in skip_words and not w.startswith('(')]
        
        if len(valid_words) >= 2:
            genus = valid_words[0].replace(' ', '_')
            species = valid_words[1].replace(' ', '_')
            return (genus, species)
        elif len(valid_words) == 1:
            genus = valid_words[0].replace(' ', '_')
            return (genus, 'unclassified')
    
    # No taxonomy found
    return (None, None)


class BatchProcessor:
    """Process genome download batches with progress tracking."""
    
    def __init__(
        self,
        database: str,
        batch_dir: Path,
        output_dir: Path,
        skip_existing: bool = True,
        compress: bool = True,
        max_retries: int = 3,
        group_by_taxonomy: bool = False
    ):
        """
        Initialize batch processor.
        
        Args:
            database: Database name (e.g., 'gtdb_bacteria')
            batch_dir: Directory containing batches and manifest
            output_dir: Output directory for downloads
            skip_existing: Skip already downloaded genomes
            compress: Compress output files
            max_retries: Maximum retry attempts for failures
            group_by_taxonomy: Organize genomes by genus/species (hybrid mode)
        """
        self.database = database
        self.batch_dir = Path(batch_dir)
        self.output_dir = Path(output_dir)
        self.skip_existing = skip_existing
        self.compress = compress
        self.max_retries = max_retries
        self.group_by_taxonomy = group_by_taxonomy
        
        # Create output directories
        self.genomes_dir = self.output_dir / database / 'genomes'
        self.progress_dir = self.output_dir / database / '.progress'
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.progress_dir.mkdir(parents=True, exist_ok=True)
        
        # Progress tracking files
        self.completed_batches_file = self.progress_dir / 'completed_batches.txt'
        self.failed_genomes_file = self.progress_dir / 'failed_genomes.jsonl'
        self.stats_file = self.progress_dir / 'stats.json'
        self.last_batch_file = self.progress_dir / 'last_batch.txt'
        
        # Load manifest
        self.manifest_path = self.batch_dir / 'manifest.json'
        self.manifest = self._load_manifest()
        
        # Get downloader
        self.downloader = get_downloader(database)
        if not self.downloader:
            raise ValueError(f"No downloader available for database: {database}")
        
        # Set NCBI email if needed
        if database.startswith('ncbi') and hasattr(self.downloader, 'email'):
            self.downloader.email = "genome_fetcher@example.com"
        
        # Load progress
        self.completed_batches = self._load_completed_batches()
        self.stats = self._load_stats()
    
    def _load_manifest(self) -> Dict:
        """Load batch manifest."""
        if not self.manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found: {self.manifest_path}")
        
        with open(self.manifest_path) as f:
            return json.load(f)
    
    def _load_completed_batches(self) -> Set[int]:
        """Load set of completed batch IDs."""
        if not self.completed_batches_file.exists():
            return set()
        
        with open(self.completed_batches_file) as f:
            return {int(line.strip()) for line in f if line.strip()}
    
    def _load_stats(self) -> Dict:
        """Load statistics."""
        if not self.stats_file.exists():
            return {
                'total_genomes': self.manifest.get('total_genomes', 0),
                'completed': 0,
                'failed': 0,
                'skipped': 0,
                'total_size_bytes': 0,
                'start_time': None,
                'last_update': None
            }
        
        with open(self.stats_file) as f:
            return json.load(f)
    
    def _save_stats(self):
        """Save statistics."""
        self.stats['last_update'] = datetime.now().isoformat()
        with open(self.stats_file, 'w') as f:
            json.dump(self.stats, f, indent=2)
    
    def _mark_batch_completed(self, batch_id: int):
        """Mark a batch as completed."""
        with open(self.completed_batches_file, 'a') as f:
            f.write(f"{batch_id}\n")
        
        self.completed_batches.add(batch_id)
        
        # Update last batch marker
        with open(self.last_batch_file, 'w') as f:
            f.write(str(batch_id))
    
    def _record_failed_genome(self, genome_id: str, batch_id: int, error: str, attempts: int = 1):
        """Record a failed genome download."""
        record = {
            'genome_id': genome_id,
            'batch_id': batch_id,
            'error': error,
            'attempts': attempts,
            'timestamp': datetime.now().isoformat()
        }
        
        with open(self.failed_genomes_file, 'a') as f:
            f.write(json.dumps(record) + '\n')
    
    def _get_output_path(self, genome_id: str, genome_data: Dict = None) -> Path:
        """
        Get output path for a genome.
        
        Supports two modes:
        1. Prefix-based (default): genomes/RS_G/RS_GCF_034719275.1/
        2. Taxonomy-based (hybrid): genomes/g__Escherichia/s__Escherichia_coli/RS_G/RS_GCF_034719275.1/
        
        In taxonomy mode, uses sub-prefix for large species to avoid filesystem issues.
        
        Args:
            genome_id: Genome identifier
            genome_data: Genome metadata (optional, needed for taxonomy mode)
        
        Returns:
            Path to genome directory
        """
        if self.group_by_taxonomy and genome_data:
            # Try to parse taxonomy
            genus, species = parse_taxonomy(genome_data)
            
            if genus and species:
                # Taxonomy-based path with sub-prefix
                # Example: genomes/g__Escherichia/s__Escherichia_coli/RS_G/RS_GCF_034719275.1/
                prefix = genome_id[:4] if len(genome_id) >= 4 else genome_id
                
                # Add 'g__' and 's__' prefixes for clarity (GTDB style)
                genus_dir = f"g__{genus}" if not genus.startswith('g__') else genus
                species_dir = f"s__{species}" if not species.startswith('s__') else species
                
                genome_dir = self.genomes_dir / genus_dir / species_dir / prefix / genome_id
            else:
                # No taxonomy - fall back to prefix-based
                prefix = genome_id[:4] if len(genome_id) >= 4 else genome_id
                genome_dir = self.genomes_dir / 'unclassified' / prefix / genome_id
        else:
            # Simple prefix-based organization
            # Example: genomes/RS_G/RS_GCF_034719275.1/
            prefix = genome_id[:4] if len(genome_id) >= 4 else genome_id
            genome_dir = self.genomes_dir / prefix / genome_id
        
        return genome_dir
    
    def _is_genome_downloaded(self, genome_dir: Path) -> bool:
        """Check if genome is already downloaded."""
        # Check for genome.fasta or genome.fasta.gz
        fasta_file = genome_dir / 'genome.fasta'
        fasta_gz_file = genome_dir / 'genome.fasta.gz'
        
        if fasta_file.exists() and fasta_file.stat().st_size > 1000:
            return True
        if fasta_gz_file.exists() and fasta_gz_file.stat().st_size > 1000:
            return True
        
        return False
    
    def process_batch(self, batch_id: int) -> Dict:
        """
        Process a single batch.
        
        Args:
            batch_id: Batch ID to process
        
        Returns:
            Dictionary with batch statistics
        """
        batch_file = self.batch_dir / f"batch_{batch_id:04d}.json"
        
        if not batch_file.exists():
            logger.error(f"Batch file not found: {batch_file}")
            return {'error': 'file_not_found'}
        
        # Load batch
        with open(batch_file) as f:
            genomes = json.load(f)
        
        print(f"\n{'='*60}")
        print(f"Processing Batch {batch_id:04d}")
        print(f"{'='*60}")
        print(f"Database: {self.database}")
        print(f"Genomes in batch: {len(genomes):,}")
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*60}\n")
        
        # Track batch stats
        batch_stats = {
            'batch_id': batch_id,
            'total': len(genomes),
            'completed': 0,
            'failed': 0,
            'skipped': 0,
            'start_time': time.time()
        }
        
        # Process each genome
        for idx, (genome_id, genome_data) in enumerate(genomes.items(), 1):
            # Progress indicator
            if idx % 100 == 0:
                elapsed = time.time() - batch_stats['start_time']
                rate = idx / elapsed if elapsed > 0 else 0
                eta = (len(genomes) - idx) / rate if rate > 0 else 0
                print(f"  Progress: {idx:,}/{len(genomes):,} ({idx/len(genomes)*100:.1f}%) "
                      f"| Rate: {rate:.1f} genomes/s | ETA: {eta/60:.1f} min")
            
            genome_dir = self._get_output_path(genome_id, genome_data)
            
            # Skip if already downloaded
            if self.skip_existing and self._is_genome_downloaded(genome_dir):
                batch_stats['skipped'] += 1
                self.stats['skipped'] += 1
                continue
            
            # Ensure genome directory exists
            genome_dir.mkdir(parents=True, exist_ok=True)
            
            # Setup logging
            log_file = genome_dir / 'download.log'
            terminal_log = DownloadLogger(log_file)
            
            try:
                # Inject database field if not present
                if 'database' not in genome_data:
                    genome_data['database'] = self.database
                
                # Download genome
                success = self.downloader.download_genome(
                    genome_record=genome_data,
                    output_dir=genome_dir.parent,  # Parent directory (without genome_id)
                    required_files=['genome.fasta'],
                    compress=self.compress,
                    terminal_log=terminal_log
                )
                
                if success:
                    batch_stats['completed'] += 1
                    self.stats['completed'] += 1
                    
                    # Update size
                    fasta_file = genome_dir / ('genome.fasta.gz' if self.compress else 'genome.fasta')
                    if fasta_file.exists():
                        self.stats['total_size_bytes'] += fasta_file.stat().st_size
                else:
                    batch_stats['failed'] += 1
                    self.stats['failed'] += 1
                    self._record_failed_genome(genome_id, batch_id, "Download returned False", 1)
            
            except Exception as e:
                batch_stats['failed'] += 1
                self.stats['failed'] += 1
                self._record_failed_genome(genome_id, batch_id, str(e), 1)
                logger.error(f"Error downloading {genome_id}: {e}")
        
        # Batch summary
        batch_stats['end_time'] = time.time()
        batch_stats['duration'] = batch_stats['end_time'] - batch_stats['start_time']
        
        print(f"\n{'='*60}")
        print(f"Batch {batch_id:04d} Complete")
        print(f"{'='*60}")
        print(f"Completed: {batch_stats['completed']:,}")
        print(f"Failed: {batch_stats['failed']:,}")
        print(f"Skipped: {batch_stats['skipped']:,}")
        print(f"Duration: {batch_stats['duration']/60:.1f} minutes")
        print(f"Rate: {batch_stats['total']/batch_stats['duration']:.2f} genomes/second")
        print(f"{'='*60}\n")
        
        return batch_stats
    
    def run(self, batch_id: Optional[int] = None, single_batch: bool = False):
        """
        Run batch processing.
        
        Args:
            batch_id: Specific batch ID to process (for SLURM array jobs)
            single_batch: Only process one batch then exit
        """
        if self.stats['start_time'] is None:
            self.stats['start_time'] = datetime.now().isoformat()
        
        # Get batches to process
        all_batches = self.manifest.get('batches', [])
        
        if batch_id is not None:
            # Process specific batch
            batches_to_process = [b for b in all_batches if b['batch_id'] == batch_id]
            if not batches_to_process:
                print(f"❌ Batch {batch_id} not found in manifest")
                return 1
        else:
            # Process all incomplete batches
            batches_to_process = [
                b for b in all_batches 
                if b['batch_id'] not in self.completed_batches
            ]
        
        if not batches_to_process:
            print(f"✅ All batches already completed!")
            self._print_overall_progress()
            return 0
        
        print(f"\n{'='*60}")
        print(f"Batch Processor - {self.database}")
        print(f"{'='*60}")
        print(f"Total batches: {len(all_batches)}")
        print(f"Completed batches: {len(self.completed_batches)}")
        print(f"Batches to process: {len(batches_to_process)}")
        print(f"Skip existing: {self.skip_existing}")
        print(f"Compress: {self.compress}")
        print(f"Group by taxonomy: {self.group_by_taxonomy}")
        print(f"{'='*60}")
        
        # Process batches
        for idx, batch_info in enumerate(batches_to_process, 1):
            batch_id = batch_info['batch_id']
            
            print(f"\n[{idx}/{len(batches_to_process)}] Processing batch {batch_id:04d}...")
            
            batch_stats = self.process_batch(batch_id)
            
            # Mark as completed
            self._mark_batch_completed(batch_id)
            
            # Save stats
            self._save_stats()
            
            # Print progress
            self._print_overall_progress()
            
            # Exit if single batch mode
            if single_batch:
                print(f"\n✅ Single batch mode - exiting after batch {batch_id}")
                return 0
        
        # Final summary
        print(f"\n{'='*60}")
        print(f"✅ All batches completed!")
        print(f"{'='*60}")
        self._print_overall_progress()
        
        return 0
    
    def retry_failed(self):
        """Retry failed genome downloads."""
        if not self.failed_genomes_file.exists():
            print("No failed genomes to retry")
            return 0
        
        # Load failed genomes
        failed_genomes = []
        with open(self.failed_genomes_file) as f:
            for line in f:
                if line.strip():
                    failed_genomes.append(json.loads(line))
        
        # Filter by max retries
        to_retry = [
            g for g in failed_genomes 
            if g.get('attempts', 0) < self.max_retries
        ]
        
        if not to_retry:
            print(f"No genomes to retry (all exceeded {self.max_retries} attempts)")
            return 0
        
        print(f"\n{'='*60}")
        print(f"Retrying Failed Downloads")
        print(f"{'='*60}")
        print(f"Failed genomes: {len(failed_genomes)}")
        print(f"To retry: {len(to_retry)}")
        print(f"{'='*60}\n")
        
        # TODO: Implement retry logic
        # This would require loading the original genome metadata
        # For now, just show what would be retried
        
        for genome_info in to_retry[:10]:  # Show first 10
            print(f"  Would retry: {genome_info['genome_id']} (batch {genome_info['batch_id']})")
        
        if len(to_retry) > 10:
            print(f"  ... and {len(to_retry) - 10} more")
        
        print(f"\n⚠️  Retry functionality not yet implemented")
        print(f"To retry, re-run specific batches with: --batch-id <id>")
        
        return 0
    
    def _print_overall_progress(self):
        """Print overall progress."""
        total = self.stats['total_genomes']
        completed = self.stats['completed']
        failed = self.stats['failed']
        skipped = self.stats['skipped']
        
        processed = completed + failed + skipped
        pct = (processed / total * 100) if total > 0 else 0
        
        print(f"\n{'='*60}")
        print(f"Overall Progress")
        print(f"{'='*60}")
        print(f"Total genomes: {total:,}")
        print(f"Completed: {completed:,}")
        print(f"Failed: {failed:,}")
        print(f"Skipped: {skipped:,}")
        print(f"Progress: {processed:,}/{total:,} ({pct:.1f}%)")
        
        if self.stats['total_size_bytes'] > 0:
            size_gb = self.stats['total_size_bytes'] / (1024**3)
            print(f"Total size: {size_gb:.2f} GB")
        
        print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description='Process genome download batches',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all batches
  python batch_processor.py \\
      --database gtdb_bacteria \\
      --batch-dir batches/gtdb_bacteria \\
      --output-dir downloads

  # Process specific batch (for SLURM array jobs)
  python batch_processor.py \\
      --database gtdb_bacteria \\
      --batch-dir batches/gtdb_bacteria \\
      --output-dir downloads \\
      --batch-id 5 \\
      --single-batch

  # Retry failed downloads
  python batch_processor.py \\
      --database gtdb_bacteria \\
      --batch-dir batches/gtdb_bacteria \\
      --output-dir downloads \\
      --retry-failed
        """
    )
    
    parser.add_argument(
        '--database',
        required=True,
        help='Database name (e.g., gtdb_bacteria)'
    )
    
    parser.add_argument(
        '--batch-dir',
        type=Path,
        required=True,
        help='Directory containing batches and manifest'
    )
    
    parser.add_argument(
        '--output-dir',
        type=Path,
        required=True,
        help='Output directory for downloads'
    )
    
    parser.add_argument(
        '--batch-id',
        type=int,
        help='Process specific batch ID (for SLURM array jobs)'
    )
    
    parser.add_argument(
        '--single-batch',
        action='store_true',
        help='Process only one batch then exit (for SLURM jobs)'
    )
    
    parser.add_argument(
        '--retry-failed',
        action='store_true',
        help='Retry previously failed downloads'
    )
    
    parser.add_argument(
        '--no-skip-existing',
        action='store_true',
        help='Re-download existing files (default: skip)'
    )
    
    parser.add_argument(
        '--no-compress',
        action='store_true',
        help='Do not compress output files (default: compress)'
    )
    
    parser.add_argument(
        '--group-by-taxonomy',
        action='store_true',
        help='Organize genomes by genus/species (hybrid mode with prefix sub-directories)'
    )
    
    parser.add_argument(
        '--max-retries',
        type=int,
        default=3,
        help='Maximum retry attempts (default: 3)'
    )
    
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Create processor
    try:
        processor = BatchProcessor(
            database=args.database,
            batch_dir=args.batch_dir,
            output_dir=args.output_dir,
            skip_existing=not args.no_skip_existing,
            compress=not args.no_compress,
            max_retries=args.max_retries,
            group_by_taxonomy=args.group_by_taxonomy
        )
    except Exception as e:
        print(f"❌ Error initializing processor: {e}")
        return 1
    
    # Run
    if args.retry_failed:
        return processor.retry_failed()
    else:
        return processor.run(
            batch_id=args.batch_id,
            single_batch=args.single_batch
        )


if __name__ == '__main__':
    sys.exit(main())
