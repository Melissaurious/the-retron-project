# other imports
from genome_fetcher.config import setup_logging, REQUIRED_FILES, OPTIONAL_FILES,NCBI_EMAIL
from genome_fetcher.downloaders import BVBRCDownloader, NCBIDownloader
from genome_fetcher.extractors import SequenceExtractor
from genome_fetcher.utils import TSVParser, CacheManager, DownloadLogger, is_file_complete


logger = logging.getLogger(__name__)

REQUIRED_FILES = [
    'protein_aminoacid.fasta',
    'protein_nucleotide.fasta',
    'genome.fasta'
]

class GenomeFetcher:
    """Main orchestrator for genome fetching with enhanced NCBI support"""
    
    def __init__(
        self,
        input_file: Path,
        output_dir: Path,
        node_column: str = 'Node',
        ncbi_email: str = None,
        skip_existing: bool = True,
        try_bvbrc_fallback: bool = True,
        max_workers: int = None,  # ← ADD THIS
        parallel_downloads: bool = False  # ← ADD THIS
    ):

        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        self.node_column = node_column
        self.skip_existing = skip_existing
        self.try_bvbrc_fallback = try_bvbrc_fallback
        self.parallel_downloads = parallel_downloads


        # Calculate optimal workers
        if max_workers is None:
            # Auto-calculate based on CPU cores and expected workload
            available_cores = multiprocessing.cpu_count()
            self.max_workers = min(available_cores // 2, 10)  # Use half cores, max 10
        else:
            self.max_workers = max_workers
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        
        # Initialize components
        self.parser = TSVParser(self.input_file)
        self.cache = CacheManager(self.output_dir)
        self.bvbrc = BVBRCDownloader()
        self.ncbi = NCBIDownloader(email=ncbi_email) if ncbi_email else None
        self.extractor = SequenceExtractor()
        
        # Statistics
        self.stats = {
            'total': 0,
            'success': 0,
            'failed': 0,
            'skipped': 0
        }

        logger.info(f"Parallel downloads: {'enabled' if parallel_downloads else 'disabled'}")
        if parallel_downloads:
            logger.info(f"Max workers: {self.max_workers}")


    def _calculate_optimal_workers(self, total_records: int) -> int:
        # assume is complete
    
    def run(self):
        """Run the genome fetcher"""
# assume is complete        
        # Parse input file
        try:
            records = self.parser.parse(node_column=self.node_column)
        except Exception as e:
            logger.error(f"Failed to parse input file: {e}")
            return False
        
        self.stats['total'] = len(records)
        logger.info(f"Processing {self.stats['total']} records")
        
        # Decide whether to use parallel processing
        if self.parallel_downloads and len(records) > 1:
            optimal_workers = self._calculate_optimal_workers(len(records))
            logger.info(f"Using parallel processing with {optimal_workers} workers")
            logger.info("")
            success = self._run_parallel(records, optimal_workers)
        else:
            logger.info("Using sequential processing")
            logger.info("")
            success = self._run_sequential(records)
        
        # Print summary
        self.print_summary()

        self.generate_detailed_report()
        
        return self.stats['failed'] == 0
    
    def _run_sequential(self, records: List[Dict]) -> bool:
        """Sequential processing (original behavior)"""
        for idx, record in enumerate(records, 1):
            logger.info(f"[{idx}/{self.stats['total']}] Processing record {idx}")
            success = self.process_record(record, idx)
            
            if success:
                self.stats['success'] += 1
            else:
                self.stats['failed'] += 1
            
            logger.info("")
        
        return True
    
    def _run_parallel(self, records: List[Dict], num_workers: int) -> bool:
        """
        Parallel processing of multiple genomes
        
        Args:
            records: List of records to process
            num_workers: Number of parallel workers
        """
        from threading import Lock
        
        # Thread-safe counter for progress
        progress_lock = Lock()
        completed = {'count': 0}
        
        def process_with_progress(record_tuple):
            """Wrapper to track progress"""
            idx, record = record_tuple
            
            # Process the record
            success = self.process_record(record, idx)
            
            # Update stats (thread-safe)
            with progress_lock:
                completed['count'] += 1
                if success:
                    self.stats['success'] += 1
                else:
                    self.stats['failed'] += 1
                
                logger.info(f"Progress: {completed['count']}/{self.stats['total']} completed")
                logger.info("")
            
            return success
        
        # Create list of (index, record) tuples
        indexed_records = [(idx, record) for idx, record in enumerate(records, 1)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            # Submit all tasks
            futures = [executor.submit(process_with_progress, record) 
                      for record in indexed_records]
            
            # Wait for all to complete
            concurrent.futures.wait(futures)
        
        return True
    
    def _get_files_info(self, terminal_dir: Path) -> Dict:
        """Get information about downloaded files"""

        files_info = {}
        for filename in REQUIRED_FILES + OPTIONAL_FILES:
            filepath = terminal_dir / filename
            if filepath.exists():
                size = filepath.stat().st_size
                status = "success" if size > 50 else "failed"
            else:
                size = 0
                status = "missing"
            
            files_info[filename] = {
                "size_bytes": size,
                "status": status
            }
        
        return files_info


    def process_record(self, record: dict, record_num: int) -> bool:
        """
        Process a single record
        
        Args:
            record: Record dictionary from TSV
            record_num: Record number for logging
        
        Returns:
            True if successful
        """


        start_time = time.time()  # Track processing time
        files_info = {}
        errors_list = []
        warnings_list = []


        # Get key information
        node_value = record.get(self.node_column, '')
        accession = record.get('Accesion', '')
        species = record.get('Species_strain', '')
        acc_type = record.get('_accession_type', 'unknown')
        genome_id = record.get('_genome_id')
        
        # Create terminal directory name
        terminal_name = f"terminal_{node_value}"
        terminal_dir = self.output_dir / terminal_name
        terminal_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize terminal logger
        log_file = terminal_dir / "download.log"
        terminal_log = DownloadLogger(log_file)
        
        terminal_log.info(f"{'='*60}")
        terminal_log.info(f"Processing Record {record_num}")
        terminal_log.info(f"{'='*60}")
        terminal_log.info(f"Node: {node_value}")
        terminal_log.info(f"Accession: {accession}")
        terminal_log.info(f"Species: {species}")
        terminal_log.info(f"Accession Type: {acc_type}")
        if genome_id:
            terminal_log.info(f"Genome ID: {genome_id}")
        terminal_log.info("")
        
        logger.info(f"  Node: {node_value}")
        logger.info(f"  Accession: {accession}")
        logger.info(f"  Type: {acc_type}")
        
        # Check if already complete
        if self.skip_existing and self.cache.is_complete(terminal_name, REQUIRED_FILES):
            logger.info(f"  ✓ Already complete - skipping")
            terminal_log.info("Terminal already complete - skipping")
            self.stats['skipped'] += 1
            return True
        
        # Process based on accession type
        try:
            if acc_type == 'bvbrc_fig':
                success = self.process_bvbrc_record(
                    accession, genome_id, terminal_dir, terminal_log
                )
            elif acc_type == 'ncbi_protein':
                # Create sequence_id from terminal name and accession
                sequence_id = f"{terminal_name}_{accession}"
                success = self.process_ncbi_protein(
                    accession, terminal_dir, terminal_log, sequence_id
                )
            elif acc_type == 'ncbi_refseq':
                terminal_log.warning("NCBI RefSeq genomes not yet supported")
                logger.warning(f"  ⚠ NCBI RefSeq not yet supported: {accession}")
                success = False
            else:
                terminal_log.error(f"Unknown accession type: {acc_type}")
                logger.error(f"  ✗ Unknown accession type: {acc_type}")
                success = False
            

            if success:
                processing_time = time.time() - start_time
                
                # Gather file info
                files_info = self._get_files_info(terminal_dir)
                
                logger.info(f"  ✓ Successfully processed ({processing_time:.1f}s)")
                terminal_log.success(f"Processing completed successfully in {processing_time:.1f}s")
                
                self.cache.mark_complete(
                    terminal_name, 
                    metadata={
                        'accession': accession,
                        'species': species,
                        'type': acc_type,
                        'genome_id': genome_id
                    },
                    files_info=files_info,
                    processing_time=processing_time,
                    download_source='bvbrc_web_api' if acc_type == 'bvbrc_fig' else 'ncbi',
                    warnings=warnings_list
                )

            else:
                processing_time = time.time() - start_time
                files_info = self._get_files_info(terminal_dir)
                
                logger.info(f"  ✗ Processing failed ({processing_time:.1f}s)")
                terminal_log.error("Processing failed")
                
                self.cache.mark_failed(
                    terminal_name, 
                    error="Download or extraction failed",
                    metadata={'accession': accession, 'species': species, 'type': acc_type},
                    files_info=files_info,
                    processing_time=processing_time,
                    error_type="ProcessingError"
                )

            return success
            
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            terminal_log.error(f"Unexpected error: {e}")
            import traceback
            terminal_log.error(f"Traceback: {traceback.format_exc()}")
            self.cache.mark_failed(terminal_name, str(e))
            return False
    
    def process_bvbrc_record(
        self,
        accession: str,
        genome_id: str,
        terminal_dir: Path,
        terminal_log: DownloadLogger
    ) -> bool:
        """Process a BV-BRC (fig|) accession"""
        terminal_log.info("Processing BV-BRC accession")
        
        # Download genome files from BV-BRC
        download_success = self.bvbrc.download_genome(
            genome_id, terminal_dir, terminal_log
        )
        
        if not download_success:
            terminal_log.error("Failed to download genome from BV-BRC")
            return False
        
        # Extract specific protein sequences
        sequence_id = f"{terminal_dir.name}_{accession}"
        extract_success = self.extractor.extract_protein(
            accession, terminal_dir, sequence_id, terminal_log
        )
        
        if not extract_success:
            terminal_log.error("Failed to extract protein sequences")
            return False
        
        return True
    
    def process_ncbi_protein(
        self,
        accession: str,
        terminal_dir: Path,
        terminal_log: DownloadLogger,
        sequence_id: str = None
    ) -> bool:
        """
        Process an NCBI protein accession - ENHANCED VERSION
        
        Args:
            accession: NCBI protein accession (e.g., "WP_000111473.1")
            terminal_dir: Terminal output directory
            terminal_log: Logger for this terminal
            sequence_id: Optional custom sequence ID
            
        Returns:
            True if successful
        """
        terminal_log.info("Processing NCBI protein accession")
        
        if not self.ncbi:
            terminal_log.error("NCBI downloader not initialized (missing email)")
            logger.error("  ✗ NCBI email not provided")
            return False
        
        if sequence_id is None:
            sequence_id = accession
        
        # Step 1: Resolve protein to genome
        terminal_log.info("Resolving protein to source genome...")
        genome_info = self.ncbi.resolve_protein_to_genome(accession, terminal_log)
        
        if not genome_info:
            terminal_log.error("Failed to resolve protein to genome")
            return False
        
        genome_acc, taxid = genome_info
        terminal_log.info(f"Resolved to genome {genome_acc} (taxid: {taxid})")
        
        # Step 2: Try to download genome from NCBI
        terminal_log.info("Downloading genome from NCBI...")
        genome_success = self.ncbi.download_genome_files(
            genome_acc, 
            terminal_dir, 
            terminal_log
        )
        
        if not genome_success:
            terminal_log.warning("NCBI genome download failed")
            
            # Step 3: Fallback to BV-BRC using taxid
            if self.try_bvbrc_fallback and taxid:
                terminal_log.info(f"Attempting BV-BRC fallback with taxid {taxid}...")
                bvbrc_success = self._try_bvbrc_by_taxid(
                    taxid, 
                    terminal_dir, 
                    terminal_log
                )
                if not bvbrc_success:
                    terminal_log.error("BV-BRC fallback also failed")
                    return False
                # Update genome_success flag for later checks
                genome_success = True
            else:
                return False
        
        # Step 4: Download protein sequences directly
        terminal_log.info("Downloading protein sequences...")
        protein_success = self.ncbi.download_protein_sequences(
            accession,
            terminal_dir,
            terminal_log,
            sequence_id=sequence_id
        )
        
        if not protein_success:
            terminal_log.warning("Direct protein download failed")
            
            # Try to extract from genome if we have it
            if genome_success:
                terminal_log.info("Attempting to extract protein from genome files...")
                extract_success = self._extract_protein_from_genome(
                    accession,
                    terminal_dir,
                    terminal_log,
                    sequence_id
                )
                if extract_success:
                    terminal_log.success("Successfully extracted protein from genome")
                    return True
                else:
                    terminal_log.warning("Could not extract protein from genome")
            
            # As last resort, return partial success if we have genome
            if genome_success:
                terminal_log.warning("Genome downloaded but protein extraction incomplete")
                terminal_log.info("You can extract the protein manually from genome files")
                return True
            else:
                terminal_log.error("Failed to download both genome and protein")
                return False
        
        terminal_log.success("Successfully downloaded protein sequences")
        return True
    
    def _try_bvbrc_by_taxid(
        self,
        taxid: str,
        terminal_dir: Path,
        terminal_log: DownloadLogger
    ) -> bool:
        """
        Try to find and download genome from BV-BRC using taxonomy ID
        
        Args:
            taxid: NCBI taxonomy ID
            terminal_dir: Output directory
            terminal_log: Logger
            
        Returns:
            True if successful
        """
        terminal_log.info(f"Searching BV-BRC for genomes with taxid {taxid}...")
        
        try:
            # Query BV-BRC API for genomes with this taxid
            import requests
            url = f"https://www.bv-brc.org/api/genome/?eq(taxon_id,{taxid})&select(genome_id,genome_name)&limit(10)"
            
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                data = response.json()
                
                if data and len(data) > 0:
                    # Use the first genome found
                    genome_id = data[0].get('genome_id')
                    genome_name = data[0].get('genome_name', 'Unknown')
                    
                    terminal_log.info(f"Found BV-BRC genome: {genome_id} ({genome_name})")
                    
                    # Download from BV-BRC
                    success = self.bvbrc.download_genome(
                        genome_id,
                        terminal_dir,
                        terminal_log
                    )
                    
                    if success:
                        terminal_log.success(f"Successfully downloaded genome from BV-BRC")
                        return True
                    else:
                        terminal_log.warning("BV-BRC genome download failed")
                        return False
                else:
                    terminal_log.warning(f"No genomes found in BV-BRC for taxid {taxid}")
                    return False
            else:
                terminal_log.warning(f"BV-BRC API query failed (status {response.status_code})")
                return False
                
        except Exception as e:
            terminal_log.error(f"Error querying BV-BRC: {e}")
            logger.error(f"BV-BRC fallback error: {e}")
            return False
    
    def _extract_protein_from_genome(
        self,
        protein_accession: str,
        terminal_dir: Path,
        terminal_log: DownloadLogger,
        sequence_id: str
    ) -> bool:
        """
        Try to extract protein from downloaded genome files
        
        Args:
            protein_accession: Protein accession to find
            terminal_dir: Directory with genome files
            terminal_log: Logger
            sequence_id: Output sequence ID
            
        Returns:
            True if successful
        """
        try:
            # Check for genome annotation files
            genome_faa = terminal_dir / "genome.faa"
            genome_ffn = terminal_dir / "genome.ffn"
            
            # If we have FAA/FFN files (from BV-BRC), try to extract
            if genome_faa.exists() or genome_ffn.exists():
                terminal_log.info("Found genome annotation files, searching for protein...")
                
                success = self.extractor.extract_protein(
                    protein_accession,
                    terminal_dir,
                    sequence_id,
                    terminal_log
                )
                
                return success
            
            # If we have GenBank files, we could parse them (more complex)
            genome_gbk = terminal_dir / "genome.gbk"
            genome_gbff = terminal_dir / "genome.gbff.gz"
            
            if genome_gbk.exists() or genome_gbff.exists():
                terminal_log.info("GenBank files available but parsing not implemented")
                terminal_log.info("Consider using BioPython to extract from GenBank")
            
            return False
            
        except Exception as e:
            terminal_log.error(f"Error extracting protein: {e}")
            return False
    
    def print_summary(self):
        # assume is complete

    def generate_detailed_report(self, output_file: Path = None):
        """Generate detailed report from cache"""
        # assume is complete

    def _write_failure_details(self, f, terminal, info):
        """Write details for a failed download"""
    # assume is complete

    def _write_partial_details(self, f, terminal, info):
        """Write details for a partial download"""
    # assume is complete


def main():
    """Main entry point"""
    # assume is complete
    
    # Use NCBI_EMAIL from config if not provided via command line
    ncbi_email = args.ncbi_email if args.ncbi_email else NCBI_EMAIL
    
    if not ncbi_email:
        logger.warning("No NCBI email provided - NCBI protein downloads will be skipped")
        logger.info("Set NCBI_EMAIL in config or use --ncbi-email argument")
    
    # Create fetcher and run
    fetcher = GenomeFetcher(
        input_file=input_file,
        output_dir=Path(args.output),
        node_column=args.node_column,
        ncbi_email=ncbi_email,
        skip_existing=args.skip_existing,
        try_bvbrc_fallback=args.try_bvbrc_fallback,
        parallel_downloads=True,  # ← Enable parallel
        max_workers=None   # ← Use 20 workers (or None for auto)
    )
    
    success = fetcher.run()
    

if __name__ == '__main__':
    main()