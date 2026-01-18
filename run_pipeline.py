#!/usr/bin/env python3
"""
Retron Pipeline Orchestrator

This script orchestrates the complete retron analysis pipeline for bacterial genomes.
It processes genome directories through multiple analysis tools and integrates the results.

Usage:
    python retron_pipeline_orchestrator.py -i <input_genome_dir> -o <output_base_dir> [options]
    
    Or for batch processing:
    python retron_pipeline_orchestrator.py -b <batch_input_dir> -o <output_base_dir> [options]
"""

import os
import sys
import argparse
import subprocess 
from pathlib import Path
from typing import Optional, List
import json
from datetime import datetime
import logging
import multiprocessing
from datetime import datetime
from pathlib import Path
from difflib import SequenceMatcher
from typing import Dict, List, Optional
import json

def get_optimal_threads(requested_threads=None, utilization=0.8):
    """
    Get optimal thread count for sequential tool execution.
    
    Args:
        requested_threads: User-requested thread count (None for auto)
        utilization: Fraction of available CPUs to use (default: 0.8 = 80%)
    
    Returns:
        int: Number of threads to use
    """
    available_cpus = multiprocessing.cpu_count()
    
    if requested_threads is None:
        # Auto-detect: use 80% of available CPUs
        optimal = max(1, int(available_cpus * utilization))
    else:
        # Use requested, but cap at available
        optimal = min(requested_threads, available_cpus)
    
    return optimal

def calculate_similarity(seq1: str, seq2: str) -> float:
    """Calculate sequence similarity ratio between two sequences."""
    return SequenceMatcher(None, seq1.upper(), seq2.upper()).ratio()

# Try to import the RNA annotation module
try:  
    sys.path.append('/ibex/user/rioszemm/the-retron-project/pipeline')
    # Change this import to use the new class
    from one_sequence_bprna_annotation import annotate_rna_sequence
    BPRNA_AVAILABLE = True
except ImportError:
    BPRNA_AVAILABLE = False
    print("Warning: rna_annotator module not found. bpRNA annotation will be skipped.")


class RetronPipelineOrchestrator:
    """Orchestrates the retron analysis pipeline."""
    
    def __init__(self, args):
        self.args = args
        self.setup_logging()

        self.genome_used_by_tools = None  # Add this
        
        # Define script paths (customize these to your environment)
        self.script_paths = {
            'step1': '/ibex/user/rioszemm/the-retron-project/pipeline_for_git/TEST_SCRIPT_31_nov_v5.py',
            'extract_myrt': '/ibex/user/rioszemm/the-retron-project/pipeline_for_git/extract_mrt_output_2_dic.py',
            'extract_padloc': '/ibex/user/rioszemm/the-retron-project/pipeline_for_git/extract_padloc_outputs_1_dic.py',
            'parse_defensefinder': '/ibex/user/rioszemm/the-retron-project/pipeline_for_git/parse_defense_finder_1_dic.py',
            'integrate': '/ibex/user/rioszemm/the-retron-project/pipeline_for_git/rt_integration_corrected.py',
            'annotate_ncrna': '/ibex/user/rioszemm/the-retron-project/pipeline_for_git/ANNOTATE_bprna.py'
        }
        

        available_cpus = multiprocessing.cpu_count()
        self.logger.info(f"System CPUs detected: {available_cpus}")
        self.logger.info(f"Using {self.args.threads} threads for tool execution")
        
        if self.args.threads == available_cpus:
            self.logger.warning(
                f"Using all {available_cpus} CPUs. System may become unresponsive. "
                f"Consider using --threads {int(available_cpus * 0.8)} instead."
            )


    def setup_logging(self):
        """Setup logging configuration."""
        log_level = logging.DEBUG if self.args.verbose else logging.INFO
        log_format = '%(asctime)s - %(levelname)s - %(message)s'
        
        # Setup file logging if specified
        if self.args.log:
            logging.basicConfig(
                level=log_level,
                format=log_format,
                handlers=[
                    logging.FileHandler(self.args.log),
                    logging.StreamHandler(sys.stdout)
                ]
            )
        else:
            logging.basicConfig(
                level=log_level,
                format=log_format,
                handlers=[logging.StreamHandler(sys.stdout)]
            )
        
        self.logger = logging.getLogger(__name__)
    
    def validate_genome_directory(self, genome_dir: Path) -> bool:
        """
        Validate that the genome directory contains required files.
        
        Args:
            genome_dir: Path to genome directory
            
        Returns:
            bool: True if valid, False otherwise
        """
        # required_file = genome_dir / "genome.fasta"
        required_files = list(genome_dir.glob("genome.fasta*"))

        
        # if not required_file.exists():
        #     self.logger.error(f"Required file not found: {required_file}")
        #     return False

        if not required_files:
            self.logger.error(f"Required file not found: genome.fasta or compressed equivalent in {genome_dir}")
            return False

        
        self.logger.info(f"Validated genome directory: {genome_dir}")
        return True
    
    def check_ground_truth(self, genome_dir: Path) -> bool:
        """
        Check if this is ground truth data (has ncRNA_DNA.fasta).
        
        Args:
            genome_dir: Path to genome directory
            
        Returns:
            bool: True if ground truth data exists
        """
        gt_file = genome_dir / "ncRNA_DNA.fasta"
        is_gt = gt_file.exists()
        
        if is_gt:
            self.logger.info(f"Ground truth data detected in {genome_dir}")
        
        return is_gt
    
    def run_command(self, cmd: List[str], step_name: str) -> bool:
        """
        Run a command and handle errors.
        
        Args:
            cmd: Command to run as list of strings
            step_name: Name of the step for logging
            
        Returns:
            bool: True if successful, False otherwise
        """
        self.logger.info(f"Starting {step_name}...")
        self.logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            if result.stdout:
                self.logger.debug(f"STDOUT: {result.stdout}")
            
            self.logger.info(f"✓ {step_name} completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"✗ {step_name} failed with return code {e.returncode}")
            self.logger.error(f"STDERR: {e.stderr}")
            if e.stdout:
                self.logger.error(f"STDOUT: {e.stdout}")
            return False
        except Exception as e:
            self.logger.error(f"✗ {step_name} failed with exception: {str(e)}")
            return False
    
    def step0b_annotate_gt_ncrna(self, genome_dir: Path, output_dir: Path) -> bool:
        """
        STEP 0b: Annotate ground truth ncRNA with bpRNA/RNAfold.
        Only runs if ground_truth_metadata.json exists.
        """
        gt_metadata_path = output_dir / "ground_truth_metadata.json"
        
        if not gt_metadata_path.exists():
            self.logger.info("No ground truth metadata found, skipping ncRNA annotation")
            return True
        
        # Check if ncRNA exists in the metadata
        try:
            with open(gt_metadata_path, 'r') as f:
                metadata = json.load(f)
            
            if 'ncRNA' not in metadata or metadata['ncRNA'] is None:
                self.logger.info("No ncRNA in ground truth metadata, skipping annotation")
                return True
        except Exception as e:
            self.logger.warning(f"Could not check ground truth metadata: {e}")
            return True
        
        # Import the annotation module
        try:
            import sys
            pipeline_path = '/ibex/user/rioszemm/the-retron-project/pipeline'
            if pipeline_path not in sys.path:
                sys.path.insert(0, pipeline_path)
            
            from one_sequence_bprna_annotation import annotate_rna_sequence
            self.logger.info("Successfully imported RNA annotator")
            
        except ImportError as e:
            self.logger.warning(f"Could not import rna_annotator: {e}")
            return True
        
        # SET THE BPRNA SCRIPT PATH HERE
        bprna_script_path = "/ibex/user/rioszemm/the-retron-project/src/bpRNA/bpRNA.pl"  # <-- UPDATE THIS!
        
        # Check if bpRNA script exists
        if not os.path.exists(bprna_script_path):
            self.logger.warning(f"bpRNA script not found at {bprna_script_path}")
            self.logger.warning("Will proceed with RNAfold-only annotation")
            bprna_script_path = None
        
        self.logger.info("Annotating ground truth ncRNA structure...")
        
        # Annotate single or multiple ncRNA sequences
        ncrna = metadata['ncRNA']
        
        try:
            if isinstance(ncrna, dict):
                # Single sequence
                sequence = ncrna.get('sequence')
                seq_id = ncrna.get('id', 'gt_ncrna')
                
                if sequence:
                    self.logger.info(f"Annotating sequence: {seq_id} ({len(sequence)} nt)")
                    
                    # PASS THE BPRNA SCRIPT PATH HERE
                    result = annotate_rna_sequence(
                        sequence, 
                        seq_id,
                        bprna_script_path=bprna_script_path  # <-- THIS IS THE KEY!
                    )
                    
                    ncrna['bprna_annotation'] = result.get('bprna_annotation', {})
                    
                    status = ncrna['bprna_annotation'].get('annotation_status', 'unknown')
                    self.logger.info(f"✓ Annotated ncRNA: {seq_id} (status: {status})")
            
            elif isinstance(ncrna, list):
                # Multiple sequences
                for i, ncrna_item in enumerate(ncrna):
                    sequence = ncrna_item.get('sequence')
                    seq_id = ncrna_item.get('id', f'gt_ncrna_{i}')
                    
                    if sequence:
                        self.logger.info(f"Annotating sequence {i+1}: {seq_id} ({len(sequence)} nt)")
                        
                        # PASS THE BPRNA SCRIPT PATH HERE TOO
                        result = annotate_rna_sequence(
                            sequence, 
                            seq_id,
                            bprna_script_path=bprna_script_path  # <-- AND HERE!
                        )
                        
                        ncrna_item['bprna_annotation'] = result.get('bprna_annotation', {})
                        
                        status = ncrna_item['bprna_annotation'].get('annotation_status', 'unknown')
                        self.logger.info(f"✓ Annotated ncRNA {i+1}: {seq_id} (status: {status})")
            
            # Save updated metadata
            metadata['ncRNA'] = ncrna
            with open(gt_metadata_path, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            self.logger.info("✓ Ground truth ncRNA annotation completed")
            return True
            
        except Exception as e:
            self.logger.error(f"Error annotating ground truth ncRNA: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False
        
    # def step1_run_analysis_tools(self, genome_fasta: Path, output_dir: Path) -> bool:
    #     """
    #     STEP 1: Run initial analysis with TEST_SCRIPT_31_nov_v5.py
    #     This runs Prodigal, MyRT, PADLOC, DefenseFinder, and optionally Infernal.
    #     """
    #     cmd = [
    #         'python',
    #         self.script_paths['step1'],
    #         '-i', str(genome_fasta),
    #         '-o', str(output_dir)
    #     ]
        
    #     if self.args.run_infernal:
    #         cmd.append('--run-infernal')
        
    #     if self.args.threads:
    #         cmd.extend(['--threads', str(self.args.threads)])
        
    #     if self.args.evalue:
    #         cmd.extend(['--evalue', str(self.args.evalue)])
        
    #     return self.run_command(cmd, "STEP 1: Run analysis tools")


    # def step1_run_analysis_tools(self, genome_fasta: Path, output_dir: Path) -> bool:
    #     """
    #     STEP 1: Run initial analysis with TEST_SCRIPT_31_nov_v5.py
    #     This runs Prodigal, MyRT, PADLOC, DefenseFinder, and optionally Infernal.
    #     """
    #     self.logger.info(f"Running analysis tools with {self.args.threads} threads")
        
    #     cmd = [
    #         'python',
    #         self.script_paths['step1'],
    #         '-i', str(genome_fasta),
    #         '-o', str(output_dir),
    #         '--threads', str(self.args.threads)  # Pass threads to the tool
    #     ]
        
    #     if self.args.run_infernal:
    #         cmd.append('--run-infernal')
        
    #     if self.args.evalue:
    #         cmd.extend(['--evalue', str(self.args.evalue)])
        
    #     return self.run_command(cmd, "STEP 1: Run analysis tools")


    def step1_run_analysis_tools(self, genome_fasta: Path, output_dir: Path) -> bool:
        """
        STEP 1: Run initial analysis with TEST_SCRIPT_31_nov_v5.py
        """
        self.logger.info(f"Running analysis tools with {self.args.threads} threads")
        
        cmd = [
            'python',
            self.script_paths['step1'],
            '-i', str(genome_fasta),
            '-o', str(output_dir),
            '--threads', str(self.args.threads)
        ]
        
        if self.args.run_infernal:
            cmd.append('--run-infernal')
        if self.args.evalue:
            cmd.extend(['--evalue', str(self.args.evalue)])
        
        success = self.run_command(cmd, "STEP 1: Run analysis tools")
        
        # After Step1 completes, determine which genome was used
        if success:
            prodigal_dir = output_dir / "prodigal_results"
            
            # Priority 1: Sanitized genome (created by Prodigal if sanitization was needed)
            sanitized_genome = prodigal_dir / "genome_sanitized.fasta"
            if sanitized_genome.exists():
                self.genome_used_by_tools = sanitized_genome
                self.logger.info(f"Tools used sanitized genome: {sanitized_genome}")
                return success
            
            # Priority 2: Decompressed genome in INPUT directory (where Prodigal decompresses .gz files)
            if str(genome_fasta).endswith('.gz'):
                decompressed = genome_fasta.with_suffix('')  # Removes .gz -> genome.fasta
                if decompressed.exists():
                    self.genome_used_by_tools = decompressed
                    self.logger.info(f"Tools used decompressed genome: {decompressed}")
                    return success
                else:
                    self.logger.error(f"Expected decompressed genome not found: {decompressed}")
                    self.logger.error("Prodigal should have decompressed the .gz file")
                    return False
            
            # Priority 3: Original genome (already uncompressed)
            else:
                if genome_fasta.exists():
                    self.genome_used_by_tools = genome_fasta
                    self.logger.info(f"Tools used original genome: {genome_fasta}")
                    return success
                else:
                    self.logger.error(f"Original genome not found: {genome_fasta}")
                    return False
        
        return success
    
    
    def step2_extract_myrt(self, output_dir: Path) -> bool:
        """
        STEP 2a: Extract MyRT output to JSON.
        """
        myrt_output = output_dir / "myrt_output"
        
        if not myrt_output.exists():
            self.logger.warning(f"MyRT output directory not found: {myrt_output}")
            return False
        
        cmd = [
            'python',
            self.script_paths['extract_myrt'],
            str(myrt_output)
        ]
        
        return self.run_command(cmd, "STEP 2a: Extract MyRT results")


    def step2_extract_padloc(self, output_dir: Path) -> bool:
        """
        STEP 2b: Extract PADLOC output to JSON.
        """
        cmd = [
            'python',
            self.script_paths['extract_padloc'],
            str(output_dir),
            '--retron-only'  # Add this flag
        ]
        
        return self.run_command(cmd, "STEP 2b: Extract PADLOC results")
    

    def step2_parse_defensefinder(self, output_dir: Path) -> bool:
        """
        STEP 2c: Parse DefenseFinder output to JSON.
        """
        cmd = [
            'python',
            self.script_paths['parse_defensefinder'],
            '-i', str(output_dir),
            '--system-type', 'Retron'  # Add this option
        ]
        
        return self.run_command(cmd, "STEP 2c: Parse DefenseFinder results")
    
    # def step3_integrate_results(self, output_dir: Path, genome_fasta: Path) -> bool:
    #     """
    #     STEP 3: Integrate outputs from all tools.
    #     """
    #     integrated_output = output_dir / "integrated_results"
    #     integrated_output.mkdir(exist_ok=True)
        
    #     cmd = [
    #         'python',
    #         self.script_paths['integrate'],
    #         '--input-dir', str(output_dir),
    #         '--fasta', str(genome_fasta),
    #         '--output', str(integrated_output)
    #     ]
        
    #     return self.run_command(cmd, "STEP 3: Integrate tool results")


    # def step3_integrate_results(self, output_dir: Path, genome_fasta: Path) -> bool:
    #     """
    #     STEP 3: Integrate outputs from all tools.
    #     """
    #     integrated_output = output_dir / "integrated_results"
    #     integrated_output.mkdir(exist_ok=True)
        
    #     # ✅ Check for sanitized genome in Prodigal directory
    #     sanitized_genome = output_dir / "prodigal_results" / "genome_sanitized.fasta"
        
    #     if sanitized_genome.exists():
    #         genome_to_use = sanitized_genome
    #         self.logger.info(f"Using sanitized genome: {sanitized_genome}")
    #     else:
    #         genome_to_use = genome_fasta
    #         self.logger.info(f"Using original genome: {genome_fasta}")
        
    #     cmd = [
    #         'python',
    #         self.script_paths['integrate'],
    #         '--input-dir', str(output_dir),
    #         '--fasta', str(genome_to_use),
    #         '--output', str(integrated_output)
    #     ]
    #     return self.run_command(cmd, "STEP 3: Integrate tool results")

    # def step3_integrate_results(self, output_dir: Path, genome_fasta: Path) -> bool:
    #     """
    #     STEP 3: Integrate outputs from all tools.
    #     """
    #     integrated_output = output_dir / "integrated_results"
    #     integrated_output.mkdir(exist_ok=True)
        
    #     # ✅ USE SANITIZED GENOME
    #     sanitized_genome = output_dir / "prodigal_results" / "genome_sanitized.fasta"
    #     genome_to_use = sanitized_genome if sanitized_genome.exists() else genome_fasta
        
    #     self.logger.info(f"Integration using: {genome_to_use}")
        
    #     cmd = [
    #         'python',
    #         self.script_paths['integrate'],
    #         '--input-dir', str(output_dir),
    #         '--fasta', str(genome_to_use),  # ← Use sanitized genome!
    #         '--output', str(integrated_output)
    #     ]
    #     return self.run_command(cmd, "STEP 3: Integrate tool results")



    # def step3_integrate_results(self, output_dir: Path) -> bool:
    #     """Integration using the genome that tools actually used"""
        
    #     if not hasattr(self, 'genome_used_by_tools') or not self.genome_used_by_tools:
    #         self.logger.error("Genome path not available - Step1 must run first")
    #         return False
        
    #     genome_fasta = self.genome_used_by_tools
        
    #     # Validate it exists and isn't compressed
    #     if not genome_fasta.exists():
    #         self.logger.error(f"Genome file not found: {genome_fasta}")
    #         return False
        
    #     if str(genome_fasta).endswith('.gz'):
    #         self.logger.error(f"Genome is still compressed: {genome_fasta}")
    #         self.logger.error("Integration requires uncompressed genome")
    #         return False
        
    #     integrated_output = output_dir / "integrated_results"
    #     integrated_output.mkdir(exist_ok=True)
        
    #     self.logger.info(f"Integration using: {genome_fasta}")
        
    #     cmd = [
    #         'python',
    #         self.script_paths['integrate'],
    #         '--input-dir', str(output_dir),
    #         '--fasta', str(genome_fasta),
    #         '--output', str(integrated_output)
    #     ]
    #     return self.run_command(cmd, "STEP 3: Integrate tool results")



    def step3_integrate_results(self, output_dir: Path) -> bool:
        """Integration using the genome that tools actually used"""
        
        if not hasattr(self, 'genome_used_by_tools') or not self.genome_used_by_tools:
            self.logger.error("Genome path not available - Step1 must run first")
            return False
        
        genome_fasta = self.genome_used_by_tools
        
        # Validate it exists and isn't compressed
        if not genome_fasta.exists():
            self.logger.error(f"Genome file not found: {genome_fasta}")
            return False
        
        if str(genome_fasta).endswith('.gz'):
            self.logger.error(f"Genome is still compressed: {genome_fasta}")
            self.logger.error("Integration requires uncompressed genome")
            return False
        
        integrated_output = output_dir / "integrated_results"
        integrated_output.mkdir(exist_ok=True)
        
        json_file = integrated_output / "tool_integrated_results.json"
        
        # Checkpoint: skip if valid integration already exists
        if json_file.exists() and json_file.stat().st_size > 0:
            try:
                with open(json_file) as f:
                    data = json.load(f)
                if isinstance(data, list) and len(data) > 0:
                    self.logger.info(f"✓ Integration already complete ({len(data)} systems)")
                    return True
            except:
                pass  # If JSON is invalid, re-run integration
        
        self.logger.info(f"Integration using: {genome_fasta}")
        
        cmd = [
            'python',
            self.script_paths['integrate'],
            '--input-dir', str(output_dir),
            '--fasta', str(genome_fasta),
            '--output', str(integrated_output)
        ]
        return self.run_command(cmd, "STEP 3: Integrate tool results")
        
    def step0_extract_ground_truth(self, genome_dir: Path, output_dir: Path) -> bool:
        """
        STEP 0: Extract ground truth metadata if available.
        Creates a JSON file with RT protein and ncRNA sequences from ground truth files.
        Also annotates ncRNA with bpRNA if available.
        """
        gt_files = {
            'ncRNA_DNA': genome_dir / "ncRNA_DNA.fasta",
            'ncRNA_RNA': genome_dir / "ncRNA_RNA.fasta",
            'protein_aa': genome_dir / "protein_aminoacid.fasta"
        }
        
        # Check if any ground truth files exist
        existing_files = {k: v for k, v in gt_files.items() if v.exists()}
        
        if not existing_files:
            self.logger.info("No ground truth files found, skipping ground truth extraction")
            return True
        
        self.logger.info(f"Found ground truth files: {list(existing_files.keys())}")
        
        # Initialize ground truth metadata
        gt_metadata = {
            'genome_id': genome_dir.name,
            'source_directory': str(genome_dir),
            'ground_truth_files': {},
            'retron_RT_protein': None,
            'ncRNA': None
        }
        
        try:
            # Extract protein sequence
            if 'protein_aa' in existing_files:
                protein_file = existing_files['protein_aa']
                gt_metadata['ground_truth_files']['protein_aminoacid'] = str(protein_file)
                
                sequences = self._parse_fasta(protein_file)
                if sequences:
                    # Take the first sequence (or all if multiple)
                    if len(sequences) == 1:
                        gt_metadata['retron_RT_protein'] = {
                            'id': sequences[0]['id'],
                            'sequence': sequences[0]['sequence'],
                            'length': len(sequences[0]['sequence'])
                        }
                    else:
                        # Multiple sequences
                        gt_metadata['retron_RT_protein'] = [
                            {
                                'id': seq['id'],
                                'sequence': seq['sequence'],
                                'length': len(seq['sequence'])
                            }
                            for seq in sequences
                        ]
                    self.logger.info(f"Extracted {len(sequences)} protein sequence(s)")
            
            # Extract ncRNA sequence (prioritize RNA over DNA)
            ncrna_sequences = None
            ncrna_type = None
            
            if 'ncRNA_RNA' in existing_files:
                ncrna_file = existing_files['ncRNA_RNA']
                gt_metadata['ground_truth_files']['ncRNA_RNA'] = str(ncrna_file)
                ncrna_sequences = self._parse_fasta(ncrna_file)
                ncrna_type = 'RNA'
                self.logger.info(f"Extracted {len(ncrna_sequences)} ncRNA sequence(s) (RNA)")
            
            elif 'ncRNA_DNA' in existing_files:
                ncrna_file = existing_files['ncRNA_DNA']
                gt_metadata['ground_truth_files']['ncRNA_DNA'] = str(ncrna_file)
                ncrna_sequences = self._parse_fasta(ncrna_file)
                ncrna_type = 'DNA'
                self.logger.info(f"Extracted {len(ncrna_sequences)} ncRNA sequence(s) (DNA)")
            
            # Process ncRNA sequences and add bpRNA annotation
            if ncrna_sequences:
                if len(ncrna_sequences) == 1:
                    # Single sequence
                    seq_data = {
                        'id': ncrna_sequences[0]['id'],
                        'sequence': ncrna_sequences[0]['sequence'],
                        'length': len(ncrna_sequences[0]['sequence']),
                        'type': ncrna_type
                    }
                    
                    # Add bpRNA annotation if available
                    if BPRNA_AVAILABLE and self.args.annotate_ncrna:
                        self.logger.info("Running bpRNA annotation on ncRNA sequence...")
                        bprna_result = self._annotate_ncrna_with_bprna(
                            ncrna_sequences[0]['sequence'],
                            ncrna_sequences[0]['id']
                        )
                        if bprna_result:
                            seq_data['bprna_annotation'] = bprna_result['bprna_annotation']
                            self.logger.info(f"✓ bpRNA annotation completed: {bprna_result['bprna_annotation']['annotation_status']}")
                    
                    gt_metadata['ncRNA'] = seq_data
                    
                else:
                    # Multiple sequences
                    ncrna_list = []
                    for seq in ncrna_sequences:
                        seq_data = {
                            'id': seq['id'],
                            'sequence': seq['sequence'],
                            'length': len(seq['sequence']),
                            'type': ncrna_type
                        }
                        
                        # Add bpRNA annotation if available
                        if BPRNA_AVAILABLE and self.args.annotate_ncrna:
                            self.logger.info(f"Running bpRNA annotation on ncRNA sequence: {seq['id']}...")
                            bprna_result = self._annotate_ncrna_with_bprna(
                                seq['sequence'],
                                seq['id']
                            )
                            if bprna_result:
                                seq_data['bprna_annotation'] = bprna_result['bprna_annotation']
                                self.logger.info(f"✓ bpRNA annotation completed: {bprna_result['bprna_annotation']['annotation_status']}")
                        
                        ncrna_list.append(seq_data)
                    
                    gt_metadata['ncRNA'] = ncrna_list
            
            # Save ground truth metadata
            gt_json_path = output_dir / "ground_truth_metadata.json"
            with open(gt_json_path, 'w') as f:
                json.dump(gt_metadata, f, indent=2)
            
            self.logger.info(f"✓ Ground truth metadata saved to: {gt_json_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error extracting ground truth metadata: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False
    
    def _annotate_ncrna_with_bprna(self, sequence: str, seq_id: str) -> Optional[dict]:
        """
        Annotate an ncRNA sequence using the bpRNA annotation script.
        
        Args:
            sequence: RNA sequence string
            seq_id: Sequence identifier
            
        Returns:
            Dictionary with bpRNA annotation results or None if failed
        """
        try:
            # Call the annotate_rna_sequence function from one_sequence_bprna_annotation.py
            result = annotate_rna_sequence(sequence, seq_id)
            return result
        except Exception as e:
            self.logger.error(f"Error running bpRNA annotation: {str(e)}")
            return None
    
    def _parse_fasta(self, fasta_file: Path) -> List[dict]:
        """
        Parse a FASTA file and return sequences.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            List of dictionaries with 'id' and 'sequence' keys
        """
        sequences = []
        current_id = None
        current_seq = []
        
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_id is not None:
                            sequences.append({
                                'id': current_id,
                                'sequence': ''.join(current_seq)
                            })
                        
                        # Start new sequence
                        current_id = line[1:].strip()
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Save last sequence
                if current_id is not None:
                    sequences.append({
                        'id': current_id,
                        'sequence': ''.join(current_seq)
                    })
        
        except Exception as e:
            self.logger.error(f"Error parsing FASTA file {fasta_file}: {str(e)}")
            return []
        
        return sequences
    
    def step4_annotate_ncrna(self, output_dir: Path) -> bool:
        """
        STEP 4: Annotate ncRNA in integrated results.
        """
        integrated_json = output_dir / "integrated_results" / "tool_integrated_results.json"
        
        if not integrated_json.exists():
            self.logger.warning(f"Integrated results JSON not found: {integrated_json}")
            return False
        
        cmd = [
            'python',
            self.script_paths['annotate_ncrna'],
            str(integrated_json)
        ]
        
        return self.run_command(cmd, "STEP 4: Annotate ncRNA structures")
    
    def safe_nested_get(data, *keys, default=None):
        """
        Safely get nested dictionary values, handling None at any level.
        
        Example:
            safe_nested_get(system, "rt_gene", "tool_metadata", "PADLOC", "ncrnas", default=[])
        """
        result = data
        for key in keys:
            if result is None:
                return default
            if isinstance(result, dict):
                result = result.get(key)
            else:
                return default
        return result if result is not None else default



    def step5_validate_results(self, genome_dir: Path, output_dir: Path) -> bool:
        """
        STEP 5: Validate pipeline results against ground truth.
        
        UPDATED LOGIC:
        1. Find system with best protein match
        2. Check ncRNA ONLY in that same system (not across all systems)
        3. Use alignment-based comparison for sequences with length differences
        4. HANDLE ncRNA-anchored systems (where rt_gene is None)
        """
        gt_metadata_path = output_dir / "ground_truth_metadata.json"
        integrated_results_path = output_dir / "integrated_results" / "tool_integrated_results_ncRNA_annotated.json"
        
        # Fallback to non-annotated version if ncRNA annotated doesn't exist
        if not integrated_results_path.exists():
            integrated_results_path = output_dir / "integrated_results" / "tool_integrated_results.json"
            self.logger.info("Using tool_integrated_results.json (no ncRNA annotation found)")
        
        # Check if ground truth exists
        if not gt_metadata_path.exists():
            self.logger.info("No ground truth metadata found, skipping validation")
            return True
        
        if not integrated_results_path.exists():
            self.logger.warning("No integrated results found for validation")
            return False
        
        try:
            self.logger.info("Starting validation against ground truth...")
            
            # Load ground truth
            with open(gt_metadata_path, 'r') as f:
                ground_truth = json.load(f)
            
            # Load integrated results
            with open(integrated_results_path, 'r') as f:
                integrated_results = json.load(f)
            
            # Handle multiple ground truth proteins/ncRNAs
            gt_proteins = ground_truth.get("retron_RT_protein", {})
            gt_ncrnas = ground_truth.get("ncRNA", {})
            
            # Normalize to list format
            if isinstance(gt_proteins, dict):
                gt_proteins = [gt_proteins]
            if isinstance(gt_ncrnas, dict):
                gt_ncrnas = [gt_ncrnas]
            elif gt_ncrnas is None:
                gt_ncrnas = []
            
            validation_results = {
                "ground_truth_id": ground_truth.get("genome_id") or (gt_proteins[0].get("id", "unknown") if gt_proteins else "unknown"),
                "total_systems_checked": len(integrated_results),
                "matches_found": [],
                "best_protein_match": None,
                "best_complete_match": None,  # NEW: Best system with BOTH protein + ncRNA
                "summary": {
                    "any_protein_match": False,
                    "any_ncrna_match": False,
                    "complete_match_systems": []
                }
            }
            
            best_protein_similarity = 0.0
            best_protein_system = None
            best_complete_similarity = 0.0  # Combined protein + ncRNA score
            best_complete_system = None
            
            # First pass: Find best protein match across all RT-anchored systems
            for system in integrated_results:
                rt_system_id = system.get("rt_system_id", "")
                anchor_type = system.get("anchor_type", "RT")  # NEW: Check anchor type
                
                # Skip ncRNA-anchored systems for protein matching
                rt_gene = system.get("rt_gene")
                if rt_gene is None:
                    continue  # Skip ncRNA-anchored systems for protein comparison
                    
                system_protein_seq = rt_gene.get("sequence", "")
                
                # Remove stop codons for fair comparison
                if system_protein_seq.endswith('*'):
                    system_protein_seq = system_protein_seq[:-1]
                
                # Find best protein match for this system
                best_gt_protein_similarity = 0.0
                matched_gt_protein_id = None
                
                for gt_protein in gt_proteins:
                    gt_protein_seq = gt_protein.get("sequence", "")
                    if gt_protein_seq.endswith('*'):
                        gt_protein_seq = gt_protein_seq[:-1]
                    
                    protein_similarity = calculate_similarity(gt_protein_seq, system_protein_seq)
                    
                    if protein_similarity > best_gt_protein_similarity:
                        best_gt_protein_similarity = protein_similarity
                        matched_gt_protein_id = gt_protein.get("id", "unknown")
                
                # Track best protein match across all systems
                if best_gt_protein_similarity > best_protein_similarity:
                    best_protein_similarity = best_gt_protein_similarity
                    best_protein_system = {
                        "rt_system_id": rt_system_id,
                        "similarity": round(best_gt_protein_similarity, 4),
                        "matched_gt_protein_id": matched_gt_protein_id
                    }
            
            # Set best protein match in results
            validation_results["best_protein_match"] = best_protein_system
            
            # Second pass: Validate each system
            for system in integrated_results:
                rt_system_id = system.get("rt_system_id", "")
                anchor_type = system.get("anchor_type", "RT")  # NEW: Check anchor type
                
                # Handle rt_gene safely - it can be None for ncRNA-anchored systems
                rt_gene = system.get("rt_gene")
                
                # Initialize match info
                match_info = {
                    "rt_system_id": rt_system_id,
                    "anchor_type": anchor_type,  # NEW: Include anchor type
                    "contig": system.get("contig", ""),
                    "protein_match": False,
                    "protein_similarity": 0.0,
                    "matched_gt_protein_id": None,
                    "ncrna_match": False,
                    "ncrna_similarity": 0.0,
                    "matched_gt_ncrna_id": None,
                    "matched_elements": [],
                    "ncrna_info": None
                }
                
                # =====================================================================
                # PROTEIN MATCHING (only for RT-anchored systems)
                # =====================================================================
                if rt_gene is not None:
                    system_protein_seq = rt_gene.get("sequence", "")
                    
                    # Remove stop codons for fair comparison
                    if system_protein_seq.endswith('*'):
                        system_protein_seq = system_protein_seq[:-1]
                    
                    best_gt_protein_similarity = 0.0
                    matched_gt_protein_id = None
                    
                    for gt_protein in gt_proteins:
                        gt_protein_seq = gt_protein.get("sequence", "")
                        if gt_protein_seq.endswith('*'):
                            gt_protein_seq = gt_protein_seq[:-1]
                        
                        protein_similarity = calculate_similarity(gt_protein_seq, system_protein_seq)
                        
                        if protein_similarity > best_gt_protein_similarity:
                            best_gt_protein_similarity = protein_similarity
                            matched_gt_protein_id = gt_protein.get("id", "unknown")
                    
                    protein_match = best_gt_protein_similarity >= self.args.protein_similarity_threshold
                    
                    match_info["protein_match"] = protein_match
                    match_info["protein_similarity"] = round(best_gt_protein_similarity, 4)
                    match_info["matched_gt_protein_id"] = matched_gt_protein_id if protein_match else None
                    match_info["protein_length_gt"] = len(gt_proteins[0].get("sequence", "").rstrip('*')) if gt_proteins else 0
                    match_info["protein_length_detected"] = len(system_protein_seq)
                    
                    if protein_match:
                        match_info["matched_elements"].append("protein")
                        validation_results["summary"]["any_protein_match"] = True
                else:
                    # ncRNA-anchored system: no protein to compare
                    match_info["protein_match"] = None  # Special value to indicate no protein
                    self.logger.debug(f"System {rt_system_id} is ncRNA-anchored (no RT gene), skipping protein validation")
                
                # =====================================================================
                # UPDATED: ncRNA MATCHING - ONLY CHECK IF THIS IS THE BEST PROTEIN SYSTEM
                # OR IF IT'S AN ncRNA-ANCHORED SYSTEM
                # =====================================================================
                should_check_ncrna = False
                
                # Check ncRNA if:
                # 1. This is the best protein system AND we have ground truth ncRNAs
                # 2. OR this is an ncRNA-anchored system AND we have ground truth ncRNAs
                if gt_ncrnas:
                    if rt_gene is None:  # ncRNA-anchored system
                        should_check_ncrna = True
                    elif best_protein_system and rt_system_id == best_protein_system["rt_system_id"]:
                        should_check_ncrna = True


                if should_check_ncrna:
                    # Get ncRNAs from THIS system only
                    ncrnas = system.get("ncrnas", [])
                    
                    self.logger.info(f"DEBUG: Checking ncRNA for system {rt_system_id}")
                    self.logger.info(f"DEBUG: Found {len(ncrnas)} ncRNAs in system.ncrnas")
                    
                    # Fallback to legacy PADLOC location
                    if not ncrnas:
                        ncrnas = self.safe_nested_get(system, "rt_gene", "tool_metadata", "PADLOC", "ncrnas", default=[])
                        self.logger.info(f"DEBUG: Fallback to PADLOC found {len(ncrnas)} ncRNAs")
                    
                    if not ncrnas:
                        self.logger.warning(f"DEBUG: No ncRNAs found in system {rt_system_id} despite should_check_ncrna=True")
                    

                
                    
                    best_system_ncrna_similarity = 0.0
                    best_ncrna_info = None
                    matched_gt_ncrna_id = None
                    
                    # Only check ncRNAs if there are any in THIS system
                    if ncrnas:
                        self.logger.info(f"DEBUG: System has {len(ncrnas)} ncRNAs, comparing against {len(gt_ncrnas)} GT ncRNAs")
                        
                        for i, ncrna in enumerate(ncrnas):
                            ncrna_seq = ncrna.get("sequence", "").upper()
                            self.logger.info(f"DEBUG: System ncRNA {i}: length={len(ncrna_seq)}, type={ncrna.get('type')}")
                            
                            for j, gt_ncrna in enumerate(gt_ncrnas):
                                gt_ncrna_seq = gt_ncrna.get("sequence", "")
                                gt_type = gt_ncrna.get("type", "RNA")
                                
                                self.logger.info(f"DEBUG: GT ncRNA {j}: length={len(gt_ncrna_seq)}, type={gt_type}")
                                self.logger.info(f"DEBUG: GT sequence (first 50bp): {gt_ncrna_seq[:50]}")
                                


                             
                                # Convert DNA to RNA if needed
                                ncrna_seq_normalized = ncrna_seq
                                if gt_type == "RNA" and 'T' in ncrna_seq and 'U' not in ncrna_seq:
                                    ncrna_seq_normalized = ncrna_seq.replace('T', 'U')
                                elif gt_type == "DNA" and 'U' in ncrna_seq:
                                    gt_ncrna_seq = gt_ncrna_seq.replace('T', 'U')


                                self.logger.info(f"DEBUG: System sequence (first 50bp): {ncrna_seq_normalized[:50]}")
                                
                                # IMPROVED: Use alignment-based comparison
                                try:
                                    # Try using Biopython alignment if available
                                    from Bio import Align
                                    aligner = Align.PairwiseAligner()
                                    aligner.mode = 'local'
                                    aligner.match_score = 2
                                    aligner.mismatch_score = -1
                                    aligner.open_gap_score = -2
                                    aligner.extend_gap_score = -0.5
                                    
                                    alignments = aligner.align(gt_ncrna_seq.upper(), ncrna_seq_normalized)
                                    if alignments:
                                        score = alignments[0].score
                                        max_score = min(len(gt_ncrna_seq), len(ncrna_seq_normalized)) * aligner.match_score
                                        ncrna_similarity = max(0.0, min(1.0, score / max_score)) if max_score > 0 else 0.0
                                    else:
                                        ncrna_similarity = 0.0
                                    
                                    self.logger.info(f"DEBUG: Alignment similarity = {ncrna_similarity}")
                                    
                                except ImportError:
                                    # Fallback to simple comparison if Biopython not available
                                    ncrna_similarity = calculate_similarity(gt_ncrna_seq, ncrna_seq_normalized)
                                    self.logger.info(f"DEBUG: Simple similarity = {ncrna_similarity}")
                                
                                if ncrna_similarity > best_system_ncrna_similarity:
                                    best_system_ncrna_similarity = ncrna_similarity
                                    matched_gt_ncrna_id = gt_ncrna.get("id", "unknown")
                                    best_ncrna_info = {
                                        "ncrna_id": ncrna.get("ncrna_id", ""),
                                        "type": ncrna.get("type", ""),
                                        "start": ncrna.get("start"),
                                        "end": ncrna.get("end"),
                                        "strand": ncrna.get("strand", "+"),
                                        "length": ncrna.get("length", len(ncrna_seq)),
                                        "intergenic_region_id": ncrna.get("intergenic_region_id"),
                                        "has_structure_annotations": "structure_annotations" in ncrna
                                    }
                        
                        self.logger.info(f"DEBUG: Best ncRNA similarity = {best_system_ncrna_similarity}, threshold = {self.args.ncrna_similarity_threshold}")


                    
                    best_system_ncrna_similarity = 0.0
                    best_ncrna_info = None
                    matched_gt_ncrna_id = None
                    
                    # Only check ncRNAs if there are any in THIS system
                    if ncrnas:
                        for ncrna in ncrnas:
                            ncrna_seq = ncrna.get("sequence", "").upper()
                            
                            for gt_ncrna in gt_ncrnas:
                                gt_ncrna_seq = gt_ncrna.get("sequence", "")
                                gt_type = gt_ncrna.get("type", "RNA")
                                
                                # Convert DNA to RNA if needed
                                ncrna_seq_normalized = ncrna_seq
                                if gt_type == "RNA" and 'T' in ncrna_seq and 'U' not in ncrna_seq:
                                    ncrna_seq_normalized = ncrna_seq.replace('T', 'U')
                                elif gt_type == "DNA" and 'U' in ncrna_seq:
                                    gt_ncrna_seq = gt_ncrna_seq.replace('T', 'U')
                                
                                # IMPROVED: Use alignment-based comparison
                                try:
                                    # Try using Biopython alignment if available
                                    from Bio import Align
                                    aligner = Align.PairwiseAligner()
                                    aligner.mode = 'local'  # Local alignment for boundary differences
                                    aligner.match_score = 2
                                    aligner.mismatch_score = -1
                                    aligner.open_gap_score = -2
                                    aligner.extend_gap_score = -0.5
                                    
                                    alignments = aligner.align(gt_ncrna_seq.upper(), ncrna_seq_normalized)
                                    if alignments:
                                        score = alignments[0].score
                                        max_score = min(len(gt_ncrna_seq), len(ncrna_seq_normalized)) * aligner.match_score
                                        ncrna_similarity = max(0.0, min(1.0, score / max_score)) if max_score > 0 else 0.0
                                    else:
                                        ncrna_similarity = 0.0
                                except ImportError:
                                    # Fallback to simple comparison if Biopython not available
                                    ncrna_similarity = calculate_similarity(gt_ncrna_seq, ncrna_seq_normalized)
                                
                                if ncrna_similarity > best_system_ncrna_similarity:
                                    best_system_ncrna_similarity = ncrna_similarity
                                    matched_gt_ncrna_id = gt_ncrna.get("id", "unknown")
                                    best_ncrna_info = {
                                        "ncrna_id": ncrna.get("ncrna_id", ""),
                                        "type": ncrna.get("type", ""),
                                        "start": ncrna.get("start"),
                                        "end": ncrna.get("end"),
                                        "strand": ncrna.get("strand", "+"),
                                        "length": ncrna.get("length", len(ncrna_seq)),
                                        "intergenic_region_id": ncrna.get("intergenic_region_id"),
                                        "has_structure_annotations": "structure_annotations" in ncrna
                                    }
                    
                    # Update match info with ncRNA results
                    if best_system_ncrna_similarity > 0:
                        match_info["ncrna_similarity"] = round(best_system_ncrna_similarity, 4)
                        match_info["ncrna_match"] = best_system_ncrna_similarity >= self.args.ncrna_similarity_threshold
                        match_info["matched_gt_ncrna_id"] = matched_gt_ncrna_id if match_info["ncrna_match"] else None
                        match_info["ncrna_info"] = best_ncrna_info
                        match_info["ncrna_length_gt"] = len(gt_ncrnas[0].get("sequence", "")) if gt_ncrnas else 0
                        match_info["ncrna_length_detected"] = best_ncrna_info.get("length", 0) if best_ncrna_info else 0
                        
                        if match_info["ncrna_match"]:
                            match_info["matched_elements"].append("ncRNA")
                            validation_results["summary"]["any_ncrna_match"] = True
                    
                    # Special case: For ncRNA-anchored systems, also check if we should add to matches
                    if rt_gene is None and match_info["ncrna_match"]:
                        # Only add ncRNA-anchored systems if they have an ncRNA match
                        validation_results["matches_found"].append(match_info)
                        continue  # Skip the complete match logic for ncRNA-anchored systems
                
                # =====================================================================
                # Track complete matches (protein + ncRNA in SAME system)
                # Only for RT-anchored systems
                # =====================================================================
                if rt_gene is not None and match_info["protein_match"] and match_info["ncrna_match"]:
                    validation_results["summary"]["complete_match_systems"].append(rt_system_id)
                    
                    # Track best complete match
                    combined_score = (match_info["protein_similarity"] + match_info["ncrna_similarity"]) / 2
                    if combined_score > best_complete_similarity:
                        best_complete_similarity = combined_score
                        best_complete_system = {
                            "rt_system_id": rt_system_id,
                            "protein_similarity": match_info["protein_similarity"],
                            "ncrna_similarity": match_info["ncrna_similarity"],
                            "combined_score": round(combined_score, 4)
                        }
                
                # Only add to matches if protein OR ncRNA matched
                # For RT-anchored systems: require at least one match
                # For ncRNA-anchored systems: already handled above
                if rt_gene is not None and (match_info["protein_match"] or match_info["ncrna_match"]):
                    validation_results["matches_found"].append(match_info)
            
            # Set best complete match
            validation_results["best_complete_match"] = best_complete_system
            
            # Add validation to ground truth
            ground_truth["validation"] = {
                "validated": len(validation_results["matches_found"]) > 0,
                "validation_timestamp": datetime.now().isoformat(),
                "integrated_results_file": str(Path(integrated_results_path).name),
                "thresholds": {
                    "protein_similarity": self.args.protein_similarity_threshold,
                    "ncrna_similarity": self.args.ncrna_similarity_threshold
                },
                "summary": validation_results["summary"],
                "matches": validation_results["matches_found"],
                "best_matches": {
                    "protein": validation_results["best_protein_match"],
                    "complete": validation_results.get("best_complete_match")
                }
            }
            
            # Save updated ground truth
            with open(gt_metadata_path, 'w') as f:
                json.dump(ground_truth, f, indent=2)
            
            # Log summary
            self.logger.info("="*80)
            self.logger.info("VALIDATION SUMMARY")
            self.logger.info("="*80)
            self.logger.info(f"Ground Truth ID: {validation_results['ground_truth_id']}")
            self.logger.info(f"Total Systems Checked: {validation_results['total_systems_checked']}")
            self.logger.info(f"Matches Found: {len(validation_results['matches_found'])}")
            
            if validation_results['best_protein_match'] and validation_results['best_protein_match']['rt_system_id']:
                self.logger.info(f"✓ Best Protein Match: {validation_results['best_protein_match']['rt_system_id']} "
                            f"(similarity: {validation_results['best_protein_match']['similarity']})")
            
            if validation_results.get('best_complete_match'):
                bcm = validation_results['best_complete_match']
                self.logger.info(f"✓ Best Complete Match: {bcm['rt_system_id']} "
                            f"(protein: {bcm['protein_similarity']}, "
                            f"ncRNA: {bcm['ncrna_similarity']}, "
                            f"combined: {bcm['combined_score']})")
            
            if validation_results['summary']['complete_match_systems']:
                self.logger.info(f"✓ Complete Matches (protein + ncRNA in same system): "
                            f"{', '.join(validation_results['summary']['complete_match_systems'])}")
            
            # Log ncRNA-anchored matches separately
            ncrna_anchored_matches = [m for m in validation_results["matches_found"] if m.get("anchor_type") == "ncRNA"]
            if ncrna_anchored_matches:
                self.logger.info(f"✓ ncRNA-anchored system matches: {len(ncrna_anchored_matches)}")
            
            self.logger.info("="*80)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error during validation: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False


    def process_genome(self, genome_dir: Path, output_base_dir: Path) -> bool:

        genome_id = genome_dir.name
        self.logger.info(f"\n{'='*80}")
        self.logger.info(f"Processing genome: {genome_id}")
        self.logger.info(f"{'='*80}\n")

        if not self.validate_genome_directory(genome_dir):
            return False

        is_ground_truth = self.check_ground_truth(genome_dir)
        genome_output_dir = output_base_dir / genome_id
        genome_output_dir.mkdir(parents=True, exist_ok=True)

        # genome_fasta = genome_dir / "genome.fasta"

        genome_fasta = next(genome_dir.glob("genome.fasta*"))


        results = {
            "genome_id": genome_id,
            "input_dir": str(genome_dir),
            "output_dir": str(genome_output_dir),
            "is_ground_truth": is_ground_truth,
            "timestamp": datetime.now().isoformat(),
            "steps": {}
        }

        steps = [
            ("step0_ground_truth",        lambda: self.step0_extract_ground_truth(genome_dir, genome_output_dir)),
            ("step0b_annotate_gt_ncrna",  lambda: self.step0b_annotate_gt_ncrna(genome_dir, genome_output_dir)),
            ("step1",                     lambda: self.step1_run_analysis_tools(genome_fasta, genome_output_dir)),
            ("step2a_myrt",               lambda: self.step2_extract_myrt(genome_output_dir)),
            ("step2b_padloc",             lambda: self.step2_extract_padloc(genome_output_dir)),
            ("step2c_defensefinder",      lambda: self.step2_parse_defensefinder(genome_output_dir)),
            # ("step3_integrate",           lambda: self.step3_integrate_results(genome_output_dir, genome_fasta)),
            ("step3_integrate",           lambda: self.step3_integrate_results(genome_output_dir)),
            ("step4_annotate",            lambda: self.step4_annotate_ncrna(genome_output_dir)),
            ("step5_validate",            lambda: self.step5_validate_results(genome_dir, genome_output_dir)),
        ]

        all_successful = True
        step2_successes = 0

        for step_name, step_func in steps:

            self.logger.info(f"\n--- Running {step_name} ---")
            success = step_func()
            results["steps"][step_name] = "success" if success else "failed"

            # Track Step2 group results
            if step_name.startswith("step2"):
                if success:
                    step2_successes += 1

            # Normal fatal behavior for non-step2
            if (not step_name.startswith("step2")) and (not success):
                all_successful = False
                if not self.args.continue_on_error:
                    self.logger.error(f"Pipeline stopped at {step_name} for {genome_id}")
                    break

            # After finishing step2c → evaluate group
            if step_name == "step2c_defensefinder":
                if step2_successes == 0:
                    self.logger.error(
                        f"\n❌ All Step2 tools failed for {genome_id}. "
                        f"No MyRT/PADLOC/DefenseFinder results. Cannot continue.\n"
                    )
                    results["steps"]["step2_group"] = "all_failed"
                    return False
                else:
                    self.logger.info(
                        f"✓ Step2 group summary: {step2_successes}/3 succeeded. Continuing.\n"
                    )
                    results["steps"]["step2_group"] = f"{step2_successes}_of_3_success"

        # Save pipeline results
        results_file = genome_output_dir / "pipeline_results.json"
        with open(results_file, "w") as f:
            json.dump(results, f, indent=2)

        if all_successful:
            self.logger.info(f"\n✓ Successfully completed all steps for {genome_id}\n")
        else:
            self.logger.warning(f"\n⚠ Some steps failed for {genome_id}\n")

        return all_successful

    
    def run(self):
        """Main execution method."""
        self.logger.info("="*80)
        self.logger.info("Retron Pipeline Orchestrator")
        self.logger.info("="*80)
        
        # Create output base directory
        output_base_dir = Path(self.args.output)
        output_base_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Output base directory: {output_base_dir}")
        
        # Determine input mode
        if self.args.batch:
            # Batch mode: process all subdirectories
            batch_dir = Path(self.args.batch)
            if not batch_dir.exists():
                self.logger.error(f"Batch directory does not exist: {batch_dir}")
                return False
            
            genome_dirs = [d for d in batch_dir.iterdir() if d.is_dir()]
            self.logger.info(f"Found {len(genome_dirs)} genome directories to process")
            
            successful = 0
            failed = 0
            
            for genome_dir in genome_dirs:
                try:
                    if self.process_genome(genome_dir, output_base_dir):
                        successful += 1
                    else:
                        failed += 1
                except Exception as e:
                    self.logger.error(f"Exception while processing {genome_dir}: {str(e)}")
                    failed += 1
            
            self.logger.info("\n" + "="*80)
            self.logger.info("Batch Processing Summary")
            self.logger.info("="*80)
            self.logger.info(f"Total genomes: {len(genome_dirs)}")
            self.logger.info(f"Successful: {successful}")
            self.logger.info(f"Failed: {failed}")
            
            return failed == 0
            
        else:
            # Single genome mode
            input_dir = Path(self.args.input)
            if not input_dir.exists():
                self.logger.error(f"Input directory does not exist: {input_dir}")
                return False
            
            return self.process_genome(input_dir, output_base_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Orchestrate the retron analysis pipeline for bacterial genomes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single genome (auto-detect threads)
  python retron_pipeline_orchestrator.py \\
    -i /path/to/genome/terminal_6 \\
    -o /path/to/output/gt_results \\
    --run-infernal \\
    --annotate-ncrna
  
  # Manually specify threads
  python retron_pipeline_orchestrator.py \\
    -i /path/to/genome/terminal_6 \\
    -o /path/to/output/gt_results \\
    --threads 16
  
  # Batch process multiple genomes
  python retron_pipeline_orchestrator.py \\
    -b /path/to/genomes/positive_data2 \\
    -o /path/to/output/gt_results \\
    --run-infernal \\
    --annotate-ncrna
        """
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-i', '--input',
        type=str,
        help='Input genome directory (single genome mode)'
    )
    input_group.add_argument(
        '-b', '--batch',
        type=str,
        help='Batch input directory containing multiple genome subdirectories'
    )
    
    # Output
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output base directory for results'
    )
    
    # Analysis options
    parser.add_argument(
        '--run-infernal',
        action='store_true',
        help='Run Infernal for ncRNA structure prediction'
    )
    parser.add_argument(
        '--annotate-ncrna',
        action='store_true',
        help='Annotate ground truth ncRNA sequences with bpRNA (requires one_sequence_bprna_annotation.py)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=None,
        metavar='N',
        help=f'Number of threads to use (default: auto = {get_optimal_threads()} threads, 80%% of {multiprocessing.cpu_count()} available CPUs)'
    )
    parser.add_argument(
        '--evalue',
        type=float,
        help='E-value threshold for similarity searches'
    )
    
    # Pipeline control
    parser.add_argument(
        '--continue-on-error',
        action='store_true',
        help='Continue to next steps even if a step fails'
    )
    
    # Logging
    parser.add_argument(
        '--log',
        type=str,
        help='Log file path (default: print to stdout)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '--protein-similarity-threshold', 
        type=float, 
        default=0.95,
        help='Minimum similarity threshold for protein match (default: 0.95)'
    )
    parser.add_argument(
        '--ncrna-similarity-threshold', 
        type=float, 
        default=0.10,
        help='Minimum similarity threshold for ncRNA match (default: 0.10)'
    )
    
    args = parser.parse_args()
    
    # Auto-detect threads if not specified
    # if args.threads is None:
    #     args.threads = get_optimal_threads()

    print("THREADS fixed to 4 for now")
    args.threads = 4
    
    # Create and run orchestrator
    orchestrator = RetronPipelineOrchestrator(args)
    success = orchestrator.run()
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()