#!/usr/bin/env python3
"""
Comprehensive Retron Operon Annotation Pipeline
Integrates MyRT, PADLOC, and DefenseFinder outputs

This specific version of the script uses Infernal together with PADLOC.

check this conversation, some useful commands to use infernal https://claude.ai/chat/3d6a7041-81e5-484f-ac3d-c68d4bedba54 

"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import pandas as pd
import subprocess
from pathlib import Path

import glob



# =============================================================================
# LOGGER SETUP
# =============================================================================

def setup_logger(log_file=None):
    """Setup logging configuration"""
    logger = logging.getLogger('retron_annotator')
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # File handler if specified
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
    return logger

# =============================================================================
# GENOME PROCESSING
# =============================================================================

def load_genome(genome_fasta, logger):
    """Load genome sequences from FASTA file"""
    logger.info(f"Loading genome from {genome_fasta}")
    
    contigs = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        contigs[record.id] = {
            'sequence': str(record.seq),
            'length': len(record.seq),
            'description': record.description
        }
    
    logger.info(f"Loaded {len(contigs)} contig(s)")
    return contigs


def run_prodigal_once(input_fasta, output_dir, prodigal_bin='prodigal'):
    """
    Run Prodigal gene prediction once for all tools to share
    Outputs go to prodigal_results subdirectory
    Handles both plain and gzipped FASTA files
    Sanitizes contig names if needed (creates a copy, doesn't modify original)
    
    Returns:
        Tuple of (faa_file, gff_file, mapping_file)
    """
    from pathlib import Path
    import subprocess
    import gzip
    import shutil
    
    output_dir = Path(output_dir)
    input_fasta = Path(input_fasta)
    
    # Create Prodigal subdirectory
    prodigal_dir = output_dir / "prodigal_results"
    prodigal_dir.mkdir(parents=True, exist_ok=True)
    
    # Output files in prodigal_results subdirectory
    gff_file = prodigal_dir / "genes.gff"
    faa_file_raw = prodigal_dir / "proteins_raw.faa"
    faa_file = prodigal_dir / "proteins.faa"
    
    # Detect if file is gzipped and decompress if needed
    actual_input = input_fasta
    if input_fasta.suffix == '.gz' or is_gzipped(input_fasta):
        print("[INFO] Input is gzipped - decompressing in place...")
        # Decompress right next to the original .gz file
        decompressed = input_fasta.with_suffix('')  # removes .gz extension
        
        with gzip.open(input_fasta, 'rb') as f_in:
            with open(decompressed, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        actual_input = decompressed
        print(f"[INFO] Decompressed to: {actual_input}")
    
    # Check if contig names need sanitization
    with open(actual_input) as f:
        first_line = f.readline().strip()
    
    if ' ' in first_line or '\t' in first_line:
        print("[INFO] Contig names contain spaces - creating sanitized copy...")
        needs_sanitization = True
        
        # Create sanitized FASTA in the output directory (not modifying original!)
        sanitized_fasta = prodigal_dir / "genome_sanitized.fasta"
        print(f"[INFO] Created sanitized copy: {sanitized_fasta}")
        
        with open(actual_input) as f_in, open(sanitized_fasta, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    # Take only the first part before any space/tab
                    contig_id = line.split()[0]
                    f_out.write(contig_id + '\n')
                else:
                    f_out.write(line)
        
        actual_input = sanitized_fasta
        print(f"[INFO] Created sanitized FASTA: {sanitized_fasta}")
        print(f"[INFO] Original files remain unchanged")
    else:
        print("[INFO] Contig names are clean - using original")
        needs_sanitization = False
    
    # Run Prodigal
    cmd = [
        prodigal_bin,
        "-i", str(actual_input),
        "-o", str(gff_file),
        "-a", str(faa_file_raw),
        "-p", "meta",
        "-f", "gff"
    ]
    
    print(f"[INFO] Running Prodigal: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # Copy raw FAA to final FAA
    shutil.copy(faa_file_raw, faa_file)
    
    print(f"[INFO] Prodigal complete:")
    print(f"       GFF: {gff_file}")
    print(f"       FAA: {faa_file}")
    print(f"       Working FASTA: {actual_input}")
    
    # return str(faa_file), str(gff_file), None
    return str(faa_file), str(gff_file), str(actual_input)



def is_gzipped(filepath):
    """Check if a file is gzipped by reading its magic number"""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'



def run_infernal(
    genome_fasta,
    output_dir,
    cpu=4,
    extract_sequences=True
):
    """
    Run Infernal to detect retron ncRNAs using PADLOC's CM database
    
    Returns:
        Tuple of (ncrna_table_path, ncrna_sequences_path)
    """
    from pathlib import Path
    import subprocess
    
    output_dir = Path(output_dir)
    
    # Create Infernal subdirectory
    infernal_dir = output_dir / "infernal_results"
    infernal_dir.mkdir(parents=True, exist_ok=True)
    
    # Use padloc2 environment's CM database
    infernal_cm_db = "/ibex/user/rioszemm/conda-environments/padloc2/data/cm/padlocdb.cm"
    
    genome_name = Path(genome_fasta).stem
    
    # All outputs go in infernal_results subdirectory
    ncrna_table_raw = infernal_dir / f"{genome_name}_ncrna_raw.tblout"
    ncrna_table = infernal_dir / f"{genome_name}_ncrna.tblout"
    ncrna_alignment = infernal_dir / f"{genome_name}_ncrna.sto"
    ncrna_sequences_all = infernal_dir / f"{genome_name}_ncrna_sequences_all.fa"
    ncrna_sequences_high_conf = infernal_dir / f"{genome_name}_ncrna_sequences_highconf.fa"
    ncrna_summary = infernal_dir / f"{genome_name}_ncrna_summary.txt"
    
    # Delete existing output if present
    if ncrna_table.exists():
        ncrna_table.unlink()
    
    # Step 1: Run cmsearch
    infernal_cmd = [
        "cmsearch",
        "--cpu", str(cpu),
        "-Z", "10",
        "--FZ", "500",
        "--acc",
        "--noali",
        "--tblout", str(ncrna_table_raw),
        str(infernal_cm_db),
        str(genome_fasta)
    ]
    
    print(f"[INFO] Running Infernal for retron ncRNA detection")
    print(f"[CMD] {' '.join(infernal_cmd)}")
    
    try:
        result = subprocess.run(
            infernal_cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print(f"[INFO] Infernal completed successfully")
        
        # Step 2: Parse and format output
        with open(ncrna_table_raw) as f:
            lines = f.readlines()
        
        header_lines = []
        data_lines = []
        high_conf_lines = []
        
        # Parse hits into structured data
        hit_metadata = {}  # key: (seqid, start, end), value: metadata dict
        
        for line in lines:
            if line.startswith('#'):
                header_lines.append(line)
            elif line.strip():
                data_lines.append(line)
                fields = line.split()
                
                if len(fields) >= 17:
                    seqid = fields[0]
                    model = fields[2]
                    seq_from = int(fields[7])
                    seq_to = int(fields[8])
                    strand = fields[9]
                    score = float(fields[14])
                    evalue = fields[15]
                    confidence = fields[16]  # '!' or '?'
                    
                    # Normalize coordinates (always start < end)
                    start = min(seq_from, seq_to)
                    end = max(seq_from, seq_to)
                    
                    key = (seqid, start, end)
                    hit_metadata[key] = {
                        'model': model,
                        'strand': strand,
                        'score': score,
                        'evalue': evalue,
                        'confidence': confidence
                    }
                    
                    if confidence == '!':
                        high_conf_lines.append(line)
        
        if data_lines:
            print(f"[INFO] Found {len(data_lines)} total ncRNA hits")
            print(f"[INFO]   - {len(high_conf_lines)} high-confidence hits (!)")
            print(f"[INFO]   - {len(data_lines) - len(high_conf_lines)} low-confidence hits (?)")
            
            # Step 3: Convert to tab-delimited and fix column name
            with open(ncrna_table, 'w') as outf:
                # Write headers - convert to tab-delimited AND fix mdl → mdl.type
                for line in header_lines:
                    if line.startswith('#target name'):
                        # First fix the column name
                        line = line.replace(' mdl mdl from', ' mdl.type mdl from', 1)
                        # Then convert to tab-delimited
                        parts = line.strip().split()
                        line = '\t'.join(parts) + '\n'
                    outf.write(line)
                
                # Write data as tab-delimited
                for line in data_lines:
                    fields = line.split()
                    
                    if len(fields) >= 17:
                        # First 17 fields as tabs
                        tab_line = '\t'.join(fields[:17])
                        
                        # Description (remaining fields) as spaces
                        if len(fields) > 17:
                            description = ' '.join(fields[17:])
                            tab_line += '\t' + description
                        
                        outf.write(tab_line + '\n')
            
            # Step 4: Create summary file
            with open(ncrna_summary, 'w') as f:
                f.write(f"Infernal ncRNA Detection Summary\n")
                f.write(f"=" * 60 + "\n\n")
                f.write(f"Genome: {genome_name}\n")
                f.write(f"Database: {infernal_cm_db}\n\n")
                f.write(f"Total hits: {len(data_lines)}\n")
                f.write(f"High-confidence hits (!): {len(high_conf_lines)}\n")
                f.write(f"Low-confidence hits (?): {len(data_lines) - len(high_conf_lines)}\n\n")
                
                if high_conf_lines:
                    f.write("High-confidence hits:\n")
                    f.write("-" * 60 + "\n")
                    for line in high_conf_lines:
                        fields = line.split()
                        if len(fields) >= 17:
                            f.write(f"  {fields[2]:20s} E-value: {fields[15]:10s} Coords: {fields[7]}-{fields[8]}\n")
                
                f.write("\nLow-confidence hits:\n")
                f.write("-" * 60 + "\n")
                for line in data_lines:
                    fields = line.split()
                    if len(fields) >= 17 and fields[16] != '!':
                        f.write(f"  {fields[2]:20s} E-value: {fields[15]:10s} Coords: {fields[7]}-{fields[8]}\n")
            
            print(f"[INFO] Summary written: {ncrna_summary}")
            
            # Step 5: Extract sequences if requested
            if extract_sequences:
                print(f"[INFO] Extracting ncRNA sequences with enhanced headers...")
                
                # Load genome sequences
                from Bio import SeqIO
                genome_seqs = {rec.id: rec for rec in SeqIO.parse(genome_fasta, "fasta")}
                
                # Extract sequences with rich metadata
                with open(ncrna_sequences_all, 'w') as f_all, \
                     open(ncrna_sequences_high_conf, 'w') as f_high:
                    
                    for (seqid, start, end), metadata in sorted(hit_metadata.items()):
                        if seqid not in genome_seqs:
                            print(f"[WARNING] Sequence {seqid} not found in genome")
                            continue
                        
                        # Extract sequence region
                        seq_record = genome_seqs[seqid]
                        # Convert to 0-based indexing
                        subseq = seq_record.seq[start-1:end]
                        
                        # Reverse complement if needed
                        if metadata['strand'] == '-':
                            subseq = subseq.reverse_complement()
                        
                        # Create enhanced FASTA header with metadata
                        confidence_str = "high_conf" if metadata['confidence'] == '!' else "low_conf"
                        header = (
                            f">{seqid}:{start}-{end}({metadata['strand']}) "
                            f"model={metadata['model']} "
                            f"evalue={metadata['evalue']} "
                            f"score={metadata['score']:.1f} "
                            f"confidence={confidence_str}"
                        )
                        
                        # Write to all sequences file
                        f_all.write(f"{header}\n")
                        # Write sequence with line breaks every 60 chars
                        seq_str = str(subseq)
                        for i in range(0, len(seq_str), 60):
                            f_all.write(f"{seq_str[i:i+60]}\n")
                        
                        # Write to high-confidence file if applicable
                        if metadata['confidence'] == '!':
                            f_high.write(f"{header}\n")
                            for i in range(0, len(seq_str), 60):
                                f_high.write(f"{seq_str[i:i+60]}\n")
                
                print(f"[INFO] All ncRNA sequences saved: {ncrna_sequences_all}")
                if high_conf_lines:
                    print(f"[INFO] High-confidence ncRNA sequences saved: {ncrna_sequences_high_conf}")
            
            # Cleanup raw file
            ncrna_table_raw.unlink()
            
            return str(ncrna_table), str(ncrna_sequences_all) if extract_sequences else None
            
        else:
            print(f"[INFO] No retron ncRNA hits found")
            # Create empty file with headers
            with open(ncrna_table, 'w') as outf:
                for line in header_lines:
                    if line.startswith('#target name'):
                        line = line.replace(' mdl mdl from', ' mdl.type mdl from', 1)
                    outf.write(line)
            
            ncrna_table_raw.unlink()
            return str(ncrna_table), None
            
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Infernal failed: {e}")
        if e.stderr:
            print(f"[ERROR] stderr: {e.stderr}")
        return None, None
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def run_padloc(
    output_dir,
    gff_file=None,
    faa_file=None,
    fna_file=None,
    mapping_file=None,
    ncrna_file=None,
    restore_names=False,
    cpu=4,
    padloc_bin=None  # Ignored - we use conda environment instead
):
    """
    Run PADLOC v2.0.0 for defense system detection via conda environment
    
    Args:
        output_dir: Output directory
        gff_file: GFF file (optional, used with faa_file)
        faa_file: Protein FASTA file (optional, used with gff_file)
        fna_file: Nucleotide FASTA file (optional, alternative to faa+gff)
        mapping_file: Contig mapping file (ignored for now)
        ncrna_file: Infernal ncRNA results file (optional)
        restore_names: Whether to restore original contig names (ignored for now)
        cpu: Number of CPUs
        padloc_bin: Path to PADLOC binary (ignored - using conda environment)
    
    Returns:
        Path to PADLOC CSV results file or None if failed
    """
    from pathlib import Path
    import subprocess
    import glob
    
    output_dir = Path(output_dir)
    padloc_output = output_dir / "padloc_results"
    padloc_output.mkdir(parents=True, exist_ok=True)
    
    # Build PADLOC command using conda environment
    padloc_cmd = [
        "conda", "run", "-n", "padloc2",
        "padloc",
        "--outdir", str(padloc_output),
        "--cpu", str(cpu)
    ]


    # # Build PADLOC command using direct binary path
    # padloc_bin_path = "/ibex/user/rioszemm/the-retron-project/src/padloc/bin/padloc"
    # padloc_cmd = [
    #     padloc_bin_path,
    #     "--outdir", str(padloc_output),
    #     "--cpu", str(cpu)
    # ]
        
    # padloc_cmd = [
    #     "padloc",  # Will use the one in current conda environment
    #     "--outdir", str(padloc_output),
    #     "--cpu", str(cpu)
    # ]


    # Add input files - prefer faa+gff if available, otherwise fna
    if faa_file and gff_file:
        padloc_cmd.extend([
            "--faa", str(faa_file),
            "--gff", str(gff_file),
            "--fix-prodigal"  # Important for Prodigal-generated files
        ])
    elif fna_file:
        padloc_cmd.extend(["--fna", str(fna_file)])
    else:
        raise ValueError("Must provide either (faa+gff) or fna file")
    
    # Add ncRNA file if provided
    if ncrna_file and Path(ncrna_file).exists():
        padloc_cmd.extend(["--ncrna", str(ncrna_file)])
        print(f"[INFO] Including retron ncRNA data from Infernal")
    
    print(f"[INFO] Running PADLOC v2.0.0")
    print(f"[CMD] {' '.join(padloc_cmd)}")
    
    try:
        result = subprocess.run(
            padloc_cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        print(f"[INFO] PADLOC completed successfully")
        print(result.stdout)
        
        # Find the output CSV file
        csv_files = list(padloc_output.glob("*_padloc.csv"))
        
        if csv_files:
            return str(csv_files[0])
        else:
            print("[WARNING] No PADLOC CSV output found")
            return None
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] PADLOC failed: {e}")
        if e.stderr:
            print(f"[ERROR] stderr: {e.stderr}")
        if e.stdout:
            print(f"[INFO] stdout: {e.stdout}")
        return None
    except Exception as e:
        print(f"[ERROR] Unexpected error running PADLOC: {e}")
        import traceback
        traceback.print_exc()
        return None

        
def run_myrt(genome_fasta, gff_file, output_dir, config, logger):
    """
    Run MyRT on nucleotide genome using pre-generated Prodigal GFF
    
    Args:
        genome_fasta: Path to nucleotide genome FASTA file (.fna/.fasta)
        gff_file: Path to Prodigal-generated GFF file
        output_dir: Output directory
        config: Configuration dict
        logger: Logger instance
    
    Returns:
        Path to MyRT results directory
    """
    logger.info("=" * 80)
    logger.info("Running MyRT for RT Detection & Classification (using Prodigal)")
    logger.info("=" * 80)
    
    # myrt_output = os.path.join(output_dir, "myrt_output")
    # os.makedirs(myrt_output, exist_ok=True)
    
    # # MyRT expects: genome_id folder with genome_id.fna and genome_id.gff
    # genome_id = Path(genome_fasta).stem
    # myrt_input_dir = os.path.join(myrt_output, genome_id)

    myrt_output = os.path.join(output_dir, "myrt_output")
    os.makedirs(myrt_output, exist_ok=True)

    # Use myrt_output directly (no subdirectory)
    genome_id = Path(genome_fasta).stem
    myrt_input_dir = myrt_output  # <-- Changed: use myrt_output directly
    os.makedirs(myrt_input_dir, exist_ok=True)
    
    # Copy nucleotide genome file
    myrt_fna = os.path.join(myrt_input_dir, f"{genome_id}.fna")
    if not os.path.exists(myrt_fna):
        import shutil
        shutil.copy(os.path.abspath(genome_fasta), myrt_fna)
        logger.info(f"Copied nucleotide genome to: {myrt_fna}")
    
    # Copy GFF file - MyRT's rerun.sh will detect it and use it
    myrt_gff = os.path.join(myrt_input_dir, f"{genome_id}.gff")
    if not os.path.exists(myrt_gff):
        import shutil
        shutil.copy(os.path.abspath(gff_file), myrt_gff)
        logger.info(f"Copied Prodigal GFF to: {myrt_gff}")
        logger.info("  MyRT will use this GFF instead of running FragGeneScan")
    
    # After copying files, add debugging
    logger.info("\n--- Debugging MyRT Input Mismatch ---")
    
    # Check FASTA headers
    logger.info("First 3 FASTA headers in MyRT input:")
    with open(myrt_fna, 'r') as f:
        count = 0
        for line in f:
            if line.startswith('>'):
                logger.info(f"  {line.strip()}")
                count += 1
                if count >= 3:
                    break
    
    # Check GFF contig names
    logger.info("First 3 GFF contig names:")
    with open(myrt_gff, 'r') as f:
        count = 0
        for line in f:
            if not line.startswith('#'):
                contig = line.split('\t')[0]
                logger.info(f"  {contig}")
                count += 1
                if count >= 3:
                    break
    
    # Check if they match
    logger.info("Checking for mismatches...")
    with open(myrt_fna, 'r') as f:
        fasta_headers = [line.split()[0][1:] for line in f if line.startswith('>')]
        
    with open(myrt_gff, 'r') as f:
        gff_contigs = set(line.split('\t')[0] for line in f if not line.startswith('#'))
    
    mismatches = gff_contigs - set(fasta_headers)
    if mismatches:
        logger.error(f"⚠ MISMATCH DETECTED! GFF has contigs not in FASTA: {list(mismatches)[:5]}")
    else:
        logger.info("✓ All GFF contigs found in FASTA")
    logger.info("=" * 80)
    
    # Verify files exist
    if not os.path.exists(myrt_fna) or os.path.getsize(myrt_fna) == 0:
        logger.error(f"MyRT input file missing or empty: {myrt_fna}")
        return None
    
    if not os.path.exists(myrt_gff) or os.path.getsize(myrt_gff) == 0:
        logger.error(f"MyRT GFF file missing or empty: {myrt_gff}")
        return None
    
    logger.info(f"MyRT input FNA size: {os.path.getsize(myrt_fna)} bytes")
    logger.info(f"MyRT input GFF size: {os.path.getsize(myrt_gff)} bytes")
    
    # Run MyRT using rerun.sh
    myrt_script_dir = config.get('myrt_script_dir', '/ibex/user/rioszemm/the-retron-project/src/myRT')
    myrt_cmd = f"cd {myrt_script_dir} && ./rerun.sh {genome_id} {os.path.abspath(myrt_output)}"
    
    logger.info(f"Running command: {myrt_cmd}")
    
    try:
        result = subprocess.run(
            myrt_cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=3600
        )
        
        if result.returncode != 0:
            logger.error(f"MyRT failed: {result.stderr}")
            logger.error(f"MyRT stdout: {result.stdout}")
            return None
        
        logger.info("MyRT completed successfully")
        
        # Debug: Check what files were created
        logger.info("\n--- MyRT Output Files ---")
        myrt_results_dir = os.path.join(myrt_output, genome_id)
        if os.path.exists(myrt_results_dir):
            for file in os.listdir(myrt_results_dir):
                file_path = os.path.join(myrt_results_dir, file)
                if os.path.isfile(file_path):
                    logger.info(f"  {file}: {os.path.getsize(file_path)} bytes")
        
        # Preview main results - MyRT typically creates various output files
        # Look for common MyRT output patterns
        for root, dirs, files in os.walk(myrt_results_dir):
            for file in files:
                if file.endswith('.tsv') or file.endswith('.txt') or 'result' in file.lower():
                    results_file = os.path.join(root, file)
                    logger.info(f"\n--- MyRT Results Preview: {file} (first 5 lines) ---")
                    with open(results_file, 'r') as f:
                        for i, line in enumerate(f):
                            if i < 5:
                                logger.info(f"  {line.rstrip()}")
                            else:
                                break
                    break  # Just show first result file
        
        return myrt_results_dir
        
    except subprocess.TimeoutExpired:
        logger.error("MyRT timed out after 1 hour")
        return None
    except Exception as e:
        logger.error(f"Error running MyRT: {e}")
        return None


def run_defensefinder(faa_file, output_dir, config, logger):
    """
    Run DefenseFinder on pre-generated Prodigal results
    
    Args:
        faa_file: Path to Prodigal-generated .faa file
        output_dir: Output directory
        config: Configuration dict
        logger: Logger instance
    
    Returns:
        Path to DefenseFinder output directory
    """
    logger.info("=" * 80)
    logger.info("Running DefenseFinder for Defense System Detection")
    logger.info("=" * 80)
    
    defensefinder_output = os.path.join(output_dir, "defensefinder_output")
    os.makedirs(defensefinder_output, exist_ok=True)
    
    genome_id = Path(faa_file).stem
    
    # DefenseFinder command
    df_cmd = [
        "defense-finder",
        "run",
        faa_file,
        "--out-dir", defensefinder_output,
        "--workers", str(config.get('threads', 4))
    ]
    
    logger.info(f"Running command: {' '.join(df_cmd)}")
    logger.info(f"Input FAA: {faa_file}")
    
    try:
        result = subprocess.run(
            df_cmd,
            capture_output=True,
            text=True,
            timeout=3600
        )
        
        if result.returncode != 0:
            logger.error(f"DefenseFinder failed: {result.stderr}")
            logger.error(f"DefenseFinder stdout: {result.stdout}")
            return None
        
        logger.info("DefenseFinder completed successfully")
        
        # Debug: Check what files were created
        logger.info("\n--- DefenseFinder Output Files ---")
        for root, dirs, files in os.walk(defensefinder_output):
            for file in files:
                file_path = os.path.join(root, file)
                rel_path = os.path.relpath(file_path, defensefinder_output)
                logger.info(f"  {rel_path}: {os.path.getsize(file_path)} bytes")
        
        # Preview main results - DefenseFinder typically creates defense_finder_systems.tsv
        results_file = None
        for root, dirs, files in os.walk(defensefinder_output):
            for file in files:
                if 'systems' in file.lower() and file.endswith('.tsv'):
                    results_file = os.path.join(root, file)
                    break
        
        if results_file and os.path.exists(results_file):
            logger.info(f"\n--- DefenseFinder Results Preview: {os.path.basename(results_file)} (first 5 lines) ---")
            with open(results_file, 'r') as f:
                for i, line in enumerate(f):
                    if i < 5:
                        logger.info(f"  {line.rstrip()}")
                    else:
                        break
            
            # Count systems
            with open(results_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    logger.info(f"\n  Total defense systems found: {len(lines) - 1}")
                else:
                    logger.info("\n  No defense systems found")
        else:
            logger.warning("Could not find main DefenseFinder results file")
        
        return defensefinder_output
        
    except subprocess.TimeoutExpired:
        logger.error("DefenseFinder timed out after 1 hour")
        return None
    except Exception as e:
        logger.error(f"Error running DefenseFinder: {e}")
        return None


def process_genome(config, logger):
    """
    Main pipeline orchestration - runs all tools with shared Prodigal results
    Includes checkpoint system to skip completed steps
    
    Args:
        config: Configuration dictionary
        logger: Logger instance
    """

    # Create output directory
    os.makedirs(config['output_dir'], exist_ok=True)
    
    # Initialize results dictionary
    results = {
        'prodigal': {'faa': None, 'gff': None, 'mapping': None},
        'infernal': {'table': None, 'sequences': None},
        'myrt': None,
        'padloc': None,
        'defensefinder': None,
    }

    from pathlib import Path
    
    genome_id = Path(config['genome_fasta']).stem

    # =========================================================================
    # STEP 0: Gene Prediction with Prodigal
    # =========================================================================
    logger.info("\n" + "=" * 80)
    logger.info("STEP 0: Gene Prediction with Prodigal")
    logger.info("=" * 80)

    # --- Paths (use Path consistently) ---
    output_dir = Path(config['output_dir'])
    prodigal_dir = output_dir / "prodigal_results"

    expected_gff = prodigal_dir / "genes.gff"
    expected_faa = prodigal_dir / "proteins_raw.faa"
    sanitized_genome = prodigal_dir / "genome_sanitized.fasta"

    # --- Checkpoint: skip Prodigal if outputs exist ---
    if (
        expected_gff.exists() and expected_faa.exists() and
        expected_gff.stat().st_size > 0 and expected_faa.stat().st_size > 0
    ):
        logger.info("✓ Prodigal outputs already exist - skipping")
        logger.info(f"  GFF: {expected_gff} ({expected_gff.stat().st_size} bytes)")
        logger.info(f"  FAA: {expected_faa} ({expected_faa.stat().st_size} bytes)")

        gff_file = str(expected_gff)
        faa_file = str(expected_faa)

        # Determine which genome to use for downstream tools
        if sanitized_genome.exists() and sanitized_genome.stat().st_size > 0:
            genome_fasta_to_use = str(sanitized_genome)
            logger.info(f"  Using sanitized genome: {genome_fasta_to_use}")
        else:
            genome_fasta_to_use = config['genome_fasta']
            logger.info(f"  Using original genome: {genome_fasta_to_use}")

    else:
        logger.info("Running Prodigal for gene prediction...")
        try:
            faa_file, gff_file, genome_fasta_to_use = run_prodigal_once(
                input_fasta=config['genome_fasta'],
                output_dir=str(output_dir),
                prodigal_bin=config.get('prodigal_bin', 'prodigal')
            )

            if not faa_file or not gff_file:
                raise RuntimeError("Prodigal did not generate expected output files")

            if not os.path.exists(faa_file) or not os.path.exists(gff_file):
                raise RuntimeError(
                    f"Prodigal output files not found at: {faa_file}, {gff_file}"
                )

            logger.info("✓ Prodigal completed successfully")
            logger.info(f"  GFF: {gff_file}")
            logger.info(f"  FAA: {faa_file}")
            
            if genome_fasta_to_use != config['genome_fasta']:
                logger.info(f"  Sanitized genome: {genome_fasta_to_use}")

        except Exception as e:
            logger.error(f"✗ Prodigal failed: {e}")
            logger.error("Cannot proceed without gene predictions")
            raise RuntimeError(f"Prodigal gene prediction failed: {e}")

    # --- Store results ---
    results['prodigal']['faa'] = faa_file
    results['prodigal']['gff'] = gff_file
    results['prodigal']['mapping'] = None  # Not used in this version
    results['prodigal']['genome_used'] = genome_fasta_to_use


    # =========================================================================
    # STEP 0.5: ncRNA Detection with Infernal (Optional)
    # =========================================================================
    if config.get('run_infernal', False):
        logger.info("\n" + "=" * 80)
        logger.info("STEP 0.5: ncRNA Detection with Infernal")
        logger.info("=" * 80)
        
        # Define expected Infernal outputs
        infernal_dir = os.path.join(config['output_dir'], 'infernal_results')
        expected_ncrna = os.path.join(infernal_dir, f"{genome_id}_ncrna.tblout")
        expected_sequences = os.path.join(infernal_dir, f"{genome_id}_ncrna_sequences_all.fa")
        expected_summary = os.path.join(infernal_dir, f"{genome_id}_ncrna_summary.txt")
        
        # Check if required Infernal output exists
        if os.path.exists(expected_ncrna) and os.path.getsize(expected_ncrna) > 0:
            logger.info(f"✓ Infernal outputs already exist - skipping")
            logger.info(f"  Directory: {infernal_dir}")
            logger.info(f"  Table: {os.path.basename(expected_ncrna)}")
            ncrna_file = expected_ncrna
            
            if os.path.exists(expected_sequences):
                logger.info(f"  Sequences: {os.path.basename(expected_sequences)}")
                ncrna_sequences = expected_sequences
            else:
                ncrna_sequences = None
            
            if os.path.exists(expected_summary):
                logger.info(f"  Summary: {os.path.basename(expected_summary)}")
                with open(expected_summary) as f:
                    summary_lines = f.readlines()
                    for line in summary_lines:
                        if line.startswith('Total hits:') or line.startswith('High-confidence'):
                            logger.info(f"    {line.strip()}")
            
            results['infernal']['table'] = ncrna_file
            results['infernal']['sequences'] = ncrna_sequences
        else:
            logger.info("Running Infernal for ncRNA detection...")
            try:
                ncrna_file, ncrna_sequences = run_infernal(
                    # genome_fasta=config['genome_fasta'],
                    genome_fasta_to_use,
                    output_dir=config['output_dir'],
                    cpu=config.get('threads', 4),
                    extract_sequences=True
                )
                
                if ncrna_file:
                    logger.info(f"✓ Infernal completed successfully")
                    logger.info(f"  Directory: {infernal_dir}")
                    logger.info(f"  Table: {os.path.basename(ncrna_file)}")
                    
                    if ncrna_sequences:
                        logger.info(f"  Sequences: {os.path.basename(ncrna_sequences)}")
                    
                    if os.path.exists(expected_summary):
                        logger.info(f"  Summary: {os.path.basename(expected_summary)}")
                        with open(expected_summary) as f:
                            summary_lines = f.readlines()
                            for line in summary_lines:
                                if line.startswith('Total hits:') or line.startswith('High-confidence'):
                                    logger.info(f"    {line.strip()}")
                    
                    results['infernal']['table'] = ncrna_file
                    results['infernal']['sequences'] = ncrna_sequences
                else:
                    logger.warning("✗ Infernal did not produce expected outputs")
                    
            except Exception as e:
                logger.warning(f"✗ Infernal failed: {e}")
                logger.warning("  Continuing without ncRNA data...")
                import traceback
                traceback.print_exc()
    else:
        logger.info("\n" + "=" * 80)
        logger.info("STEP 0.5: ncRNA Detection - SKIPPED (disabled in config)")
        logger.info("=" * 80)

    # =========================================================================
    # STEP 1: RT Detection with MyRT
    # =========================================================================
    logger.info("\n" + "=" * 80)
    logger.info("STEP 1: RT Detection with MyRT")
    logger.info("=" * 80)
    
    # # Define expected MyRT outputs (BOTH required files)
    # myrt_dir = os.path.join(config['output_dir'], 'myrt_output', genome_id)
    # expected_rvt_domtblout = os.path.join(myrt_dir, f"{genome_id}-RVT.domtblout")
    # expected_gn_list = os.path.join(myrt_dir, f"{genome_id}-gn.list")


    # Define expected MyRT outputs (BOTH required files)
    # Use the actual genome name that was passed to MyRT
    # myrt_genome_id = Path(genome_fasta_to_use).stem
    # myrt_dir = os.path.join(config['output_dir'], 'myrt_output', myrt_genome_id)
    # expected_rvt_domtblout = os.path.join(myrt_dir, f"{myrt_genome_id}-RVT.domtblout")
    # expected_gn_list = os.path.join(myrt_dir, f"{myrt_genome_id}-gn.list")



    myrt_genome_id = Path(genome_fasta_to_use).stem
    myrt_dir = os.path.join(config['output_dir'], 'myrt_output')  # <-- No subdirectory
    expected_rvt_domtblout = os.path.join(myrt_dir, f"{myrt_genome_id}-RVT.domtblout")
    expected_gn_list = os.path.join(myrt_dir, f"{myrt_genome_id}-gn.list")

    
    # Check if BOTH required MyRT outputs exist
    if (os.path.exists(expected_rvt_domtblout) and os.path.exists(expected_gn_list) and
        os.path.getsize(expected_rvt_domtblout) > 0 and os.path.getsize(expected_gn_list) > 0):
        
        logger.info(f"✓ MyRT outputs already exist - skipping")
        logger.info(f"  Directory: {myrt_dir}")
        logger.info(f"  RVT domtblout: {os.path.basename(expected_rvt_domtblout)}")
        logger.info(f"  Gene list: {os.path.basename(expected_gn_list)}")
        results['myrt'] = myrt_dir
    else:
        logger.info("Running MyRT...")
        try:
            myrt_results = run_myrt(
                # config['genome_fasta'],
                genome_fasta_to_use,
                gff_file,
                config['output_dir'],
                config,
                logger
            )
            
            if myrt_results and os.path.exists(myrt_results):
                # Verify the required files were created
                if (os.path.exists(expected_rvt_domtblout) and 
                    os.path.exists(expected_gn_list)):
                    logger.info(f"✓ MyRT completed successfully")
                    logger.info(f"  Results: {myrt_results}")
                    results['myrt'] = myrt_results
                else:
                    logger.warning("✗ MyRT completed but required files missing")
                    logger.warning(f"  Expected: {expected_rvt_domtblout}")
                    logger.warning(f"  Expected: {expected_gn_list}")
            else:
                logger.warning("✗ MyRT did not produce expected outputs")
                
        except subprocess.CalledProcessError as e:
            logger.warning(f"✗ MyRT failed with return code {e.returncode}")
            logger.warning(f"  Command: {' '.join(e.cmd) if hasattr(e, 'cmd') else 'N/A'}")
            logger.warning("  Continuing with other tools...")
        except Exception as e:
            logger.warning(f"✗ MyRT failed: {e}")
            logger.warning("  Continuing with other tools...")

    # =========================================================================
    # STEP 2: Defense System Detection with PADLOC
    # =========================================================================
    logger.info("\n" + "=" * 80)
    logger.info("STEP 2: Defense System Detection with PADLOC")
    logger.info("=" * 80)

    # Define expected PADLOC output
    padloc_dir = os.path.join(config['output_dir'], 'padloc_results')
    expected_padloc_csv = os.path.join(padloc_dir, "proteins_padloc.csv")

    # Check if required PADLOC output exists
    if os.path.exists(expected_padloc_csv) and os.path.getsize(expected_padloc_csv) > 0:
        logger.info(f"✓ PADLOC outputs already exist - skipping")
        logger.info(f"  Results: {expected_padloc_csv}")
        results['padloc'] = expected_padloc_csv
    else:
        logger.info("Running PADLOC...")
        
        # Get ncRNA file from results if it exists
        ncrna_file_to_use = results.get('infernal', {}).get('table', None)
        
        if ncrna_file_to_use and os.path.exists(ncrna_file_to_use):
            logger.info(f"  Including ncRNA data from Infernal")
            logger.info(f"  ncRNA file: {os.path.basename(ncrna_file_to_use)}")
        else:
            logger.info(f"  No ncRNA data available")
            ncrna_file_to_use = None
        
        try:
            padloc_results = run_padloc(
                output_dir=config['output_dir'],
                gff_file=gff_file,
                faa_file=faa_file,
                ncrna_file=ncrna_file_to_use,
                cpu=config.get('threads', 4)
            )
            
            if padloc_results and os.path.exists(padloc_results):
                logger.info(f"✓ PADLOC completed successfully")
                logger.info(f"  Results: {padloc_results}")
                results['padloc'] = padloc_results
            else:
                logger.warning("✗ PADLOC did not produce CSV output")
                
        except Exception as e:
            logger.warning(f"✗ PADLOC failed: {e}")
            logger.warning("  Continuing with other tools...")
            import traceback
            traceback.print_exc()
    
    # =========================================================================
    # STEP 3: Defense System Detection with DefenseFinder
    # =========================================================================
    logger.info("\n" + "=" * 80)
    logger.info("STEP 3: Defense System Detection with DefenseFinder")
    logger.info("=" * 80)
    
    # Define expected DefenseFinder output
    df_dir = os.path.join(config['output_dir'], 'defensefinder_output')
    expected_df_systems = os.path.join(df_dir, "proteins_defense_finder_systems.tsv")
    
    # Check if required DefenseFinder output exists
    if os.path.exists(expected_df_systems) and os.path.getsize(expected_df_systems) > 0:
        logger.info(f"✓ DefenseFinder outputs already exist - skipping")
        logger.info(f"  Directory: {df_dir}")
        logger.info(f"  Systems: {os.path.basename(expected_df_systems)}")
        results['defensefinder'] = df_dir
    else:
        logger.info("Running DefenseFinder...")
        try:
            defensefinder_results = run_defensefinder(
                faa_file,
                config['output_dir'],
                config,
                logger
            )
            
            if defensefinder_results and os.path.exists(defensefinder_results):
                # Verify the required file was created
                if os.path.exists(expected_df_systems):
                    logger.info(f"✓ DefenseFinder completed successfully")
                    logger.info(f"  Results: {defensefinder_results}")
                    results['defensefinder'] = defensefinder_results
                else:
                    logger.warning("✗ DefenseFinder completed but systems file missing")
                    logger.warning(f"  Expected: {expected_df_systems}")
            else:
                logger.warning("✗ DefenseFinder did not produce expected outputs")
                
        except subprocess.CalledProcessError as e:
            logger.warning(f"✗ DefenseFinder failed with return code {e.returncode}")
            logger.warning(f"  Command: {' '.join(e.cmd) if hasattr(e, 'cmd') else 'N/A'}")
            logger.warning("  Continuing with other tools...")
        except Exception as e:
            logger.warning(f"✗ DefenseFinder failed: {e}")
            logger.warning("  Continuing with other tools...")
    
    # =========================================================================
    # Final Summary
    # =========================================================================
    logger.info("\n" + "=" * 80)
    logger.info("Pipeline Execution Summary")
    logger.info("=" * 80)
    logger.info("\nAll raw output files are preserved in their respective directories:")
    
    # Prodigal
    if results['prodigal']['faa'] and results['prodigal']['gff']:
        logger.info(f"  ✓ Prodigal results:")
        logger.info(f"    - GFF: {os.path.basename(results['prodigal']['gff'])}")
        logger.info(f"    - FAA: {os.path.basename(results['prodigal']['faa'])}")
    else:
        logger.info(f"  ✗ Prodigal: Failed")
    
    # Infernal
    if results['infernal']['table']:
        logger.info(f"  ✓ Infernal results:")
        logger.info(f"    - Table: {os.path.basename(results['infernal']['table'])}")
        if results['infernal']['sequences']:
            logger.info(f"    - Sequences: {os.path.basename(results['infernal']['sequences'])}")
    elif config.get('run_infernal', False):
        logger.info(f"  ✗ Infernal: Failed")
    else:
        logger.info(f"  ⊘ Infernal: Skipped (disabled)")
    
    # MyRT
    if results['myrt']:
        logger.info(f"  ✓ MyRT results: {results['myrt']}")
    else:
        logger.info(f"  ✗ MyRT: Failed or skipped")
    
    # PADLOC
    if results['padloc']:
        logger.info(f"  ✓ PADLOC results: {results['padloc']}")
    else:
        logger.info(f"  ✗ PADLOC: Failed or skipped")
    
    # DefenseFinder
    if results['defensefinder']:
        logger.info(f"  ✓ DefenseFinder results: {results['defensefinder']}")
    else:
        logger.info(f"  ✗ DefenseFinder: Failed or skipped")
    
    # Count successes
    successful_tools = sum([
        bool(results['prodigal']['faa']),
        bool(results['myrt']),
        bool(results['padloc']),
        bool(results['defensefinder'])
    ])
    
    if config.get('run_infernal', False) and results['infernal']['table']:
        successful_tools += 1
    
    total_tools = 5 if config.get('run_infernal', False) else 4
    
    logger.info(f"\n📊 Successfully completed: {successful_tools}/{total_tools} tools")
    logger.info("✅ No integration or summarization performed")
    logger.info("✅ Analyze raw tool outputs directly for maximum detail")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive defense systems and retron annotation pipeline'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input genome FASTA file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of CPU threads to use (default: 4)'
    )
    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-5,
        help='E-value cutoff for MyRT (default: 1e-5)'
    )
    parser.add_argument(
        '--run-infernal',
        action='store_true',
        help='Run Infernal for ncRNA detection (improves retron classification)'
    )
    parser.add_argument(
        '--infernal-cm-db',
        default='/ibex/user/rioszemm/the-retron-project/src/padloc/data/cm/padlocdb.cm',
        help='Path to Infernal CM database (default: PADLOC CM database)'
    )
    parser.add_argument(
        '--log',
        help='Log file path (optional)'
    )
    parser.add_argument(
        '--myrt-script-dir',
        default='/ibex/user/rioszemm/the-retron-project/src/myRT',
        help='Path to MyRT script directory'
    )
    
    args = parser.parse_args()
    
    # Setup logger
    logger = setup_logger(args.log)
    
    # Handle gzipped input - decompress once at the start
    from pathlib import Path
    import gzip
    import shutil
    
    input_path = Path(args.input)
    actual_input = args.input
    
    if input_path.suffix == '.gz' or is_gzipped(input_path):
        logger.info("=" * 80)
        logger.info("Input file is gzipped - decompressing...")
        logger.info("=" * 80)
        
        decompressed = input_path.with_suffix('')
        
        if not decompressed.exists():
            with gzip.open(input_path, 'rb') as f_in:
                with open(decompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"Decompressed: {input_path} → {decompressed}")
        else:
            logger.info(f"Decompressed file already exists: {decompressed}")
        
        actual_input = str(decompressed)
    
    # Setup config dictionary with decompressed path
    config = {
        'genome_fasta': actual_input,  # Use decompressed path
        'output_dir': args.output,
        'evalue_cutoff': args.evalue,
        'threads': args.threads,
        'myrt_script_dir': args.myrt_script_dir,
        'run_infernal': args.run_infernal,
        'infernal_cm_db': args.infernal_cm_db,
    }
    
    try:
        # Run pipeline
        logger.info("=" * 80)
        logger.info("Starting Defense Systems Analysis Pipeline")
        logger.info("=" * 80)
        logger.info(f"Input genome: {actual_input}")
        logger.info(f"Output directory: {args.output}")
        logger.info(f"Threads: {args.threads}")
        logger.info(f"Infernal ncRNA detection: {'Enabled' if args.run_infernal else 'Disabled'}")
        
        process_genome(config, logger)
        
        logger.info("=" * 80)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


def is_gzipped(filepath):
    """Check if a file is gzipped by reading its magic number"""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


if __name__ == '__main__':
    main()





# =============================================================================
# ENTRY POINT
# =============================================================================

# def main():
#     parser = argparse.ArgumentParser(
#         description='Comprehensive defense systems and retron annotation pipeline'
#     )
#     parser.add_argument(
#         '-i', '--input',
#         required=True,
#         help='Input genome FASTA file'
#     )
#     parser.add_argument(
#         '-o', '--output',
#         required=True,
#         help='Output directory'
#     )
#     parser.add_argument(
#         '--threads',
#         type=int,
#         default=4,
#         help='Number of CPU threads to use (default: 4)'
#     )
#     parser.add_argument(
#         '--evalue',
#         type=float,
#         default=1e-5,
#         help='E-value cutoff for MyRT (default: 1e-5)'
#     )
#     parser.add_argument(
#         '--run-infernal',  # ← NEW
#         action='store_true',
#         help='Run Infernal for ncRNA detection (improves retron classification)'
#     )
#     parser.add_argument(
#         '--infernal-cm-db',  # ← NEW
#         default='/ibex/user/rioszemm/the-retron-project/src/padloc/data/cm/padlocdb.cm',
#         help='Path to Infernal CM database (default: PADLOC CM database)'
#     )
#     parser.add_argument(
#         '--log',
#         help='Log file path (optional)'
#     )
#     parser.add_argument(
#         '--myrt-script-dir',
#         default='/ibex/user/rioszemm/the-retron-project/src/myRT',
#         help='Path to MyRT script directory'
#     )
    
#     args = parser.parse_args()
    
#     # Setup logger
#     logger = setup_logger(args.log)
    
#     # Setup config dictionary
#     config = {
#         'genome_fasta': args.input,
#         'output_dir': args.output,
#         'evalue_cutoff': args.evalue,
#         'threads': args.threads,
#         'myrt_script_dir': args.myrt_script_dir,
#         'run_infernal': args.run_infernal,  # ← NEW
#         'infernal_cm_db': args.infernal_cm_db,  # ← NEW
#     }
    
#     try:
#         # Run pipeline
#         logger.info("=" * 80)
#         logger.info("Starting Defense Systems Analysis Pipeline")
#         logger.info("=" * 80)
#         logger.info(f"Input genome: {args.input}")
#         logger.info(f"Output directory: {args.output}")
#         logger.info(f"Threads: {args.threads}")
#         # logger.info(f"MyRT E-value cutoff: {args.evalue}")
#         logger.info(f"Infernal ncRNA detection: {'Enabled' if args.run_infernal else 'Disabled'}")  # ← NEW
        
#         process_genome(config, logger)
        
#         logger.info("=" * 80)
#         logger.info("Pipeline completed successfully!")
#         logger.info("=" * 80)
        
#     except Exception as e:
#         logger.error(f"Pipeline failed: {e}", exc_info=True)
#         sys.exit(1)

# if __name__ == '__main__':
#     main()

