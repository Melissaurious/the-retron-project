#!/usr/bin/env python3
"""
Batch Splitter - Split large genome metadata files into manageable batches.

Features:
- Splits large JSON files (millions of genomes) into batches
- Optional quality-based sorting (high-quality genomes first)
- Creates manifest for tracking
- Generates SLURM job scripts for parallel processing

Usage:
    # Basic splitting
    python batch_splitter.py \\
        --input gtdb_bacteria_parsed.json \\
        --output batches/gtdb_bacteria \\
        --batch-size 10000

    # With quality sorting and SLURM jobs
    python batch_splitter.py \\
        --input gtdb_bacteria_parsed.json \\
        --output batches/gtdb_bacteria \\
        --batch-size 10000 \\
        --sort-by-quality \\
        --create-slurm-jobs \\
        --slurm-time 24:00:00 \\
        --slurm-mem 16G
"""

import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
import sys

def calculate_quality_score(genome_data: Dict) -> float:
    """
    Calculate quality score for genome.
    
    Uses GTDB standard formula: completeness - (5 × contamination)
    
    Quality Tiers (GTDB/MIMAG standards):
    - High quality: ≥90 (completeness >95%, contamination <1%)
    - Medium quality: ≥50 (completeness >70%, contamination <4%)
    - Low quality: <50 (may be excluded from GTDB)
    
    Reference: Parks et al. (2020) Nature Biotechnology
    https://gtdb.ecogenomic.org/methods
    
    Args:
        genome_data: Genome metadata dictionary
    
    Returns:
        Quality score (-50 to 100)
    """
    # Try different field names for completeness
    completeness = (
        genome_data.get('Completeness') or
        genome_data.get('completeness') or
        genome_data.get('checkm_completeness') or
        genome_data.get('checkm2_completeness') or
        0
    )
    
    # Try different field names for contamination
    contamination = (
        genome_data.get('Contamination') or
        genome_data.get('contamination') or
        genome_data.get('checkm_contamination') or
        genome_data.get('checkm2_contamination') or
        0
    )
    
    # GTDB standard quality formula
    return float(completeness) - (5 * float(contamination))


def sort_genomes_by_quality(genomes: Dict) -> List[Tuple[str, Dict]]:
    """
    Sort genomes by quality score (high to low).
    
    Args:
        genomes: Dictionary of genome_id -> metadata
    
    Returns:
        List of (genome_id, metadata) tuples, sorted by quality
    """
    genome_list = [(gid, data) for gid, data in genomes.items()]
    
    # Sort by quality score (descending)
    genome_list.sort(key=lambda x: calculate_quality_score(x[1]), reverse=True)
    
    return genome_list


def create_batches(
    genomes: Dict,
    batch_size: int,
    sort_by_quality: bool = False
) -> List[Dict]:
    """
    Split genomes into batches.
    
    Args:
        genomes: Dictionary of genome metadata
        batch_size: Number of genomes per batch
        sort_by_quality: Whether to sort by quality first
    
    Returns:
        List of batch dictionaries
    """
    # Sort if requested
    if sort_by_quality:
        print(f"Sorting {len(genomes):,} genomes by quality...")
        genome_list = sort_genomes_by_quality(genomes)
    else:
        genome_list = list(genomes.items())
    
    # Split into batches
    batches = []
    for i in range(0, len(genome_list), batch_size):
        batch_genomes = genome_list[i:i + batch_size]
        batch = {gid: data for gid, data in batch_genomes}
        batches.append(batch)
    
    return batches


def save_batches(
    batches: List[Dict],
    output_dir: Path,
    database_name: str
) -> List[Dict]:
    """
    Save batches to disk.
    
    Args:
        batches: List of batch dictionaries
        output_dir: Output directory
        database_name: Database name
    
    Returns:
        List of batch info dictionaries for manifest
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    batch_info_list = []
    
    for batch_id, batch_genomes in enumerate(batches):
        batch_file = f"batch_{batch_id:04d}.json"
        batch_path = output_dir / batch_file
        
        # Save batch
        with open(batch_path, 'w') as f:
            json.dump(batch_genomes, f, indent=2)
        
        # Record info
        batch_info = {
            'batch_id': batch_id,
            'file': batch_file,
            'genome_count': len(batch_genomes),
            'size_bytes': batch_path.stat().st_size
        }
        batch_info_list.append(batch_info)
        
        print(f"  Created {batch_file} ({len(batch_genomes):,} genomes, {batch_info['size_bytes'] / 1024 / 1024:.1f} MB)")
    
    return batch_info_list


def create_manifest(
    database_name: str,
    total_genomes: int,
    batch_info_list: List[Dict],
    output_dir: Path,
    metadata: Dict = None
) -> Path:
    """
    Create manifest file.
    
    Args:
        database_name: Database name
        total_genomes: Total number of genomes
        batch_info_list: List of batch info dictionaries
        output_dir: Output directory
        metadata: Additional metadata to include
    
    Returns:
        Path to manifest file
    """
    manifest = {
        'database': database_name,
        'total_genomes': total_genomes,
        'total_batches': len(batch_info_list),
        'batch_size': batch_info_list[0]['genome_count'] if batch_info_list else 0,
        'created': datetime.now().isoformat(),
        'batches': batch_info_list
    }
    
    if metadata:
        manifest['metadata'] = metadata
    
    manifest_path = output_dir / 'manifest.json'
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"\nCreated manifest: {manifest_path}")
    return manifest_path


def create_slurm_job_scripts(
    database_name: str,
    batch_info_list: List[Dict],
    output_dir: Path,
    downloads_dir: Path,
    time_limit: str = "24:00:00",
    memory: str = "16G",
    cpus: int = 4,
    partition: str = "batch",
    email: str = None
) -> Path:
    """
    Create SLURM job scripts for parallel batch processing.
    
    Args:
        database_name: Database name
        batch_info_list: List of batch info
        output_dir: Batch directory
        downloads_dir: Where to save downloads
        time_limit: SLURM time limit
        memory: Memory allocation
        cpus: CPUs per task
        partition: SLURM partition
        email: Email for notifications
    
    Returns:
        Path to main submit script
    """
    slurm_dir = output_dir / 'slurm_jobs'
    slurm_dir.mkdir(parents=True, exist_ok=True)
    
    # Create individual job scripts for each batch
    job_scripts = []
    
    for batch_info in batch_info_list:
        batch_id = batch_info['batch_id']
        job_name = f"{database_name}_batch_{batch_id:04d}"
        job_file = slurm_dir / f"job_batch_{batch_id:04d}.sh"
        
        job_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={slurm_dir}/batch_{batch_id:04d}_%j.out
#SBATCH --error={slurm_dir}/batch_{batch_id:04d}_%j.err
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --partition={partition}
"""
        
        if email:
            job_script += f"#SBATCH --mail-user={email}\n"
            job_script += "#SBATCH --mail-type=END,FAIL\n"
        
        job_script += f"""
# Job info
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Database: {database_name}"
echo "Batch: {batch_id:04d}"
echo "Start time: $(date)"

# Activate conda environment (adjust path as needed)
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate retron_tradicional

# Run batch processor
python batch_processor.py \\
    --database {database_name} \\
    --batch-dir {output_dir.absolute()} \\
    --output-dir {downloads_dir.absolute()} \\
    --batch-id {batch_id} \\
    --single-batch

# Check exit status
if [ $? -eq 0 ]; then
    echo "Batch {batch_id:04d} completed successfully"
else
    echo "Batch {batch_id:04d} failed with exit code $?"
    exit 1
fi

echo "End time: $(date)"
"""
        
        with open(job_file, 'w') as f:
            f.write(job_script)
        
        job_file.chmod(0o755)  # Make executable
        job_scripts.append(job_file)
    
    # Create array job script (more efficient)
    array_job_file = slurm_dir / f"submit_array_{database_name}.sh"
    
    array_script = f"""#!/bin/bash
#SBATCH --job-name={database_name}_array
#SBATCH --output={slurm_dir}/batch_%a_%j.out
#SBATCH --error={slurm_dir}/batch_%a_%j.err
#SBATCH --array=0-{len(batch_info_list)-1}
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --partition={partition}
"""
    
    if email:
        array_script += f"#SBATCH --mail-user={email}\n"
        array_script += "#SBATCH --mail-type=END,FAIL\n"
    
    array_script += f"""
# Job info
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Database: {database_name}"
echo "Start time: $(date)"

# Activate conda environment (adjust path as needed)
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate retron_tradicional

# Run batch processor for this array task
python batch_processor.py \\
    --database {database_name} \\
    --batch-dir {output_dir.absolute()} \\
    --output-dir {downloads_dir.absolute()} \\
    --batch-id $SLURM_ARRAY_TASK_ID \\
    --single-batch

# Check exit status
if [ $? -eq 0 ]; then
    echo "Batch $SLURM_ARRAY_TASK_ID completed successfully"
else
    echo "Batch $SLURM_ARRAY_TASK_ID failed"
    exit 1
fi

echo "End time: $(date)"
"""
    
    with open(array_job_file, 'w') as f:
        f.write(array_script)
    
    array_job_file.chmod(0o755)
    
    # Create main submit script
    submit_script_file = slurm_dir / f"submit_{database_name}.sh"
    
    submit_script = f"""#!/bin/bash
# Submit all {database_name} batch jobs
#
# Usage:
#   bash {submit_script_file.name}
#
# Or submit as array job (recommended):
#   sbatch {array_job_file.name}

echo "Submitting {len(batch_info_list)} batch jobs for {database_name}..."

# Option 1: Submit as array job (RECOMMENDED - more efficient)
echo "Recommended: sbatch {array_job_file.name}"

# Option 2: Submit individual jobs (use if array jobs not available)
# Uncomment below to use:
"""
    
    for job_file in job_scripts:
        submit_script += f"# sbatch {job_file.name}\n"
    
    submit_script += f"""
echo "Done! Monitor with: squeue -u $USER"
echo "Check progress: ls {slurm_dir}/*.out | wc -l"
"""
    
    with open(submit_script_file, 'w') as f:
        f.write(submit_script)
    
    submit_script_file.chmod(0o755)
    
    print(f"\n✅ Created SLURM job scripts:")
    print(f"   Array job (recommended): {array_job_file}")
    print(f"   Individual jobs: {len(job_scripts)} scripts in {slurm_dir}/")
    print(f"   Submit script: {submit_script_file}")
    print(f"\nTo submit:")
    print(f"   sbatch {array_job_file}")
    print(f"   # or")
    print(f"   bash {submit_script_file}")
    
    return submit_script_file


def main():
    parser = argparse.ArgumentParser(
        description='Split large genome metadata into manageable batches',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic splitting
  python batch_splitter.py \\
      --input gtdb_bacteria_parsed.json \\
      --output batches/gtdb_bacteria \\
      --batch-size 10000

  # With quality sorting and SLURM jobs
  python batch_splitter.py \\
      --input gtdb_bacteria_parsed.json \\
      --output batches/gtdb_bacteria \\
      --batch-size 10000 \\
      --sort-by-quality \\
      --create-slurm-jobs \\
      --downloads-dir /scratch/user/downloads
        """
    )
    
    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input JSON file (e.g., gtdb_bacteria_parsed.json)'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output directory for batches'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=10000,
        help='Number of genomes per batch (default: 10000)'
    )
    
    parser.add_argument(
        '--database',
        help='Database name (auto-detected from filename if not provided)'
    )
    
    parser.add_argument(
        '--sort-by-quality',
        action='store_true',
        help='Sort genomes by quality score before batching (high quality first)'
    )
    
    # SLURM options
    parser.add_argument(
        '--create-slurm-jobs',
        action='store_true',
        help='Create SLURM job scripts for parallel processing'
    )
    
    parser.add_argument(
        '--downloads-dir',
        type=Path,
        default=Path('./downloads'),
        help='Directory for downloads (used in SLURM scripts, default: ./downloads)'
    )
    
    parser.add_argument(
        '--slurm-time',
        default='24:00:00',
        help='SLURM time limit (default: 24:00:00)'
    )
    
    parser.add_argument(
        '--slurm-mem',
        default='16G',
        help='SLURM memory allocation (default: 16G)'
    )
    
    parser.add_argument(
        '--slurm-cpus',
        type=int,
        default=4,
        help='SLURM CPUs per task (default: 4)'
    )
    
    parser.add_argument(
        '--slurm-partition',
        default='batch',
        help='SLURM partition (default: batch)'
    )
    
    parser.add_argument(
        '--slurm-email',
        help='Email for SLURM notifications'
    )
    
    args = parser.parse_args()
    
    # Validate input
    if not args.input.exists():
        print(f"❌ Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Auto-detect database name from filename
    if not args.database:
        # gtdb_bacteria_parsed.json -> gtdb_bacteria
        database_name = args.input.stem.replace('_parsed', '')
    else:
        database_name = args.database
    
    print(f"{'='*60}")
    print(f"Batch Splitter")
    print(f"{'='*60}")
    print(f"Database: {database_name}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Batch size: {args.batch_size:,}")
    print(f"Sort by quality: {args.sort_by_quality}")
    print(f"Create SLURM jobs: {args.create_slurm_jobs}")
    print(f"{'='*60}\n")
    
    # Load input JSON
    print(f"Loading {args.input}...")
    with open(args.input) as f:
        genomes = json.load(f)
    
    total_genomes = len(genomes)
    print(f"Loaded {total_genomes:,} genomes\n")
    
    # Create batches
    print(f"Creating batches...")
    batches = create_batches(genomes, args.batch_size, args.sort_by_quality)
    print(f"Created {len(batches)} batches\n")
    
    # Save batches
    print(f"Saving batches to {args.output}...")
    batch_info_list = save_batches(batches, args.output, database_name)
    
    # Create manifest
    create_manifest(
        database_name=database_name,
        total_genomes=total_genomes,
        batch_info_list=batch_info_list,
        output_dir=args.output,
        metadata={
            'source_file': str(args.input),
            'batch_size': args.batch_size,
            'sorted_by_quality': args.sort_by_quality
        }
    )
    
    # Create SLURM jobs if requested
    if args.create_slurm_jobs:
        print(f"\nCreating SLURM job scripts...")
        create_slurm_job_scripts(
            database_name=database_name,
            batch_info_list=batch_info_list,
            output_dir=args.output,
            downloads_dir=args.downloads_dir,
            time_limit=args.slurm_time,
            memory=args.slurm_mem,
            cpus=args.slurm_cpus,
            partition=args.slurm_partition,
            email=args.slurm_email
        )
    
    # Summary
    print(f"\n{'='*60}")
    print(f"✅ Batch creation complete!")
    print(f"{'='*60}")
    print(f"Total genomes: {total_genomes:,}")
    print(f"Total batches: {len(batches)}")
    print(f"Batches directory: {args.output}")
    print(f"\nNext steps:")
    if args.create_slurm_jobs:
        print(f"  1. Review SLURM scripts in {args.output}/slurm_jobs/")
        print(f"  2. Submit jobs: sbatch {args.output}/slurm_jobs/submit_array_{database_name}.sh")
    else:
        print(f"  1. Run: python batch_processor.py --database {database_name} --batch-dir {args.output}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
