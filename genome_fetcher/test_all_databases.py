#!/usr/bin/env python3
"""
Simple Test Batch Generator
Extracts first N genomes from each database JSON for testing.
"""

import json
import argparse
from pathlib import Path

def create_test_batch(input_file: Path, output_file: Path, num_genomes: int = 10):
    """Extract first N genomes from a database JSON file."""
    
    print(f"Reading: {input_file}")
    
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    # Handle both dict and list formats
    if isinstance(data, dict):
        # Take first N items
        items = list(data.items())[:num_genomes]
        batch_data = {k: v for k, v in items}
    elif isinstance(data, list):
        batch_data = data[:num_genomes]
    else:
        raise ValueError(f"Unexpected data format in {input_file}")
    
    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Write batch file
    with open(output_file, 'w') as f:
        json.dump(batch_data, f, indent=2)
    
    print(f"  Created batch with {len(batch_data)} genomes")
    print(f"  Output: {output_file}")
    
    return len(batch_data)

def main():
    parser = argparse.ArgumentParser(description='Create test batches from database JSONs')
    # parser.add_argument('metadata-dir', required=True, help='Directory with parsed JSON files')
    # parser.add_argument('--output-dir', required=True, help='Output directory for test batches')
    # parser.add_argument('--num-genomes', type=int, default=10, help='Number of genomes per batch')

    parser.add_argument('--metadata-dir', default="/ibex/user/rioszemm/the-retron-project/huge_databases/second_try/parsed_metadafa", help='Directory with parsed JSON files')
    parser.add_argument('--output-dir', default="/ibex/user/rioszemm/the-retron-project/download_data_/genome_fetcher", help='Output directory for test batches')
    parser.add_argument('--num-genomes', type=int, default=10, help='Number of genomes per batch')


    
    args = parser.parse_args()
    
    metadata_dir = Path(args.metadata_dir)
    output_dir = Path(args.output_dir)
    
    databases = [
        'gtdb_bacteria',
        'gtdb_archaea',
        'ncbi_bacteria',
        'ncbi_archaea',
        'mgnify_human_gut',
        'mgnify_marine',
        'mgnify_soil',
        'uhgg'
    ]
    
    print("="*60)
    print("CREATING TEST BATCHES")
    print("="*60)
    print(f"Metadata directory: {metadata_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Genomes per batch: {args.num_genomes}")
    print("")
    
    total_genomes = 0
    
    for db in databases:
        print(f"\n{db}:")
        
        input_file = metadata_dir / f"{db}_parsed.json"
        output_batch_dir = output_dir / db / "tier_1_high_quality"
        output_file = output_batch_dir / "batch_0000.json"
        
        if not input_file.exists():
            print(f"  WARNING: Input file not found: {input_file}")
            continue
        
        try:
            count = create_test_batch(input_file, output_file, args.num_genomes)
            total_genomes += count
        except Exception as e:
            print(f"  ERROR: {e}")
    
    print("")
    print("="*60)
    print(f"COMPLETE: Created batches for {total_genomes} total genomes")
    print("="*60)

if __name__ == '__main__':
    main()