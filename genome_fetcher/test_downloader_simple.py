#!/usr/bin/env python3
"""
Simple direct test - run this from the genome_fetcher directory.

Usage:
    cd /path/to/genome_fetcher
    python test_downloaders_simple.py
"""

import json
import sys
from pathlib import Path

# Setup paths
script_dir = Path(__file__).resolve().parent
project_root = script_dir.parent if script_dir.name == 'genome_fetcher' else script_dir
sys.path.insert(0, str(project_root))

# Imports
from genome_fetcher.downloaders import get_downloader
from genome_fetcher.utils import DownloadLogger
from genome_fetcher.config import setup_logging
import logging

setup_logging(level=logging.INFO)
logger = logging.getLogger(__name__)


# Configuration
METADATA_DIR = "/ibex/user/rioszemm/the-retron-project/huge_databases/second_try/parsed_metadafa"
OUTPUT_DIR = "./test_output"
SAMPLE_SIZE = 1  # Test just 1 genome per database

DATABASE_FILES = {
    'gtdb_bacteria': 'gtdb_bacteria_parsed.json',
    'gtdb_archaea': 'gtdb_archaea_parsed.json',
    'mgnify_human_gut': 'mgnify_human_gut_parsed.json',
    'mgnify_marine': 'mgnify_marine_parsed.json',
    'mgnify_soil': 'mgnify_soil_parsed.json',
    'uhgg': 'uhgg_parsed.json',
    'ncbi_bacteria': 'ncbi_bacteria_parsed.json',
    'ncbi_archaea': 'ncbi_archaea_parsed.json',
}


def main():
    metadata_dir = Path(METADATA_DIR)
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    print("=" * 80)
    print("DOWNLOADER TEST - Testing all databases")
    print("=" * 80)
    print(f"Metadata: {metadata_dir}")
    print(f"Output:   {output_dir}")
    print(f"Sample:   {SAMPLE_SIZE} genome per database")
    print("=" * 80)
    print()
    
    for db_name, json_file in DATABASE_FILES.items():
        print(f"\n{'='*80}")
        print(f"Testing: {db_name}")
        print('='*80)
        
        # Load genome
        json_path = metadata_dir / json_file
        if not json_path.exists():
            print(f"❌ File not found: {json_file}")
            results[db_name] = 'FILE_NOT_FOUND'
            continue
        
        print(f"Loading: {json_path.name}")
        
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Get first genome
        if isinstance(data, dict):
            genome_id = list(data.keys())[0]
            genome_record = data[genome_id].copy()
            genome_record['genome_id'] = genome_id
        else:
            genome_record = data[0]
            genome_id = genome_record.get('genome_id', genome_record.get('accession', 'unknown'))
        
        # Add database field if not present
        if 'database' not in genome_record:
            genome_record['database'] = db_name
        
        print(f"Test genome: {genome_id}")
        
        # Get downloader
        downloader = get_downloader(db_name)
        if not downloader:
            print(f"❌ No downloader available")
            results[db_name] = 'NO_DOWNLOADER'
            continue
        
        print(f"Downloader: {downloader.__class__.__name__}")
        
        # Set NCBI email if needed
        if db_name.startswith('NCBI') and hasattr(downloader, 'email'):
            downloader.email = "test@example.com"
        
        # Setup output
        db_output_dir = output_dir / db_name
        db_output_dir.mkdir(parents=True, exist_ok=True)
        
        genome_dir = db_output_dir / genome_id
        genome_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = genome_dir / "download.log"
        terminal_log = DownloadLogger(log_file)
        
        # Download
        print("Downloading...")
        try:
            success = downloader.download_genome(
                genome_record=genome_record,
                output_dir=db_output_dir,
                required_files=['genome.fasta'],
                compress=True,
                terminal_log=terminal_log
            )
            
            if success:
                # Check if file exists
                fasta_gz = genome_dir / "genome.fasta.gz"
                fasta = genome_dir / "genome.fasta"
                
                if fasta_gz.exists():
                    size_mb = fasta_gz.stat().st_size / 1024 / 1024
                    print(f"✅ SUCCESS - Downloaded {size_mb:.2f} MB")
                    results[db_name] = 'PASSED'
                elif fasta.exists():
                    size_mb = fasta.stat().st_size / 1024 / 1024
                    print(f"✅ SUCCESS - Downloaded {size_mb:.2f} MB")
                    results[db_name] = 'PASSED'
                else:
                    print(f"⚠️  Download succeeded but file not found")
                    results[db_name] = 'FILE_MISSING'
            else:
                print(f"❌ FAILED - Check log: {log_file}")
                results[db_name] = 'FAILED'
        
        except Exception as e:
            print(f"❌ ERROR: {e}")
            results[db_name] = f'ERROR: {str(e)[:50]}'
            import traceback
            traceback.print_exc()
    
    # Summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for r in results.values() if r == 'PASSED')
    total = len(results)
    
    for db_name, status in results.items():
        symbol = "✅" if status == "PASSED" else "❌"
        print(f"{symbol} {db_name:25} {status}")
    
    print("=" * 80)
    print(f"Passed: {passed}/{total}")
    print("=" * 80)
    
    if passed < total:
        print("\n⚠️  Some tests failed. Check the logs in:", output_dir)
    else:
        print("\n🎉 All tests passed!")
    
    return passed == total


if __name__ == '__main__':
    try:
        success = main()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\nTest interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
