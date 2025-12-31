#!/usr/bin/env python3
"""
Quick GTDB Downloader - Based on your old working code
Uses the simple filenames without version numbers
"""

import requests
import gzip
from pathlib import Path
import sys

def download_gtdb(output_dir="/ibex/user/rioszemm/the-retron-project/huge_databases/second_try/metadata_files"):
    """Download GTDB using the old working method"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Try multiple base URLs (mirrors)
    base_urls = [
        "https://data.gtdb.ecogenomic.org/releases/latest/",
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/",
        "https://data.gtdb.aau.ecogenomic.org/releases/latest/"
    ]
    
    # Simple filenames (no version numbers - this is what worked before!)
    file_mapping = {
        'bacteria': 'bac120_metadata.tsv.gz',
        'archaea': 'ar53_metadata.tsv.gz'
    }
    
    session = requests.Session()
    session.headers.update({'User-Agent': 'GTDB-Downloader/1.0'})
    
    for organism, filename in file_mapping.items():
        output_file = output_dir / f"gtdb_{organism}_metadata.tsv.gz"
        
        if output_file.exists():
            print(f"✓ File exists: {output_file.name}")
            # Count genomes
            try:
                with gzip.open(output_file, 'rt') as f:
                    count = sum(1 for _ in f) - 1  # Subtract header
                print(f"  {count:,} genomes")
            except:
                pass
            continue
        
        print(f"\n📥 Downloading GTDB {organism}...")
        
        downloaded = False
        for base_url in base_urls:
            url = base_url + filename
            print(f"  Trying: {url}")
            
            try:
                response = session.get(url, stream=True, timeout=60)
                
                if response.status_code == 200:
                    total_size = int(response.headers.get('content-length', 0))
                    downloaded_size = 0
                    
                    with open(output_file, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                downloaded_size += len(chunk)
                                if total_size > 0:
                                    percent = 100 * downloaded_size / total_size
                                    print(f"\r  Progress: {percent:.1f}%", end='', flush=True)
                    
                    print()  # New line after progress
                    
                    # Verify and count
                    size_mb = output_file.stat().st_size / (1024 * 1024)
                    print(f"  ✓ Downloaded: {size_mb:.2f} MB")
                    
                    try:
                        with gzip.open(output_file, 'rt') as f:
                            count = sum(1 for _ in f) - 1
                        print(f"  ✓ {count:,} genomes")
                    except:
                        pass
                    
                    downloaded = True
                    break
                
                else:
                    print(f"  ✗ HTTP {response.status_code}")
            
            except Exception as e:
                print(f"  ✗ Error: {e}")
        
        if not downloaded:
            print(f"  ❌ Failed to download {organism}")
            print(f"  Manual download:")
            print(f"    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/{filename} \\")
            print(f"      -O {output_file}")
            return False
    
    print("\n✅ GTDB download complete!")
    return True


if __name__ == '__main__':
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "/ibex/user/rioszemm/the-retron-project/huge_databases/second_try/metadata_files"
    
    print("=" * 80)
    print("GTDB Quick Downloader")
    print("=" * 80)
    print(f"Output: {output_dir}\n")
    
    success = download_gtdb(output_dir)
    
    if success:
        print("\n🎉 All done! GTDB metadata files ready to use.")
    else:
        print("\n⚠️ Some downloads failed. See manual commands above.")
        sys.exit(1)