#!/bin/bash
# GTDB Manual Download Script
# This uses the simple filenames that worked in your old script

cd /ibex/user/rioszemm/the-retron-project/huge_databases/second_try/metadata_files

echo "================================"
echo "Downloading GTDB Metadata Files"
echo "================================"
echo ""

echo "Method 1: Using Australian mirror (recommended)"
echo "------------------------------------------------"

# Bacteria
echo "Downloading bacteria metadata..."
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz \
  -O gtdb_bacteria_metadata.tsv.gz

# Archaea
echo "Downloading archaea metadata..."
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz \
  -O gtdb_archaea_metadata.tsv.gz

echo ""
echo "If the above didn't work, try Method 2:"
echo "---------------------------------------"
echo ""
echo "# Method 2: Using ecogenomic.org mirror"
echo "wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz -O gtdb_bacteria_metadata.tsv.gz"
echo "wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz -O gtdb_archaea_metadata.tsv.gz"
echo ""
echo "# Method 3: Using Danish mirror (if both above fail)"
echo "wget https://data.gtdb.aau.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz -O gtdb_bacteria_metadata.tsv.gz"
echo "wget https://data.gtdb.aau.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz -O gtdb_archaea_metadata.tsv.gz"
echo ""

# Verify downloads
echo "Verifying downloads..."
if [ -f "gtdb_bacteria_metadata.tsv.gz" ]; then
    SIZE=$(du -h gtdb_bacteria_metadata.tsv.gz | cut -f1)
    echo "✓ Bacteria metadata: $SIZE"
    
    # Count genomes
    COUNT=$(zcat gtdb_bacteria_metadata.tsv.gz | wc -l)
    COUNT=$((COUNT - 1))  # Subtract header
    echo "  Genomes: $COUNT"
else
    echo "✗ Bacteria metadata download failed"
fi

if [ -f "gtdb_archaea_metadata.tsv.gz" ]; then
    SIZE=$(du -h gtdb_archaea_metadata.tsv.gz | cut -f1)
    echo "✓ Archaea metadata: $SIZE"
    
    # Count genomes
    COUNT=$(zcat gtdb_archaea_metadata.tsv.gz | wc -l)
    COUNT=$((COUNT - 1))  # Subtract header
    echo "  Genomes: $COUNT"
else
    echo "✗ Archaea metadata download failed"
fi

echo ""
echo "================================"
echo "Download complete!"
echo "================================"