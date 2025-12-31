#!/bin/bash
################################################################################
# Full Database Test Script
# Tests genome downloads from all 8 databases (10 genomes each)
################################################################################

set -e  # Exit on error

# ============================================================================
# CONFIGURATION - CHANGE THESE PATHS
# ============================================================================

METADATA_DIR="/ibex/user/rioszemm/the-retron-project/huge_databases/second_try/parsed_metadafa"
PROJECT_DIR="/ibex/user/rioszemm/the-retron-project/download_data_"
TEST_OUTPUT_DIR="${PROJECT_DIR}/test_output_full"
TEST_BATCHES_DIR="${PROJECT_DIR}/test_batches_full"
# CONDA_ENV="/ibex/user/rioszemm/conda-environments/retron_tradicional"
NCBI_EMAIL="melissa.rioszertuche@kaust.edu.sa"  # ⚠️ CHANGE THIS!

# Test configuration
GENOMES_PER_DATABASE=10
MAX_WORKERS=2

# ============================================================================
# SETUP
# ============================================================================

echo "========================================"
echo "FULL DATABASE TEST SUITE"
echo "========================================"
echo ""
echo "Configuration:"
echo "  Metadata dir: $METADATA_DIR"
echo "  Test output: $TEST_OUTPUT_DIR"
echo "  Genomes per DB: $GENOMES_PER_DATABASE"
echo "  NCBI email: $NCBI_EMAIL"
echo ""

# Check NCBI email is set
if [ "$NCBI_EMAIL" = "your_email@example.com" ]; then
    echo "❌ ERROR: Please set your NCBI_EMAIL in the script!"
    echo "   Edit this script and change NCBI_EMAIL to your actual email"
    exit 1
fi


# Create directories
mkdir -p "$TEST_OUTPUT_DIR"
mkdir -p "$TEST_BATCHES_DIR"

# Change to project directory
cd "$PROJECT_DIR"

# ============================================================================
# DATABASES TO TEST
# ============================================================================

databases=(
    "gtdb_bacteria"
    "gtdb_archaea"
    "ncbi_bacteria"
    "ncbi_archaea"
    "mgnify_human_gut"
    "mgnify_marine"
    "mgnify_soil"
    "uhgg"
)

# ============================================================================
# STEP 1: CREATE TEST BATCHES
# ============================================================================

echo "========================================"
echo "STEP 1: Creating Test Batches"
echo "========================================"
echo ""

for db in "${databases[@]}"; do
    echo "Creating batch for: $db"
    
    INPUT_FILE="${METADATA_DIR}/${db}_parsed.json"
    OUTPUT_DIR="${TEST_BATCHES_DIR}/${db}"
    
    # Check if input file exists
    if [ ! -f "$INPUT_FILE" ]; then
        echo "  ⚠️  WARNING: Metadata file not found: $INPUT_FILE"
        echo "  Skipping $db"
        echo ""
        continue
    fi
    
    # Create batch
    # python prepare_download_jobs.py \
    #     --input "$INPUT_FILE" \
    #     --output-dir "$OUTPUT_DIR" \
    #     --database "$db" \
    #     --batch-size "$GENOMES_PER_DATABASE" \
    #     --metadata-dir /ibex/user/rioszemm/the-retron-project/huge_databases/second_try \
    #     --filter-quality high


    # Create batch
    python prepare_download_jobs.py \
        --output-dir "$OUTPUT_DIR" \
        --metadata-dir /ibex/user/rioszemm/the-retron-project/huge_databases/second_try/parsed_metadafa \

        # --batch-size "$GENOMES_PER_DATABASE" \
        # --filter-quality high
    
    if [ $? -eq 0 ]; then
        echo "  ✅ Batch created successfully"
    else
        echo "  ❌ Failed to create batch"
    fi
    
    echo ""
done

# ============================================================================
# STEP 2: DOWNLOAD GENOMES
# ============================================================================

echo "========================================"
echo "STEP 2: Downloading Test Genomes"
echo "========================================"
echo ""

download_start_time=$(date +%s)

for db in "${databases[@]}"; do
    echo "========================================"
    echo "Testing: $db"
    echo "========================================"
    
    BATCH_FILE="${TEST_BATCHES_DIR}/${db}/tier_1_high_quality/batch_0000.json"
    
    # Check if batch file exists
    if [ ! -f "$BATCH_FILE" ]; then
        echo "⚠️  Batch file not found: $BATCH_FILE"
        echo "Skipping $db"
        echo ""
        continue
    fi
    
    # Count genomes in batch
    genome_count=$(python3 -c "import json; print(len(json.load(open('$BATCH_FILE'))))" 2>/dev/null)
    echo "Batch contains: $genome_count genomes"
    echo "Starting download..."
    echo ""
    
    # Download
    db_start_time=$(date +%s)
    
    python -m genome_fetcher.fetch_genomes \
        --input "$BATCH_FILE" \
        --output "$TEST_OUTPUT_DIR" \
        --input-type json \
        --required-files genome.fasta \
        --compress \
        --ncbi-email "$NCBI_EMAIL" \
        --parallel \
        --max-workers "$MAX_WORKERS"
    
    download_exit_code=$?
    db_end_time=$(date +%s)
    db_duration=$((db_end_time - db_start_time))
    
    if [ $download_exit_code -eq 0 ]; then
        echo "✅ $db completed successfully in ${db_duration}s"
    else
        echo "❌ $db failed with exit code: $download_exit_code"
    fi
    
    echo ""
done

download_end_time=$(date +%s)
total_duration=$((download_end_time - download_start_time))

# ============================================================================
# STEP 3: ANALYZE RESULTS
# ============================================================================

echo "========================================"
echo "STEP 3: Analyzing Results"
echo "========================================"
echo ""

# Create results summary
RESULTS_FILE="${TEST_OUTPUT_DIR}/test_results_summary.txt"

{
    echo "========================================"
    echo "FULL DATABASE TEST RESULTS"
    echo "========================================"
    echo "Date: $(date)"
    echo "Total duration: ${total_duration}s"
    echo ""
    echo "========================================"
    echo "PER-DATABASE RESULTS"
    echo "========================================"
    echo ""
} > "$RESULTS_FILE"

total_genomes=0
total_successful=0
total_failed=0
total_size=0

for db in "${databases[@]}"; do
    echo "Analyzing: $db"
    
    db_dir="$TEST_OUTPUT_DIR/$db"
    
    if [ -d "$db_dir" ]; then
        # Count downloaded genomes
        genome_count=$(find "$db_dir" -name "genome.fasta.gz" 2>/dev/null | wc -l)
        
        # Get total size
        if command -v du &> /dev/null; then
            db_size=$(du -sh "$db_dir" 2>/dev/null | cut -f1)
            db_size_bytes=$(du -sb "$db_dir" 2>/dev/null | cut -f1)
        else
            db_size="N/A"
            db_size_bytes=0
        fi
        
        # Check manifest
        manifest_file="$db_dir/manifest.json"
        if [ -f "$manifest_file" ]; then
            successful=$(python3 -c "import json; data=json.load(open('$manifest_file')); print(data.get('successful', 0))" 2>/dev/null || echo "0")
            failed=$(python3 -c "import json; data=json.load(open('$manifest_file')); print(data.get('failed', 0))" 2>/dev/null || echo "0")
        else
            successful=$genome_count
            failed=0
        fi
        
        # Update totals
        total_genomes=$((total_genomes + genome_count))
        total_successful=$((total_successful + successful))
        total_failed=$((total_failed + failed))
        total_size=$((total_size + db_size_bytes))
        
        # Print to console
        echo "  Genomes: $genome_count"
        echo "  Successful: $successful"
        echo "  Failed: $failed"
        echo "  Size: $db_size"
        
        # Write to file
        {
            echo "Database: $db"
            echo "  Downloaded: $genome_count genomes"
            echo "  Successful: $successful"
            echo "  Failed: $failed"
            echo "  Total size: $db_size"
            echo ""
        } >> "$RESULTS_FILE"
        
        # Sample file check
        sample_file=$(find "$db_dir" -name "genome.fasta.gz" 2>/dev/null | head -1)
        if [ -n "$sample_file" ]; then
            file_size=$(ls -lh "$sample_file" | awk '{print $5}')
            echo "  Sample file: $(basename $(dirname $sample_file)) ($file_size)"
            
            # Try to validate gzip
            if gunzip -t "$sample_file" 2>/dev/null; then
                echo "  ✅ Sample file is valid gzip"
            else
                echo "  ⚠️  Sample file gzip validation failed"
            fi
        fi
    else
        echo "  ⚠️  No output directory found"
        
        {
            echo "Database: $db"
            echo "  Status: NOT TESTED (no output directory)"
            echo ""
        } >> "$RESULTS_FILE"
    fi
    
    echo ""
done

# ============================================================================
# STEP 4: SUMMARY
# ============================================================================

# Convert total size to human readable
total_size_mb=$((total_size / 1024 / 1024))

{
    echo "========================================"
    echo "OVERALL SUMMARY"
    echo "========================================"
    echo "Total genomes downloaded: $total_genomes"
    echo "Successful: $total_successful"
    echo "Failed: $total_failed"
    echo "Total size: ${total_size_mb} MB"
    echo "Total time: ${total_duration}s"
    
    if [ $total_genomes -gt 0 ]; then
        avg_time=$((total_duration / total_genomes))
        echo "Average time per genome: ${avg_time}s"
    fi
    
    echo ""
    echo "========================================"
    echo "SUCCESS RATE"
    echo "========================================"
    
    if [ $total_genomes -gt 0 ]; then
        success_rate=$((total_successful * 100 / total_genomes))
        echo "Success rate: ${success_rate}%"
        
        if [ $success_rate -ge 90 ]; then
            echo "Status: ✅ EXCELLENT"
        elif [ $success_rate -ge 70 ]; then
            echo "Status: ⚠️  ACCEPTABLE"
        else
            echo "Status: ❌ NEEDS ATTENTION"
        fi
    else
        echo "Status: ❌ NO GENOMES DOWNLOADED"
    fi
    
} | tee -a "$RESULTS_FILE"

echo ""
echo "========================================"
echo "TEST COMPLETE"
echo "========================================"
echo ""
echo "Results saved to: $RESULTS_FILE"
echo "Output directory: $TEST_OUTPUT_DIR"
echo ""
echo "Next steps:"
echo "  1. Review results: cat $RESULTS_FILE"
echo "  2. Check downloaded files: ls -lh $TEST_OUTPUT_DIR/*/"
echo "  3. If successful, proceed with full-scale downloads"
echo ""

# ============================================================================
# STEP 5: DETAILED FILE LISTING
# ============================================================================

echo "Sample of downloaded files:"
find "$TEST_OUTPUT_DIR" -name "genome.fasta.gz" -exec ls -lh {} \; 2>/dev/null | head -20

exit 0