#!/bin/bash
# Master pipeline for ARCHS4 Direction 3
# Run this after mouse_gene_v2.5.h5 is fully downloaded
set -e

PROJECT_ROOT="/home/scroll/zzhang/transcriptome-aging"
H5_FILE="$PROJECT_ROOT/data/archs4/mouse_gene_v2.5.h5"
EXPECTED_SIZE=38960132574

echo "========================================"
echo "ARCHS4 Direction 3 Master Pipeline"
echo "========================================"

# Check H5 file
if [ ! -f "$H5_FILE" ]; then
    echo "ERROR: H5 file not found: $H5_FILE"
    exit 1
fi

ACTUAL_SIZE=$(stat -c%s "$H5_FILE")
if [ "$ACTUAL_SIZE" -lt "$EXPECTED_SIZE" ]; then
    echo "ERROR: H5 file incomplete ($ACTUAL_SIZE / $EXPECTED_SIZE bytes)"
    echo "Wait for download to complete, then rerun."
    exit 1
fi

echo "H5 file verified: $ACTUAL_SIZE bytes"

# Step 1: Explore H5 and extract metadata
echo ""
echo "Step 1: Extract metadata"
python3 "$PROJECT_ROOT/src/python/direction3/explore_archs4_h5.py"

# Step 2: Extract counts and run tAge prediction
echo ""
echo "Step 2: Extract counts and predict tAge"
python3 "$PROJECT_ROOT/src/python/direction3/extract_and_predict_archs4.py"

# Step 3: Analysis and figures
echo ""
echo "Step 3: Analysis and figures"
python3 "$PROJECT_ROOT/src/python/direction3/analyze_archs4_predictions.py"

# Step 4: Generate report
echo ""
echo "Step 4: Generate report"
python3 "$PROJECT_ROOT/src/python/direction3/generate_report.py"

echo ""
echo "========================================"
echo "Pipeline complete!"
echo "Results: $PROJECT_ROOT/results/direction3/"
echo "Figures: $PROJECT_ROOT/figures/direction3/"
echo "Report:  $PROJECT_ROOT/docs/REPORT_Direction3_ARCHS4_Mining.md"
echo "========================================"
