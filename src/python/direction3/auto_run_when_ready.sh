#!/bin/bash
# Watcher that auto-runs the ARCHS4 pipeline when download completes
H5="/home/scroll/zzhang/transcriptome-aging/data/archs4/mouse_gene_v2.5.h5"
EXPECTED=38960132574
PIPELINE="/home/scroll/zzhang/transcriptome-aging/src/python/direction3/run_pipeline.sh"
LOG="/home/scroll/zzhang/transcriptome-aging/logs/direction3/auto_run.log"

while true; do
    if [ -f "$H5" ]; then
        SIZE=$(stat -c%s "$H5")
        PCT=$(awk "BEGIN {printf \"%.2f\", $SIZE/$EXPECTED*100}")
        echo "$(date): $SIZE bytes ($PCT%)" >> "$LOG"
        if [ "$SIZE" -ge "$EXPECTED" ]; then
            echo "$(date): Download complete. Starting pipeline." >> "$LOG"
            bash "$PIPELINE" >> "$LOG" 2>&1
            echo "$(date): Pipeline finished with exit code $?" >> "$LOG"
            break
        fi
    fi
    sleep 300  # check every 5 minutes
done
