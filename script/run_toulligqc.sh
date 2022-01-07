#BSUB -J run_toulligqc
#BSUB -n 2
#BSUB -o %J.stdout
#BSUB -e %J.stderr

toulligqc \
    --barcoding \
    --telemetry-source          guppy_out/sequencing_telemetry.js \
    --sequencing-summary-source guppy_out/sequencing_summary.txt \
    --barcodes                  barcode01,barcode02 \
    --output-directory          ToulligQC