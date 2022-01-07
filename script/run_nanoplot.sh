#BSUB -J run_nanoplot
#BSUB -n 36
#BSUB -o %J.stdout
#BSUB -e %J.stderr

NanoPlot \
  --threads 32 \
  --summary sequencing_summary.txt \
  --loglength \
  -o nanoplot
