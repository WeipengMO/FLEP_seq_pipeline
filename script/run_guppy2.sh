#BSUB -J o_run_guppy
#BSUB -n 8
#BSUB -o %J.stdout
#BSUB -e %J.stderr
#BSUB -q q4gpu
#BSUB -gpu "num=4"


time guppy_basecaller \
  -i fast5 \
  -s guppy_out \
  -c dna_r9.4.1_450bps_hac.cfg \
  --recursive \
  --disable_pings \
  --qscore_filtering \
  --device "cuda:all:100%" 

sendWechat.py

