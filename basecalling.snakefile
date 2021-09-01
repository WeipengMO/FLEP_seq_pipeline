'''
@Author       : windz
@Date         : 2020-05-08 09:51:48
@LastEditTime : 2020-05-08 11:26:12

snakemake -j 15 -c 'bsub -J {rulename} -n {threads} -gpu "num=2" -o %J.stdout -e %J.stderr' -s basecalling.snakefile

'''

import os
import glob


DIR_NAME = [os.path.basename(fn).split('.')[0] for fn in glob.glob('fast5/*.sm')]
SAMPLE_NAME = [os.getcwd().split('/')[-1]]

rule all:
    input:
        expand('guppy_out/{dir_name}/sequencing_summary.txt', dir_name=DIR_NAME),
        'fast5/split_files_into_dirs.log',
        #expand('{sample_name}.fastq.gz', sample_name = SAMPLE_NAME)


rule split_fast5_files:
    input:
        'script/split_files_into_dirs.py'
    output:
        'fast5/split_files_into_dirs.log'
    threads: 1
    shell:
        '''
python script/split_files_into_dirs.py --dir_path fast5/ --suffix fast5  --split_num 15
        '''


rule run_guppy:
    input:
        'fast5/split_files_into_dirs.log',  # log文件，只是用来等上一步完成
        'fast5/{dir_name}.sm'  # 空文件，guppy输入是目录
    output:
        'guppy_out/{dir_name}/sequencing_summary.txt'
    threads: 8
    params:
        dir='{dir_name}',
        out='guppy_out/{dir_name}/',
        gpu=2
    shell:
        '''
guppy_basecaller -i fast5/{params.dir} -s {params.out} -c dna_r9.4.1_450bps_hac.cfg --recursive --fast5_out --disable_pings --qscore_filtering --device "cuda:all:100%"
        '''
