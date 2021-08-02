'''
Author       : windz
Date         : 2020-05-08 09:51:48
LastEditTime : 2020-12-07 09:51:37
Description  : A pipline for calling polya tail length

    You can run like this:
    snakemake -j 8 -p -c 'bsub -J {rulename} -n {threads} -o log/%J.stdout -e log/%J.stderr'
'''


# create log dir
if not os.path.exists('log'):
    os.mkdir('log')


configfile: 'config.yml'
SAMPLE_NAME=config['sample_name']
genome_data=config['genome']


rule all:
    input:
        expand('results/{sample_name}.read.info.txt', sample_name=SAMPLE_NAME),
        expand('results/{sample_name}.read.splicing_kinetics.txt', sample_name=SAMPLE_NAME),
        expand('results/{sample_name}.read.rna.ir.stat', sample_name=SAMPLE_NAME),


rule fastq_to_fasta:
    input:
        'guppy_out/pass/',
    output:
        'basecalled_data/{sample_name}.fasta'
    threads: 1
    shell:
        '''
python script/fastqdir2fasta.py --indir {input} --out {output}
        '''


rule mapping_to_genome:
    input:
        'basecalled_data/{sample_name}.fasta'
    output:
        bam='aligned_data/{sample_name}.sorted.bam'
    params:
        genome=genome_data
    threads: 32
    shell:
        '''
minimap2 -t {threads} -ax splice --secondary=no -G 12000 {params.genome} {input} | samtools sort -@ {threads} -o {output.bam} -
samtools index -@ {threads} {output.bam}
        '''


rule find_3linker:
    input:
        bam='aligned_data/{sample_name}.sorted.bam',
        fasta='basecalled_data/{sample_name}.fasta'
    output:
        'aligned_data/{sample_name}.adapter.result.txt'
    threads: 36
    shell:
        '''
python script/adapterFinder.py --inbam {input.bam} --inseq {input.fasta} --out {output} --threads {threads} --mode 1
        '''


rule polyacaller:
    input:
        adapter_result='aligned_data/{sample_name}.adapter.result.txt',
        sequencing_summary='sequencing_summary.txt',
        fast5_dir='guppy_out/workspace',
    output:
        'aligned_data/{sample_name}.polyA_tail.result.txt'
    threads: 36
    shell:
        '''
python script/PolyACaller.py --inadapter {input.adapter_result} --summary {input.sequencing_summary}  --fast5dir {input.fast5_dir} --out {output} --threads {threads}
        '''


rule identify_read_info:
    input:
        'aligned_data/{sample_name}.sorted.bam'
    output:
        'aligned_data/{sample_name}.read_info.result.txt'
    params:
        'genome_data/exon_intron_pos.repr.bed'
    threads: 1
    shell:
        '''
python script/extract_read_info.py --inbed {params} --inbam {input} --out {output}
        '''


rule merge_info:
    input:
        read_info='aligned_data/{sample_name}.read_info.result.txt',
        adapter='aligned_data/{sample_name}.adapter.result.txt',
        polya='aligned_data/{sample_name}.polyA_tail.result.txt'
    output:
        read_info_result='results/{sample_name}.read.info.txt'
    threads: 1
    params:
        'Nanopore'
    shell:
        '''
export PATH=/scem/work/mowp/miniconda3/envs/R/bin/:$PATH
Rscript script/merge_read_info.R --type {params} --inreadinfo {input.read_info} --inadapter {input.adapter} --inpolya {input.polya} --out {output.read_info_result}
        '''


rule splicing_kinetics:
    input:
        read_info='results/{sample_name}.read.info.txt'
    output:
        splicing_data='results/{sample_name}.read.intron.pos.splicing.txt',
        splicing_kinetics='results/{sample_name}.read.splicing_kinetics.txt',
        figure='results/{sample_name}.read.splicing_kinetics.pdf'
    threads: 1
    params:
        inbed='genome_data/exon_intron_pos.repr.bed',
        select_intron='genome_data/select_introns.txt'
    shell:
        '''
python script/prepare_data_for_splice_kinetics.py --inreadinfo {input.read_info} --inbed {params.inbed} --out {output.splicing_data}
export PATH=/scem/work/mowp/miniconda3/envs/R/bin/:$PATH
Rscript script/plot_intron_splicing_kinetics.R --inrelpos {output.splicing_data} --inreadinfo {input.read_info} --inintron {params.select_intron} --out {output.splicing_kinetics} --pdf {output.figure}
        '''


# Calculate intron retention ratio of polyadenylated transcripts
rule intron_retention_ratio:
    input:
        splicing_data='results/{sample_name}.read.intron.pos.splicing.txt',
        read_info='results/{sample_name}.read.info.txt',
    output:
        rna_ir='results/{sample_name}.read.rna.ir.stat',
        intron_ir='results/{sample_name}.read.intron.ir.stat',
    threads: 1
    shell:
        '''
export PATH=/scem/work/mowp/miniconda3/envs/R/bin/:$PATH
Rscript script/cal_polya_transcript_ir.R --inrelpos {input.splicing_data} --inreadinfo {input.read_info} --outrna {output.rna_ir} --outintron {output.intron_ir}
        '''

