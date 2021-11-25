#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = "s3://awsscwsbucket/seqs/*_{2,1}.fastq.gz"
params.annotation = "s3://awsscwsbucket/ref/"
params.codebase = "~"
params.baseDir = "."
log.info """\
         SCVH - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.annotation}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .view()
    .set { read_pairs_ch }

/*
 * 1. Mapping 
 */
process Map {
    cpus 12
    memory '40 GB'
    publishDir "${params.outdir}/${pair_id}", mode: "copy"

    input:
    tuple val(pair_id), file(reads) from read_pairs_ch
    path ref from params.annotation
    output:
    file("${pair_id}") into results_ch
    
    shell
    """
    kb count -x=10XV2 -g="${ref}/transcripts_to_genes.txt"  -i="${ref}/transcriptome.idx" -o="${pair_id}" --tmp="./kbtemp" -m=32  -t=12 --h5ad \
    "${params.baseDir}/${pair_id}_1.fastq.gz" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" 
    """
}