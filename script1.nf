#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.reads = "s3://awsscwsbucket/seqs/*_{2,1}.fastq.gz" 
params.refdir = "s3://awsscwsbucket/ref/"
params.codebase = "~"

log.info """\
        - N F   P I P E L I N E -
         ===================================
         references   : ${params.refdir}
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
    input:
    tuple val(SRR_id), file(reads) from read_pairs_ch
    path ref from params.refdir
    output:
    val SRR_id into id_ch
    file("alignment_results/") into results_ch

    shell
    """
    kb count -x=10XV2 -g="${ref}/transcripts_to_genes.txt"  -i="${ref}/transcriptome.idx" -o="alignment_results" --tmp="~/kbtemp" --h5ad \
    "${params.baseDir}/${SRR_id}_1.fastq.gz" \
    "${params.baseDir}/${SRR_id}_2.fastq.gz" \
    """
}
/*
 * 2. Analysis
 */
process Filter {
    publishDir "${params.outdir}", mode: "copy"
    input:
    file result from results_ch
    output:
    file('alignment_results/') into results2_ch
    
    shell
    """
    cd alignment_results
    Rscript ${params.codebase}/analysis.py "adata.h5ad"
    """
}