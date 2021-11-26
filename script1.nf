#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.reads = "s3://awsscwsbucket/seqs/SRR11537951/*_{2,1}.fastq.gz" 
params.ref = "s3://awsscwsbucket/ref/"
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
    input:
    tuple val(SRR_id), file(reads) from read_pairs_ch
    path ref from params.ref
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
process Analysis {
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}", mode: "copy"
    input:
    file result from results_ch
    output:
    file('alignment_results/') into results2_ch
    
    shell
    """
    ls
    cd alignment_results
    ls
    python ${params.codebase}/analysis.py "counts_unfiltered/adata.h5ad"
    cd ../
    """
}