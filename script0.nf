#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.baseDir = "."
params.reads = "$baseDir/data/*_{2,1}.fastq.gz" 
params.outdir = "$baseDir/outputs"
params.refdir = "$baseDir/ref"
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
    
    publishDir params.outdir, mode: 'copy', overwrite: false

    input:
    tuple val(pair_id), file(reads) from read_pairs_ch
    path ref from params.refdir
    output:
    val pair_id into id_ch
    set file("${pair_id}/") into results_ch

    shell
    """
    kb count -x=10XV2 -g="${ref}/transcripts_to_genes.txt"  -i="${ref}/transcriptome.idx" -o="${pair_id}" --tmp="~/kbtemp" --h5ad \
    "${params.baseDir}/${pair_id}_1.fastq.gz" \
    "${params.baseDir}/${pair_id}_2.fastq.gz" \
    """
}