#!/bin/bash
git clone https://github.com/juychen/aws-workshop.git ;
cd aws-workshop ;
cat << \EOF > $HOME/.nextflow/config
docker.enabled = true
EOF

aws s3 cp s3://awsscwsbucket/ref/transcripts_to_genes.txt ref/transcripts_to_genes.txt & \
aws s3 cp s3://awsscwsbucket/ref/transcriptome.idx ref/transcriptome.idx & \
aws s3 cp s3://awsscwsbucket/seqs/SRR11181959_2.fastq.gz data/SRR11537951_1.fastq.gz &\
aws s3 cp s3://awsscwsbucket/seqs/SRR11181959_1.fastq.gz data/SRR11537951_2.fastq.gz;

kb count -x=10XV2 -g="ref/transcripts_to_genes.txt"  -i="ref/transcriptome.idx" -o="SRR11537951" --tmp="~/kbtemp" --h5ad \
"data/SRR11537951_1.fastq.gz" \
"data/SRR11537951_2.fastq.gz" ;

aws s3 cp "SRR11537951" s3://${BUCKET_NAME_RESULTS}/outputs/batch/