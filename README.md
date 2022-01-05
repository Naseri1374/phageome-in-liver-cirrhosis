# phageome-in-liver-cirrhosis
This repository provide details for the analyses performed this study:
Investigation and characterization of human gut phageome in liver cirrhosis of defined etiologies




Here, special attention has been paid to detecting the role of gut phageome in liver diseases; hence, future studies need to delve into the characteristics of the gut phageome in various liver disorders. Moreover, access to a reference gene catalog is inevitable to deeply investigate microbial environments such as the human gut. In this regard, a better understanding of the bacteriophage content of the human gut depends on access to a comprehensive gene catalog of the gut phageome. In the present study, we adopted three bioinformatic strategies to identify the large scaffolds of phage origin to investigate the structural and functional composition of gut phageome in several etiologies of liver diseases. This is the first study to correlate the gut phageome with liver cirrhosis to provide valuable insights into the design of novel phage-related markers in monitoring liver disease progression and its treatment.

## Tools and Dependencies

SRAtoolkit\
FastViromeExplorer\
megahit\
edger\
phyloseq

## Input dataset
In this study, three previously published metagenomic datasets were employed for the analyses under the accession number of PRJEB6337, PRJNA373901, and PRJEB18041.

## Phage catalog construction
Here, we proposed a method for developing phage gene-catalog. The pipeline used three strategies to construct the gut phage catalog to facilitate the taxonomic and functional analysis of the gut phageome in liver cirrhosis.

![image](https://user-images.githubusercontent.com/39089097/147871774-87f8e140-cdc7-43db-b668-e1f1e1a15ca3.png)


## Analysis pipeline

###### Download Data:

$ ./prefetch $(<SraAccList.txt)

$ ./fastq-dump --split-3 --gzip $(<SraAccList.txt)

$ ./fastqc *.fq

$ bowtie2 -p 8 -x GRCh38_noalt_as -1 SAMPLE_R1.fastq.gz -2 SAMPLE_R2.fastq.gz --un-conc SAMPLE_host_removed > SAMPLE_mapped_and_unmapped.sam

$ for i in *_1.fastq; do java -jar Trimmomatic-0.36/trimmomatic-0.36.jar PE $i ${i/_1/_2} -baseout ${i/_sm/.fastq ILLUMINACLIP:Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:7 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:80; done


###### FastViromeExplorer:

$ for i in *_1.fastq; do java -cp bin FastViromeExplorer -1 $i -2 ${i/_1/_2} -i ncbi-virus-kallisto-index-k31.idx -o output${i/_*.fastq/}; done

$ ncbi-acc-download --format fasta NC_016158 NC_011357

$ csplit -z NC_011043,NC_029014.fa '/>/' {*}


###### NCBI Genome and Accession Download Scripts:

$ ncbi-genome-download --format fasta,assembly-report --genus T2D_Meta_ACC_genome.txt bacteria

$ ncbi-genome-download -l complete --parallel 3 --retries 10 -v -F protein-fasta --genus T2D_Meta_ACC_genome.txt bacteria


###### NCBI ftp genome downzload:
$ grep -E 'GCF_000312305|GCF_000312645' assembly_summary_refseq.txt | cut -f 20 > ftp_folder.txt
$ awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_fna_files.sh
$ source download_fna_files.sh


###### megahit:
$ for i in ERR*_1.fastq.gz; do ./megahit -1 $i -2 ${i/_1/_2} -m 0.9  -t 12 -o assembly_${i/_*.fastq.gz/} ; done


###### convert Acssesion number to taxid:
$ esearch -db nucleotide -query NC_031915.1,NC_031918.1|esummary|xtract -pattern TaxId -element TaxId >nnn.txt

###### get full linage of taxid
$ taxonkit lineage --data-dir /media/biocool/aaaaaaaaa/taxdump taxid_control_cirr.txt| tee lineage_control_cirr.txt
