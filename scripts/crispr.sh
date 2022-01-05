#!/bin/bash
while read i; do

cd $i

ls |awk '/genomic/ && !/rna/ && !/cds/ && !/gff/ && !/gbff/ && !/gtf/ '> cc.txt


while read line ; do 
cd $i
cp ../sel392v2.so .

CRISPRCasFinder.pl -cf CasFinder-2.0.2 General -cas -in $line -out Results -ccc 20000 -ccvRep -keep -html -rcfowce -def S -cpuM 4 -copyCSS

cat NC_* > virus_database.fa


cd Results/CRISPRFinderProperties/NC*/Spacers/
cat *fasta > spacers.fa
cd ../
cp Results/CRISPRFinderProperties/NC*/Spacers/spacers.fa .

blastn -query spacers.fa -subject virus_database.fa -word_size 7 -evalue 1 -reward 1 -penalty -1 -outfmt 6 -gapopen 10 -gapextend 2 -out Blast_crispr.txt

cd ..

done < cc.txt
done < ss.txt
