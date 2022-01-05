#!/bin/bash
while read i; do 

cd $i
cat id22.txt > a

ncbi-genome-download --format fasta --taxid a viral --section genbank

cp genbank/viral/G*/G* .
rm a
cd ..

sed -i '1d' id.txt
done < id.txt 
