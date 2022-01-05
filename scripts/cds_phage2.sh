#!/bin/bash
while read i; do 

mkdir $i 
cd $i

head -n1 ../id.txt > id22.txt

grep -w $(cat id22.txt) ../assembly_summary.txt | cut -f 20 > ddd.txt

awk 'BEGIN{FS=OFS="/";filesuffix="*"} {ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ddd.txt > download_files.sh

source download_files.sh

cp GCF* $i.fna.gz

cd ..
sed -i '1d' id.txt
done < id.txt 
