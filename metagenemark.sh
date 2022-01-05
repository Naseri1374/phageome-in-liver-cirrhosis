#!/bin/bash
for i in * ; do

cd $i
cp ../MetaGeneMark_v1.mod .

cat GC*_genomic.fna > allVirus

gmhmmp -o ${i}_cds_from_genomic.fna -f L -a -d -A ${i}pro -D ${i}nt -m MetaGeneMark_v1.mod allVirus

rm allVirus
cd ..

done 
