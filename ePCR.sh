#!/bin/bash

echo "print out files containing primer ID, Sequences, and product size ..."

for i in `ls *.csv`; 
do 
cat $i|awk -F, '{print$1"\t"$3"\t"$8}' >$i.temp.txt; 
done

echo "removing header, and print Primer_F,Primer_R and combine them into one line ..."

for i in `ls *.temp.txt`; 
do 
sed -i '1d' $i; 
cat $i|sed -n 'p;n' >$i.1; 
cat $i|sed -n 'n;p' >$i.2; 
paste -d "\t" $i.1 $i.2|perl -p -i -e "s/\r//g"|sort|uniq >$i.temp2.txt;
rm $i.1 $i.2; 
done  

echo "generating new ePCR file ..."

for i in `ls *.temp2.txt`; 
do 
cat header >$i.temp3;
sort $i|uniq |awk -F '\t' '{print$1$4"\t"$2"\t"$5"\t"$6}' >>$i.temp3; 
mv -f $i.temp3 `echo $i.temp3|sed 's/csv.temp.txt.temp2.txt.temp3/ePCR/g'`; 
done
rm *.temp.txt *.temp2.txt

echo "OK! Let's do some ePCR ..."
echo "please input your Database:"
read Database

for i in `ls *.ePCR`; 
do 
../e-PCR-2.3.12/e-PCR -n 5 -g 1 -o $i.out -t 3 $i $Database; 
done

mv *.ePCR ePCR_out
mv *.csv Soltis_et_al_2010.*.fasta Primer_file
