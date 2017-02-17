#!/bin/bash

#this script is coded by Miao. It's used to summarize the ePCR results to see primer efficiency
# need to prepare a gene list before hand, for which gene the primers are designed

#mkdir ePCR_file

for i in `cat gene.list`; do
sed -i '1d' ${i}.ePCR;
#cp $i.ePCR ./ePCR_file/ && sed -i '1d' $i.ePCR;
awk '{$0=v"_"$0}1' v="$i" ${i}.ePCR|sort|uniq >>summarypart1.txt;

done

for i in `cat gene.list`; do

awk -F '\t' '{print$2"\t"$1}' ${i}.ePCR.out|sort|uniq|sed 's/^/'$i'_/g' >>ePCR_count.tmp;
done

cut -f1 ePCR_count.tmp| sort|uniq -c|awk -F ' ' '{print$2"\t"$1}' >>ePCR_Count.txt

join --nocheck-order summarypart1.txt ePCR_Count.txt >Primer_ePCR_summary.txt.tmp

cat header2 Primer_ePCR_summary.txt.tmp >Primer_ePCR_summary.txt

#sed -i 's/ /\t/g' Primer_ePCR_summary.txt

#rm *.tmp
