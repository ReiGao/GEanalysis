#########################################################################
#	File Name: run.sh
#	> Author: QiangGao
#	> Mail: qgao@genetics.ac.cn 
# Created Time: Mon 30 Oct 2017 02:16:22 PM CST
#########################################################################
#!/bin/bash
i=$1
perl /public-supool/home/gaolab/scripts/splitbarcodeCLEANdataPE.pl $i $i.csv ./cleandata/$i/$i\_clean_r1.fq.gz
#perl /public-supool/home/gaolab/scripts/splitbarcodeCLEANdataSE.pl $i $i.csv ./cleandata/$i/$i\_clean_r1.fq.gz
