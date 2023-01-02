#!/bin/sh

#### -------------------------------- Script header -------------------------------- ####
# Date:         09.12.2022                                                              #
# Author:       Siri Naerland Skodvin                                                   #
# Filename:     MFchr.sh                                                                #
# Description:  Sequential runs of PLINK, to recode genetic data, and R, to analyze     #
#               and plot results.                                                       #
#### --------------------------------------------------------------------------------####

# Input chromosome - CHANGE HERE
chrchr=05       # chromosome number in two characters
chrnum=5        # chromosome number with minimum characters
startkbp=0      # start kbp of chromosome
stopkbp=190000  # stop kbp of chromosome

# Input chunks - CAN USE DEFAULT
bykbp=200       # number of kbp to include in each gen recode file
gengrp=6        # number of gen recode files per chunk

# Input pruning - CAN USE DEFAULT
window=50
step=5
r2=0.95

# Preprocess input
cpath=/PATH_TO_CODE_AND_RESULTS/
dpath=/PATH_TO_TEMPORARY_STORAGE_AND_MAIN_R_SCRIPTS/
seq=`expr $stopkbp - $startkbp`
mod=`expr $seq % $bykbp`
if [ $mod -ne 0 ]
then
	stopkbp=`expr $stopkbp - $mod + $bykbp`
	seq=`expr $stopkbp - $startkbp`
fi
genno=`expr $seq / $bykbp`
genchunks=`expr \( $seq + \( \( $bykbp \* $gengrp \) - 1 \) \) / \( $bykbp \* $gengrp \)`
chunkno=1

# Initiate log file
now=`date`
cat <<EOF > $dpath/log_chr$chrchr.txt
INPUT
Chromosome: $chrnum
Start-kbp: $startkbp
Stop-kbp: $stopkbp
By-kbp: $bykbp
Grouped by: $gengrp
Sequences: $genchunks
Prune window-step-r2: $window - $step - $r2
ANALYSIS STARTED
$now
ANALYZED SEQUENCES
EOF

# Main analysis - Plink and R in sequence
i=1
while [ $i -le $genno ]
do
	if [ `expr $i + $gengrp - 1` -le $genno ]
	then
		last=`expr $i + $gengrp - 1`
	else
		last=$genno
	fi
	
	for file in `ls $dpath`
	do
		filepre=`expr "$file" | cut -c1-3`
		if [ $filepre = "gen" ]
		then
			rm $dpath/$file
		fi
		if [ $filepre = "ld." ]
		then
			rm $dpath/$file
		fi
	done
	
	j=1
	while [ $i -le $last ]
	do
		echo ---- PLINK RECODE PART $i of $genno ----
		filename=$dpath/gen$j
		fromkb=`expr $startkbp + \( $bykbp \* \( $i - 1 \) \)`
		if [ $j -eq 1 ]
		then
			fromkbini=$fromkb
		fi
		tokb=`expr $fromkb + $bykbp`
		
		/opt/plink --fam GENETICS.fam --bim GENETICS.bim --bed GENETICS.bed --maf 0.01 --chr $chrnum --from-kb $fromkb --to-kb $tokb --indep-pairwise $window $step $r2 --out $dpath/ld
		
		/opt/plink --fam GENETICS.fam --bim GENETICS.bim --bed GENETICS.bed --maf 0.01 --chr $chrnum --extract $dpath/ld.prune.in --recode --out $filename

		for file in `ls $dpath`
		do
			filepre=`expr "$file" | cut -c1-3`
			if [ $filepre = "ld." ]
			then
				rm $dpath/$file
			fi
		done
		
		i=`expr $i + 1`
		j=`expr $j + 1`
	done
	echo ---- R INTERACTION ANALYSIS SEQUENCE $chunkno of $genchunks ----
		fi
	done
	
	k=0
	for file in `ls $dpath`
	do
		filesuf=`expr "$file" | cut -c6-8`
		if [ $filesuf = "ped" ]
		then
			k=`expr $k + 1`
	
	if [ $k -eq 0 ]
	then
		newline="seq $chunkno: ${fromkbini} to ${tokb} kbp => no variants"
	else
		Rscript --vanilla $dpath/MFinter.R
		newline="seq $chunkno: ${fromkbini} to ${tokb} kbp"
	fi
	
	echo $newline >>$dpath/log_chr$chrchr.txt
	chunkno=`expr $chunkno + 1`
done

Rscript --vanilla $dpath/MFmplot.R

lastline="ANALYSIS FINISHED"
now=`date`
vartotline="NUMBER OF VARIANTS"
vartot=`cat $cpath/numvar$chrchr.txt`
echo $lastline >>$dpath/log_chr$chrchr.txt 
echo $now >>$dpath/log_chr$chrchr.txt
echo $vartotline >>$dpath/log_chr$chrchr.txt
echo $vartot >>$dpath/log_chr$chrchr.txt
cp $dpath/log_chr$chrchr.txt $cpath
rm $cpath/numvar$chrchr.txt

for file in `ls $cpath`
do
	resfile=`expr "$file" | cut -c1-4`
	if [ $resfile = "2022" ]
	then
		rm $cpath/$file
	fi
done

echo ---- ANALYSIS FINISHED ----
