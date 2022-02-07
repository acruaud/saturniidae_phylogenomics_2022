#!/bin/bash
#SBATCH -p workq
#SBATCH -J shuffle
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

#purpose : random selection of x columns in a phylip file (without repetition) to generate y data sets
#x and y are provided by users

#usage
#create datadir
#copy this script and your phylip file to datadir
#submit job to Slurm :
#sbatch shuffle_sites_phylip.sh file.phy NbOfSitesToPick NbOfDatasetsToGenerate
#output : NbOfDatasetsToGenerate data sets called file.phy_random_$i with $i ranging from 1 to NbOfSitesToPick
#example :
#sbatch shuffle_sites_phylip.sh merge.phy 10000 10
#will create 10 data sets of 10,000 random sites sampled from the merge.phy file and data sets will be named merge.phy_random_1; merge.phy_random_2; ... ; merge.phy_random_10

#record ntaxa from phylip file
ntax=$(head -n 1 $1 |awk '{print $1}')
echo "nb of taxa in your data set=" $ntax

#record nbp from phylip file
nbp=$(head -n 1 $1 |awk '{print $2}')
echo "nb of bp in your alignment=" $nbp

#record input filename
file=$1
echo "your input file is "$file

#record nb of nt to include in final data sets
nsite=$2
echo "you want me to keep "$nsite" random sites in each output data set"

#record nb of data sets to be created
ndataset=$3
echo "you want me to create "$ndataset" data set(s) with "$nsite" random sites"

#extract taxa label from phylip file
awk '{if(NR>1) {print $1}}' $1 > tmp_mytaxa

#extract sequences from phylip file
awk '{if(NR>1) {print $2}}' $1 > tmp_myseq

#insert space between nt
sed 's/\([^ ]\)/\1 /g' tmp_myseq > tmp_myseq2

#cleaning
rm tmp_myseq

#create liste of sites(=column nb = position in the alignment) to retain in each data set
#one file is created per data set
for i in `seq 1 $ndataset`
do
seq 1 $nbp |shuf |head -n $nsite |tr '\n' ' ' |sed 's/\<[0-9]*\>/$&/g' > tmp_sitetokeep_$i
done

#if subsets are composed of more than 10,000 sites the next cmd will return a "list of argument too long" error
#I need to split the  tmp_sitetokeep_* into list of 10,000 sites or less
for i in tmp_sitetokeep_*
do
sed -i s/" "/"\n"/g $i
done

for i in tmp_sitetokeep_*
do
split -l 10000 $i $i"_split"
done

for i in tmp_sitetokeep_*_split*
do
cat $i |tr '\n' ' ' > $i".tmp"
done

#create and execute cmd to extract sites
for i in tmp_sitetokeep_*_split*.tmp
do
awk '{print "awk \047{print "$0"}\047 tmp_myseq2 > "FILENAME"_seq"}' $i |bash
done

#create cmd to paste data subsets of 10,000 bp
for i in `seq 1 $ndataset`
do
ls "tmp_sitetokeep_"$i"_split"*".tmp_seq" |tr '\n' ' ' > tmp_mylist_$i
done

for i in `seq 1 $ndataset`
do
awk '{print "paste -d\042\\0\042 "$0" > tmp_sitetokeep_"'${i}'"_merge"}' tmp_mylist_$i |bash
done

#paste with taxa name
for i in tmp_sitetokeep_*_merge
do
paste -d " " tmp_mytaxa $i > $i"_def"
done

#add ntax and nbp to file
ls tmp_sitetokeep_*_merge_def > tmp_myfiles

awk -v ntax=$ntax -v bp=$nsite -v fname=$file  '{print "awk \047BEGIN{print \042"ntax" "bp"\042} {print $0}\047 "$0" > "fname"_"$0"_def2"}' tmp_myfiles |bash

#rename final files
rename 'tmp_sitetokeep' 'random' *_tmp_sitetokeep_*_merge_def_def2
rename '_merge_def_def2' '' *_merge_def_def2

#convert unique space to 4 spaces to fit with MCMCtree requirements (at least 2 spaces between label and sequence)
sed -i '2,$s/ /    /g' *_random_*

#cleaning
rm tmp_*

