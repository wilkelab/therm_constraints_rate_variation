#!/bin/bash

#Create Working Directory
#$WDIR=$SCRATCH/$JOB_NAME-$JOB_ID-$1
WDIR=$SCRATCH/rosetta_ddG_calculations/$2

mkdir -p $WDIR
if [ ! -d $WDIR ] 
then
	echo $WDIR not created
	exit
fi

#Move into working directory
cd $WDIR

#Copy pdb
cp $WORK/MyProjects/ddG_files/renumbered_pdbs/$3 .

#Copy flag file for minimization
cp $WORK/MyProjects/ddG_files/minimize_flag_files/$4 .

#Copy flag file for ddG monomer calculations
cp $WORK/MyProjects/ddG_files/ddg_flag_files/$1/$5 .

#Copy resfile for ddG monomer calculations
cp $WORK/MyProjects/ddG_files/resfiles/$1/$2_resfile .

#Copy pdb_list
cp $WORK/MyProjects/ddG_files/pdb_lists/$1_pdb_list.txt .

#Copy constraint maker
cp $WORK/MyProjects/ddG_files/convert_to_cst_file.sh .

#Copy weights file
cp $WORK/MyProjects/ddG_files/sp2_paper_talaris2013_scaled.wts .

date
hostname

#Put science related commands here

#Minimize the structure
$WORK/rosetta_2014wk05_bundle/main/source/bin/minimize_with_cst.linuxiccrelease @$4 > $2_mincst.log

#Generate the constraints file
./convert_to_cst_file.sh $2_mincst.log > $1_input.cst

#Run ddG Monomer
$WORK/rosetta_2014wk05_bundle/main/source/bin/ddg_monomer.linuxiccrelease @$5 > $2_out.txt

date

#Copy Results Back to the Home Directory

RDIR=$WORK/MyProjects/rosetta_ddg_results/$2
#RDIR=$HOME/MyProjects/rosetta_ddg_results/$2
mkdir -p $RDIR
cp $WDIR/* $RDIR

rm $RDIR/$2/mut_*.pdb
rm $RDIR/$2/repacked_wt*.pdb

#Cleanup
#rm -rf $WDIR
