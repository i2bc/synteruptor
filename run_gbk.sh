#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

# Stop with any error
set -e

DIR=""
NAME="syntebase"
BLOCKS_TOLERANCE=2

# Optionals
descrip=""
author=""
JOBS=0

function usage {
	if [ -n "$1" ]; then
		echo "[ $1 ]"
	fi
	read -d '' help << '_EOF_' || true
	usage: run_migenis.sh -i path/to/gbk/files/ -n dbname
	
	Files:
	-i <path> : path to the gbk files directory
	-n <str>  : name of the database
	
	Optional:
	-p <num>  : blocks tolerance (default: 2)
	-N <str>  : database description
	-A <str>  : database authors
	-j <num>  : number of threads (default: 1)
_EOF_
	echo "$help"
	exit
}
export -f usage

####################################################
# Run with all parameters
while getopts "i:n:p:N:A:j:h" option
do
	case $option in
		i)
			DIR=$OPTARG
			;;
		n)
			NAME=$OPTARG
			;;
		p)
			BLOCKS_TOLERANCE=$OPTARG
			;;
		N)
			descrip=$OPTARG
			;;
		A)
			author=$OPTARG
			;;
		j)
			JOBS=$OPTARG
			;;
		h)
			usage
			;;
	esac
done
if [ -z "$author" ] ; then author=""; fi
if [ -z "$descrip" ] ; then descrip=""; fi
if [ -z "$JOBS" ] ; then JOBS=1; fi
if [ -z "$BLOCKS_TOLERANCE" ] ; then BLOCKS_TOLERANCE=2; fi
if [ -z "$DIR" ] ; then usage "Gbk directory needed (-i)"; fi
DIR=$(realpath $DIR)

function echo_log {
	echo "[$(date +'%F %T')] $1" 1>&2
}

WORK_DIR=$DIR/temp
mkdir -p $WORK_DIR

# Delete any previous files with the same name
cd $WORK_DIR
rm ./*.tmp ./*.faa* -f
cp $DIR/*.* $WORK_DIR/

# Check usable files
n_files=$(ls *.gb* *.dat *.txt *.embl 2> /dev/null | wc -l)

if [ $n_files -lt 2 ]
then
	echo_log "Not enough files to build a Synteruptor database ($n_files). Verify the extension used"
	exit 1
fi

# Remove special characters from file names
rename 's/ /_/g' *.*
rename 's/[^A-Za-z0-9_\-.]+//g' *.*

# Prepare genes file
echo_log "Begin Migenis database creation for $NAME with $n_files files"
shopt -s extglob
list=`ls +(*.gb*|*.dat|*.txt|*.embl)`

# Blast
echo_log "Extract fasta files"
parallel --jobs $JOBS gbk_parser.pl -i {} -f {.}.faa ::: $list

# Check that all files have sequences
empty_fasta="0"
for fasta in $(ls *.faa); do
	nseqs=$(grep '>' $fasta | wc -l)
	if [ "$nseqs" == "0" ]
	then
		echo "Fasta file $fasta has no sequences."
		empty_fasta=1
	fi
done
if [ "$empty_fasta" == "1" ]
then
	echo "Some fasta files had no sequences. Check the input files."
	exit 1
fi

echo_log "Blast all vs all"
BLAST_FILE=$NAME"_blast.txt.tmp"
rm -f $BLAST_FILE
blaster_local.sh -n $JOBS >&2 || exit 1
cat *.blast > $BLAST_FILE

echo_log "Prepare genes data"
GENES_FILE=$NAME"_genes.txt.tmp"
GENOMES_FILE=$NAME"_genomes.txt.tmp"
BLASTDB=$NAME".faa"
OPTP=""
if [ -n "$BLOCKS_TOLERANCE" ]; then
	OPTP="-p $BLOCKS_TOLERANCE"
fi
gbk_parser.pl -i "*.gb* *.dat *.txt *.embl" -o $GENES_FILE -g $GENOMES_FILE -f $BLASTDB

# Run the breaks search
echo_log "Search for breaks and create the database"
DATABASE_FILE=$NAME".sqlite"
rm -f $DATABASE_FILE
run_migenis.sh -i $BLAST_FILE -g $GENES_FILE -d $DATABASE_FILE -G $GENOMES_FILE $OPTP -A "$author" -N "$descrip" >&2
if [ $? -eq 0 ]; then
	cp $WORK_DIR/$DATABASE_FILE $DIR/
	cp $WORK_DIR/$BLAST_DB $DIR/
	echo_log "Database created: $DIR/$DATABASE_FILE"
	echo_log "Blast database also created: $DIR/$BLASTDB"
else
	echo_log "Error: run_migenis.sh failed at some point"
	exit $!
fi
