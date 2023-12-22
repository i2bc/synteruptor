#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
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
	-j <num>  : number of threads
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
if [ -z "$BLOCKS_TOLERANCE" ] ; then BLOCKS_TOLERANCE=2; fi
if [ -z "$DIR" ] ; then usage "Gbk directory needed (-i)"; fi

cd $DIR

# Delete any previous files with the same name
rm ./$NAME* -f
rm ./*.fa* -f

# Prepare genes file
echo "[$(date +'%F %T')] Begin Migenis database creation for $NAME"
shopt -s extglob
list=`ls +(*.gb*|*.genbank|*.dat)`

# Blast
echo "[$(date +'%F %T')] Extract fasta files" >&2
parallel --jobs $JOBS gbk_parser.pl -i {} -f {.}.faa ::: $list

echo "[$(date +'%F %T')] Blast all vs all" >&2
BLAST_FILE=$NAME"_blast.txt"
rm -f $BLAST_FILE
time -p blaster_local.sh -n $JOBS >&2 || exit 1
cat *.blast > $BLAST_FILE

echo "[$(date +'%F %T')] Prepare genes data" >&2
GENES_FILE=$NAME"_genes.txt"
GENOMES_FILE=$NAME"_genomes.txt"
BLASTDB=$NAME".faa"
OPTP=""
if [ -n "$BLOCKS_TOLERANCE" ]; then
	OPTP="-p $BLOCKS_TOLERANCE"
fi
gbk_parser.pl -i "*.gb* *.genbank *.dat" -o $GENES_FILE -g $GENOMES_FILE -f $BLASTDB

# Run the breaks search
echo "[$(date +'%F %T')] Search the breaks and create the database" >&2
DATABASE_FILE=$NAME".sqlite"
rm -f $DATABASE_FILE
echo "[$(date +'%F %T')] time -p run_migenis.sh -i $BLAST_FILE -g $GENES_FILE -d $DATABASE_FILE -G $GENOMES_FILE $OPTP -A '$author' -N '$descrip' >&2";
time -p run_migenis.sh -i $BLAST_FILE -g $GENES_FILE -d $DATABASE_FILE -G $GENOMES_FILE $OPTP -A "$author" -N "$descrip" >&2
if [ $? -eq 0 ]; then
	echo "[$(date +'%F %T')] Database created: $DIR/$DATABASE_FILE" >&2
	echo "[$(date +'%F %T')] Blast database also created: $DIR/$BLASTDB" >&2
else
	echo "[$(date +'%F %T')] Error: run_migenis.sh failed at some point"
	exit $!
fi

