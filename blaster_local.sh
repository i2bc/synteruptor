#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

# Check if parallel and blast+ are installed
if ! command -v parallel >/dev/null; then
	echo "GNU parallel is required. Abort."
	exit 1;
fi
if ! command -v makeblastdb >/dev/null; then
	echo "Blast++ is required. Abort."
	exit 1;
fi

DATA_PATH=.
DATA_FILES=*.fa*
JOBS=""

usage() {
	if [ -n "$1" ]; then
		echo "[ $1 ]"
	fi
	read -d '' help << '_EOF_' || true
	usage: local_blaster.sh
	
	Files:
	-i <path> : directory path where the fasta files are
	-f <path> : file descriptor for the fasta files (default: *.fa*)
	-n <int>  : max number of threads (default: 1)
_EOF_
	echo "$help"
	exit
}
export -f usage

####################################################
# Run with all parameters
while getopts "i:f:n:h" option
do
	case $option in
		i)
			DATA_PATH=$OPTARG
			;;
		f)
			DATA_FILES=$OPTARG
			;;
		n)
			JOBS=$OPTARG
			;;
		h)
			usage
			;;
	esac
done
if [ -z "$JOBS" ] ; then JOBS=1; fi

####################################################
cd $DATA_PATH;

# Loop through files
files=`ls $DATA_FILES`

# Prepare all blastdb
(parallel -j $JOBS makeblastdb -dbtype prot -in {} ::: $files) 1> /dev/null 2> /dev/null

# Run in parallel, like 2 nested loops all files vs all files
parallel -j 1 blaster_function.sh {1} {2} $JOBS ::: $files ::: $files || exit 1
