#!/bin/sh --login

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

#PBS -N blaster
#PBS -l ncpus=4

# Task number needed
: ${PBS_ARRAY_INDEX:?"Need to set PBS_ARRAY_INDEX non-empty"}

# DATAPATH needed
: ${DATA_PATH:?"Need to set DATA_PATH non-empty"}
cd $DATA_PATH

# DATA_FILES needed
: ${DATA_FILES:?"Need to set DATA_FILES (e.g. *.faa) non-empty"}

# Choose the couple to use
i=0
db=''
query=''
files=`ls $DATA_FILES`
for f1 in $files; do
	for f2 in $files; do
		i=`expr $i + 1`
		if [ $i -eq $PBS_ARRAY_INDEX ]; then
			echo "$i\tdb=$f1\tquery=$f2"
			db=$f1
			query=$f2
			break 2
		fi
	done
done

if [ -z $db ] || [ -z $query ]; then
	echo "db and query are not defined (wrong task number maybe? Compare $PBS_ARRAY_INDEX and $i)"
	echo "Path: $DATA_PATH"
	echo "Files: $DATA_FILES"
	echo "found files: $files"
	exit
fi

# blastp="blastp-2.2.29+"		# On ebio
blastp="blastp"		# On i2bc
dbname="${db%.*}"
queryname="${query%.*}"
blastout="$dbname"-"$queryname".blast

# Don't do it again if the file already exists!
if [ -s $blastout ] ; then
	echo "File $blastout already exists: no blast done. If you want to redo the comparison, just delete the file and resubmit the job (array id: $PBS_ARRAY_INDEX)."
else
	$blastp -db $db -query $query -evalue 1E-5 -outfmt 6 -out $blastout -num_threads 4
fi
