#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

db=$1
query=$2

blastp="blastp"		# Blast+
threads=1

# Prepare file name
dbname="${db%.*}"
queryname="${query%.*}"
blastout="$dbname"-"$queryname".blast

# Do blast if file doesn't exist
if [ -s $blastout ] ; then
	echo "[$(date +'%F %T')]" "File $blastout already exists: no blast done. If you want to redo the comparison, just delete the file and restart the script."
else
	echo "[$(date +'%F %T')]" "Blastp $db vs $query"
	$blastp -db $db -query $query -evalue 1E-5 -outfmt 6 -out $blastout -num_threads $threads
fi

