#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

db=$1
query=$2
threads=$3

blastp="blastp"		# Blast+

# Prepare file name
dbname="${db%.*}"
queryname="${query%.*}"
blastout="$dbname"-"$queryname".blast

function echo_log {
	echo "[$(date +'%F %T')] $1" 1>&2
}

# Do blast if file doesn't exist
if [ -s $blastout ] ; then
	echo_log "Reusing $blastout"
else
	$blastp -db $db -query $query -evalue 1E-5 -outfmt 6 -out $blastout -num_threads $threads
	echo_log "Blastp done for $db vs $query"
fi
