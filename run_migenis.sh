#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

# Stop with any error
set -e

BLAST=""
GENES_DATA=""
DATABASE=""
GENOMES=""

# Default parameters
pairs=2
bre=6
para_id=40

# Optionals
descrip=""
author=""

function usage {
	if [ -n "$1" ]; then
		echo "[ $1 ]"
	fi
	read -d '' help << '_EOF_' || true
	usage: run_migenis.sh -b blast.txt -g genes.genes -d database.sqlite
	
	Files:
	-i <path> : blast output of every genome against each other
	-g <path> : list of all genes extracted from the gbk files with gbk_parser.pl
	-d <path> : database to be created
	-G <path> : list of genomes names
	
	Parameters:
	-p <int>  : number of gaps allowed in blocks
	-b <int>  : number of distant blocks allowed in a break
	-P <int>  : minimum identity of paralogs
	
	Optional:
	-N <str>  : database description
	-A <str>  : database authors
_EOF_
	echo "$help"
	echo "Defaults:"
	echo "gaps=$pair"
	echo "blocks_in_breaks=$bre"
	echo "paralogs_id=$para_id"
	exit
}
export -f usage

####################################################
# Run with all parameters
while getopts "i:g:d:G:p:b:P:N:A:h" option
do
	case $option in
		i)
			BLAST=$OPTARG
			;;
		g)
			GENES_DATA=$OPTARG
			;;
		d)
			DATABASE=$OPTARG
			;;
		G)
			GENOMES=$OPTARG
			;;
		p)
			pairs=$OPTARG
			;;
		b)
			bre=$OPTARG
			;;
		p)
			para_id=$OPTARG
			;;
		N)
			descrip=$OPTARG
			;;
		A)
			author=$OPTARG
			;;
		h)
			usage
			;;
	esac
done
if [ -z "$BLAST" ] ; then usage "Blast file needed (-i)"; fi
if [ -z "$GENES_DATA" ] ; then usage "Genes data needed (-g)"; fi
if [ -z "$DATABASE" ] ; then usage "Database needed (-d)"; fi
if [ -z "$author" ] ; then author=""; fi
if [ -z "$descrip" ] ; then descrip=""; fi
if [ -z "$bre" ] ; then bre=6; fi

ortho="orthos_"$GENES_DATA
paro="paras_"$GENES_DATA
errpath=$DATABASE".err"

function echo_log {
	echo "[$(date +'%F %T')] $1" 1>&2
}

if [ ! -s $DATABASE ]; then
	echo_log "Compute $DATABASE"
	echo_log "Make ortholog pairs"
	ortholog_pairs.pl -i $BLAST -g $GENES_DATA -o $ortho || exit 1
	echo_log "Find paralogs"
	paralog_pairs.pl -i $BLAST -g $GENES_DATA -o $paro -s $para_id || exit 1
	
	# Then create the database
	if [ -z "$GENOMES" ]; then
		par_genome="";
	else
		par_genome="-G $GENOMES"
	fi
	
	init_db.pl -i "$ortho" -p "$paro" -g "$GENES_DATA" $par_genome -d "$DATABASE" -A "$author" -N "$descrip" 2>> "$errpath" || exit 1
	
	echo_log "Find blocks"
	find_blocks.pl -d $DATABASE -t $pairs 2>> $errpath || exit 1
	
	echo_log "Find breaks"
	find_breaks.pl -d $DATABASE -b $bre 2>> $errpath || exit 1
	
	# Finalize
	make_genes_lists.pl -d $DATABASE 2>> $errpath || exit 1
	echo_log "After MAKE GENES LISTS"
	ranking.pl -d $DATABASE -C 2>> $errpath || exit 1
	echo_log "After RANKING"
	order_parts.pl -d $DATABASE -a 2>> $errpath || exit 1
	echo_log "Computing GOC"
	goc.py $DATABASE
	echo_log "Databse $DATABASE created"
else
	echo_log "The database already exists"
fi