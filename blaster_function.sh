#!/bin/bash

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.


# https://www.meziantou.net/retry-a-bash-command.htm
# Define the retry function
retry() {
  local retries="$1"
  local command="$2"
  local options="$-" # Get the current "set" options

  # Disable set -e
  if [[ $options == *e* ]]; then
    set +e
  fi

  # Run the command, and save the exit code
  $command
  local exit_code=$?

  # restore initial options
  if [[ $options == *e* ]]; then
    set -e
  fi

  # If the exit code is non-zero (i.e. command failed), and we have not
  # reached the maximum number of retries, run the command again
  if [[ $exit_code -ne 0 && $retries -gt 0 ]]; then
    retry $(($retries - 1)) "$command"
  else
    # Return the exit code from the command
    return $exit_code
  fi
}


db=$1
query=$2
threads=$3
threads=${threads:-"1"}

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
	retry 3 "$blastp -db $db -query $query -evalue 1E-5 -outfmt 6 -out $blastout -num_threads $threads"
	echo_log "Blastp done for $db vs $query"
fi

