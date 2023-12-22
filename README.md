# Synteruptor scripts

## Description

This is a set of scripts that create a Synteruptor database from a set of genome files.

## Installation requirements (with ubuntu packages names):

- Perl
- Bioperl
	Hint for local installs: use cpanminus
	curl -L http://cpanmin.us | perl - App::cpanminus
	curl -L http://cpanmin.us | perl - Bio::Perl
- Statistics::Basic (libstatistics-basic-perl)
- GNU parallel (parallel)
- Blast+ (ncbi-blast+)
- Python 3
- Python 3 libs:
	sqlite3

## Environment

Add the path to the src/ folder to your PATH

## Input

Synteruptor requires each genome assembly to be in a file each with both DNA sequence and annotations.

File formats supported are:
- Genbank (.gb*)
- EMBL (.dat, .embl, .txt)

### Note

Make sure you have at least 2 genome files to compare.

## Usage

- Move to a working directory
- Place genomic files in a subfolder
    NB: intermediate files will be created in that subfolder
- Run the script run_gbk.sh

$ run_gbk.sh -i /path/to/subfolder -n db_name

- It will then create a database named db_name.sqlite in the subfolder, as well as a Blast DB db_name.sqlite.faa
- Place the DB in the db folder of the Web Synteruptor to explore its data.
