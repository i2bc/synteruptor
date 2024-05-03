# Synteruptor scripts

## Description

This is a set of scripts that create a Synteruptor database from a set of genome files.

## Installation

We recommend runnning the Docker image as it is the most straightforward way to run this program.

### Docker image

See the README file in the `docker/` folder for steps to build an image.

### Installation requirements

If you want to install it in your own machine, here's the list of requirements.

#### Recommended Ubuntu packages

Assuming Ubuntu 22.04:

git
parallel
rename
bioperl
libstatistics-basic-perl
ncbi-blast+
sqlite3


#### Details

The programs and libraries required to run the scripts are as follow:
- Perl
- Bioperl

	Hint for local installs if you don't have it in a package: use cpanminus
```
	curl -L http://cpanmin.us | perl - App::cpanminus
	curl -L http://cpanmin.us | perl - Bio::Perl
```
- Statistics::Basic (libstatistics-basic-perl)
- GNU parallel (parallel)
- Blast+ (ncbi-blast+)
- Python 3
- Python 3 libs (either from a package or with manager like pip):

	sqlite3


## Environment

Add the path to the `src/` folder to your `PATH`

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

```
run_gbk.sh -i /path/to/subfolder -n db_name
```
You can run this script with multiple threads by adding the parameter `-j`. E.g. to run on 4 cores:

```
run_gbk.sh -i /path/to/subfolder -n db_name -j 4
```

- It will then create a database named db_name.sqlite in the subfolder, as well as a Blast DB db_name.sqlite.faa
- Place the DB in the `db` folder of the Web Synteruptor to explore its data. Make sure the web server has permission to read the `db` folder and the database file.

## 3 genomes example

As an example, download the gbff files for 3 genomes from NCBI:
- Streptococcus anginosus C1051 (GCF_000463465.1)
- Streptococcus anginosus C238 (GCF_000463505.1)
- Streptococcus anginosus subsp. whileyi MAS624 (GCA_000478925.1)
You can do this through [this link](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1353243,862971,862970).

Unzip the file and move the .gbff files in a subfolder, make sure to rename each file to something meaningful.

Then run the script on this subfolder as above (assuming 4 cores):
```
run_gbk.sh -i /path/to/subfolder -n db_name -j 4
```

## Database content

### Tables

The database contains the following tables at initiation:
- `genes`
- `genomes`
- `genome_parts`(separate DNA sequences)
The data in those tables come from the input files.

When running run_migenis.sh the following tables are populated:
- `orthos`(orthologs and paralogs genes pairs, from BLAST BRH)
- `info`(metadata)
- `pairs`(orthologs and paralog pairs)
- `blocks`(synteny blocks)
- `breaks`(synteny breaks)
- `breaks_genes`(all genes in a break, per species)
- `breaks_ranking`
- `breaks_graph` (to represent similar breaks among a group of species)
- `goc`

Additional tables and views created to ease the queries:
- `orthos_all`(joined pairs with genes species1 and genes species2)
- `blocks_all`(joined block, orthos pairs at start of block, end of block + genes at start and end in both species)
- `breaks_all`(joined break, blocks left and right, orthos pairs left and right, genes at start and end un both species)

The final breaks data are stored in breaks (breaks_all for more data) and breaks_ranking (contains the various attributes used for ranking).

You can use sqlite3 `.schema` command to see the details of each table.

