#!/usr/bin/env python3

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

import sqlite3, argparse, math, os, sys

parser = argparse.ArgumentParser(
    description="""This script computes GOC for pairwise genomes in the database"""
)
parser.add_argument("database", type=str, help="path to the sqlite3 database")
parser.add_argument("-v", dest="verbose", action="store_true", help="verbose mode")
args = parser.parse_args()

db_file = args.database

# Load database
if not os.path.isfile(db_file):
    sys.exit("exit")

conn = sqlite3.connect(db_file)

cur = conn.cursor()

cur.execute("""DROP TABLE IF EXISTS goc""")
cur.execute("""CREATE TABLE goc(sp1 TEXT, sp2 TEXT, pos INTEGER, score REAL)""")

# Get all species in the database
cur.execute("SELECT DISTINCT sp FROM genes;")
list_species = [x[0] for x in cur.fetchall()]

c = 0
window_proportion = 3
total_comput = (len(list_species) ** 2) - len(list_species)
for ref in list_species:
    for tar in list_species:
        if ref == tar:
            continue
        c += 1
        if args.verbose:
            print(f"\r\t{c}/{total_comput} GOC computation", end="")

        # Get all CDS for the reference
        cur.execute(
            """
            SELECT g.pid FROM genes g JOIN genome_parts gp ON g.gpart = gp.gpart
            WHERE feat = 'CDS' and g.sp = ?
            ORDER BY gp.min, loc_start ASC;
            """,
            (ref,),
        )
        list_cds_ref = [x[0] for x in cur.fetchall()]

        # Get all genes for the reference
        cur.execute(
            """
            SELECT g.pid FROM genes g JOIN genome_parts gp ON g.gpart = gp.gpart
            WHERE g.sp = ?
            ORDER BY gp.min, loc_start ASC;
            """,
            (ref,),
        )
        list_gene_ref = [x[0] for x in cur.fetchall()]

        # Get all CDS for the target
        cur.execute(
            """
            SELECT g.pid FROM genes g JOIN genome_parts gp ON g.gpart = gp.gpart
            WHERE feat = 'CDS' and g.sp = ?
            ORDER BY gp.min, loc_start ASC
            """,
            (tar,),
        )
        list_cds_tar = [x[0] for x in cur.fetchall()]

        # Get all ortholog pairs
        sql_prep_cds_list = "'" + ("', '").join(list_cds_ref) + "'"
        sql_handler = f"""
            SELECT g1.pid gene_id1, g2.pid gene_id2
            FROM orthos o JOIN genes g1 ON o.pid1 = g1.pid JOIN genes g2 ON o.pid2 = g2.pid
            WHERE o.pid1 IN ({sql_prep_cds_list}) AND g2.sp = '{tar}'
            """
        cur.execute(sql_handler)
        ort = {}
        raw_sql_ort = cur.fetchall()
        for i in raw_sql_ort:
            ort[i[0]] = i[-1]

        # GOC computation
        goc_cds = []
        goc_loc = []
        limit = len(list_cds_ref)
        window_length = math.ceil(limit / 100.0) * window_proportion
        loc_start_window = math.ceil(window_length / 2.0)
        start_window = 0
        list_tar_ort = [ort[x] if x in ort else "NA" for x in list_cds_ref]
        while start_window < (limit - window_length):
            window = list_tar_ort[start_window : (start_window + window_length)]
            limit_computation = len(window) - 1
            index_tar_max = len(list_cds_tar) - 1
            new_synt_region = True
            number_cds_in_synt_region = 0
            for index_ref, cds_tar in enumerate(window):
                if cds_tar != "NA":
                    if index_ref < limit_computation:
                        index_tar = list_cds_tar.index(window[index_ref])
                        if (
                            index_tar != index_tar_max
                            and list_cds_tar[index_tar + 1] == window[index_ref + 1]
                        ):
                            number_cds_in_synt_region += 1
                            if new_synt_region:
                                number_cds_in_synt_region += 1
                                new_synt_region = False
                            continue
                        elif index_ref > 0 and index_tar != 0:
                            if list_cds_tar[index_tar - 1] == window[index_ref + 1]:
                                number_cds_in_synt_region += 1
                                if new_synt_region:
                                    number_cds_in_synt_region += 1
                                    new_synt_region = False
                            else:
                                new_synt_region = True
                                continue
                else:
                    new_synt_region = True
            goc_cds = number_cds_in_synt_region / len(window)
            # goc_loc = loc_start_window
            goc_loc = list_gene_ref.index(list_cds_ref[loc_start_window])
            start_window += 1
            loc_start_window += 1

            # Insert a row of data
            handler = f"INSERT INTO goc VALUES ('{ref}','{tar}', '{goc_loc}','{goc_cds}')"
            cur.execute(handler)
        # Save the changes
        conn.commit()

if args.verbose:
    print("\n")

cur.close()
conn.close()
