#!/usr/bin/env python3

import argparse
import polars as pl
import gzip
import io


def main(args):
    # Read both input files
    df1 = pl.read_csv(args.bac, separator="\t", infer_schema_length=1000)
    df2 = pl.read_csv(args.arc, separator="\t", infer_schema_length=1000)

    # Combine into one DataFrame
    df = pl.concat([df1, df2])

    # Remove the "RS_" or "GB_" prefix from the accession column
    df = df.with_columns(
        pl.col("accession")
        .str.replace(r"^(RS_|GB_)", "")
        .alias("accession")
    )

    # Check necessary columns exist
    if "gtdb_taxonomy" not in df.columns or "accession" not in df.columns:
        raise ValueError("Input files must contain 'accession' and 'gtdb_taxonomy' columns.")
    if "ncbi_organism_name" not in df.columns:
        raise ValueError("Input files must contain 'ncbi_organism_name' column.")

    ### Output 1: taxonomy split into domain, phylum, etc.
    taxonomy_split = (
    df.select(["accession", "gtdb_taxonomy"])
    .with_columns(
        pl.col("gtdb_taxonomy")
        .str.split(";")
        .alias("taxonomy_split")
    )
    .with_columns(
        pl.col("taxonomy_split").list.get(0).alias("domain"),
        pl.col("taxonomy_split").list.get(1).alias("phylum"),
        pl.col("taxonomy_split").list.get(2).alias("class"),
        pl.col("taxonomy_split").list.get(3).alias("order"),
        pl.col("taxonomy_split").list.get(4).alias("family"),
        pl.col("taxonomy_split").list.get(5).alias("genus"),
        pl.col("taxonomy_split").list.get(6).alias("species"),
    )
    .select([
        pl.col("accession").alias("ident"),
        "domain", "phylum", "class", "order", "family", "genus", "species"
    ])
)

    csv_buffer = io.StringIO()
    taxonomy_split.write_csv(csv_buffer, separator=",")

    with gzip.open(args.output_taxonomy, 'wb') as f:
        f.write(csv_buffer.getvalue().encode('utf-8'))

    ### Output 2: species representative taxonomy
    
    # Join the representative flag with taxonomy_split
    taxonomy_with_repinfo = taxonomy_split.join(
        df.select(["accession", "gtdb_representative"]),
        left_on="ident",
        right_on="accession",
        how="left"
    )

    rep_taxonomy = taxonomy_with_repinfo.filter(
        pl.col("gtdb_representative") == "t"
    ).drop("gtdb_representative")

    csv_buffer_rep = io.StringIO()
    rep_taxonomy.write_csv(csv_buffer_rep, separator=",")

    with gzip.open(args.reps_taxonomy, 'wb') as f:
        f.write(csv_buffer_rep.getvalue().encode('utf-8'))

    ### Output 3: accession + ncbi_organism_name combined into a name
    name_file = (
        df.select(["accession", "ncbi_organism_name"])
        .with_columns(
            (pl.col("accession") + " " + pl.col("ncbi_organism_name")).alias("name")
        )
        .select([
            pl.col("accession"),
            "name"
        ])
    )

    name_file.write_csv(args.output_gbsketch, separator=",")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process two tsv.gz files and create taxonomy and name CSVs.")
    parser.add_argument("--bac", required=True, help="bacteria metadata .tsv.gz file")
    parser.add_argument("--arc", required=True, help="archaea metadata .tsv.gz file")
    parser.add_argument("--output-taxonomy", required=True, help="Output compressed CSV (.csv.gz) with taxonomy split")
    parser.add_argument("--reps-taxonomy", required=True, help="Output compressed CSV (.csv.gz) for species representative taxonomy")
    parser.add_argument("--output-gbsketch", required=True, help="Output gbsketch CSV (.csv) with accession and name")

    args = parser.parse_args()
    main(args)

