#I predict variant effects and interface sites using pre-computed data
import os, os.path, sqlite3, pandas as pd, numpy as np, functools, scipy as sp, scipy.stats, seaborn as sns, matplotlib.pyplot as plt, matplotlib
from af2genomics import *
import sqlite3
import re

##### Fetch data #####

# Functions
def query_missense_all(variants):
    with sqlite3.connect('../../../../data/biol_databases/AlphaMissense/alphamissense_hs.sqlite') as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + str(tuple(variants)), con=db)
    return df_

#Set up 
prefix = "ubnVariantsVEP"

#Get data
AAchangeDat = pd.read_csv("..//output/data/VEP/ubnVariants.csv")

#Query VEPs from my alphaMissense data (more complete)
queryOutput = query_missense_all(AAchangeDat["variant_id"].tolist())
queryOutput.to_csv(f"..//output/data/VEP/{prefix}_all.csv")

#Query  VEPs from base data
queryOutput = query_missense(AAchangeDat["variant_id"].tolist())
queryOutput.to_csv(f"..//output/data/VEP/{prefix}.csv")

#Query  VEP with interfaces
queryOutput = merge_missense(AAchangeDat, "variant_id")
queryOutput.to_csv(f"..//output/data/VEP/{prefix}_fMerge.csv")


##### Merge data #####
def _prep_common_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Add uniprot_site, uniprot, mutAA (last char of variant_id), normalize slashes."""
    df = df.copy()
    # uniprot_site: drop trailing single capital letter, then replace '/' with '_'
    df["uniprot_site"] = (
        df["variant_id"]
        .str.replace(r"[A-Z]$", "", regex=True)
        .str.replace("/", "_", regex=False)
    )
    # uniprot: everything up to first underscore
    df["uniprot"] = df["uniprot_site"].str.replace(r"_.+$", "", regex=True)
    # mutAA: last character of variant_id
    df["mutAA"] = df["variant_id"].str[-1]
    return df

def _wide_by_mutAA(df: pd.DataFrame, var_cols: list[str]) -> pd.DataFrame:
    """
    For each mutAA, create a block of columns with prefix 'mut{AA}_',
    indexed by uniprot_site; then column-bind all blocks.
    """
    pieces = []
    for aa, block in df.groupby("mutAA", sort=False):
        part = (
            block[["uniprot_site"] + var_cols]
            .drop_duplicates(subset=["uniprot_site"])
            .set_index("uniprot_site")
            .rename(columns=lambda c: f"mut{aa}_{c}")
        )
        pieces.append(part)
    if pieces:
        out = pd.concat(pieces, axis=1)
        out = out.loc[:, ~out.columns.duplicated()]  # guard against dup columns
    else:
        out = pd.DataFrame(index=pd.Index([], name="uniprot_site"))
    # restore key as a column for merging
    out = out.reset_index()
    return out

def _process_vep_main(path_csv: str) -> pd.DataFrame:
    """
    Replicates the first VEP block:
      - invarCols kept once per uniprot_site
      - varCols split per mutAA into wide columns
      - returns one row per uniprot_site with both pieces combined
    """
    invar_cols = [
        "plddt","sasa","pocketscore","pocketrank","interface","interface_strict",
        "freq","xref","significance","association","clintype",
    ]
    var_cols = ["ESM1b_LLR","ESM1b_is_pathogenic","am_pathogenicity","am_class","eve","pred_ddg"]

    df = pd.read_csv(path_csv)
    if df.columns[0].lower() in {"", "unnamed: 0", "index"}:
        df = df.iloc[:, 1:]  # drop first dummy column, like R code

    df = _prep_common_cols(df)

    invar = (
        df[["uniprot_site"] + invar_cols]
        .drop_duplicates(subset=["uniprot_site"])
        .set_index("uniprot_site")
    )

    varwide = _wide_by_mutAA(df[["uniprot_site", "mutAA"] + var_cols], var_cols).set_index("uniprot_site")

    combined = pd.concat([varwide, invar], axis=1).reset_index()
    return combined  # has 'uniprot_site' column

def _process_vep_interface(path_csv: str) -> pd.DataFrame:
    """Replicates the interface subset + distinct."""
    keep = ["interface_pdockq", "interface_label", "interaction_id", "uniprot_site"]
    df = pd.read_csv(path_csv)
    if df.columns[0].lower() in {"", "unnamed: 0", "index"}:
        df = df.iloc[:, 1:]
    df = _prep_common_cols(df)
    return df[keep].drop_duplicates()

def _process_vep_all(path_csv: str) -> pd.DataFrame:
    """
    Replicates the 'all' VEP block:
      - build wide per-AA for ['am_pathogenicity','am_class']
      - then rename *_am_class -> *_am_class_all and *_am_pathogenicity -> *_am_pathogenicity_all
    """
    var_cols = ["am_pathogenicity", "am_class"]
    df = pd.read_csv(path_csv)
    if df.columns[0].lower() in {"", "unnamed: 0", "index"}:
        df = df.iloc[:, 1:]
    df = _prep_common_cols(df)

    wide = _wide_by_mutAA(df[["uniprot_site", "mutAA"] + var_cols], var_cols)
    wide = wide.rename(columns=lambda c: (
        c.replace("am_class", "am_class_all").replace("am_pathogenicity", "am_pathogenicity_all")
        if c.startswith("mut") else c
    ))
    return wide  # has 'uniprot_site'

def merge_all_vep(
    vep_main_csv: str,
    vep_interface_csv: str,
    vep_all_csv: str
) -> pd.DataFrame:
    """
    Process three VEP-like files and combine them with a FULL (outer) merge
    on 'uniprot_site' so we keep the union of keys across sources.
    Then merge once into InigoFeatProc (left join keeps original rows).
    """
    main_proc = _process_vep_main(vep_main_csv)
    int_proc  = _process_vep_interface(vep_interface_csv)
    all_proc  = _process_vep_all(vep_all_csv)

    # FULL/OUTER merge across VEP sources -> comprehensive union of uniprot_site
    merged_vep = (
        main_proc
        .merge(int_proc, on="uniprot_site", how="outer")
        .merge(all_proc, on="uniprot_site", how="outer")
    )
    return merged_vep


merged_vep = merge_all_vep(
    f"..//output/data/VEP/{prefix}.csv",
    f"..//output/data/VEP/{prefix}_fMerge.csv",
    f"..//output/data/VEP/{prefix}_all.csv",
)
merged_vep.to_csv(f"..//output/data/VEP/{prefix}_merged.csv")