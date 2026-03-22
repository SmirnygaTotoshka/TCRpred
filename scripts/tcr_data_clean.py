#!/usr/bin/env python

import os
import time
import uuid
import click
import warnings
import mhcgnomes
import traceback
import numpy as np
import pandas as pd

PROTEIN_CHECK="^[ACDEFGHIKLMNPQRSTVWY]+$"
HOST_SPECIES = ['human', 'mouse']
CDR3_CHAINS = ['alpha', 'beta']

warnings.filterwarnings("ignore")
input_option = click.option("-i","--input","input", required = True, type=str, help = "Input file")
output_option = click.option("-o","--output","output", required = True, type=str, help = "Output directory")

def is_valid_mhc(allele):
    try:
        _ = mhcgnomes.parse(allele)
    except (mhcgnomes.errors.ParseError, TypeError):
        return False
    else:
        return True

def calculate_receptor_id(data: pd.DataFrame, column: str) -> pd.DataFrame:
    cleaned_data = data.copy(deep=True)
    unique_id = cleaned_data[column].unique()
    replacement = {k:str(uuid.uuid4()) for k in unique_id}
    cleaned_data[column] = cleaned_data[column].map(replacement)
    return cleaned_data

@click.group()
def main():
    pass


@main.command()
@input_option
@output_option
@click.option("-d","--database","database", required = True, type=click.Choice(['IEDB', 'CEDAR'], case_sensitive=True))
def mysql(input, output,database):
    raw_data = pd.read_csv(input, sep = ";", header = 0)
    species_columns_to_clean = ["species_name_receptor","species_name_chain","mhc_source","host_organism"]
    print("Unificate species")
    for col in species_columns_to_clean:
        raw_data.loc[raw_data[col].str.contains("Homo sapiens"),col] = "human"
        raw_data.loc[raw_data[col].str.contains("Mus musculus"),col] = "mouse"

    table_schema = {
        "id": "id",
        "chain_type": "chain",
        "cdr3_seq": "sequence",
        "linear_peptide_seq": "epitope",
        "chain_i_name": "mhc_alpha",
        "chain_ii_name": "mhc_beta",
        "class": "mhc_class",
        "host_organism": "host_species",
        "epitope_source": "epitope_species",
        "source_antigen_accession": "epitope_source",
        "v_gene": "V",
        "d_gene": "D",
        "j_gene": "J",
        "database": "database"
    }
    print("Check mhc names")
    raw_data['valid_i_name'] = raw_data['chain_i_name'].apply(lambda x: is_valid_mhc(x))
    raw_data['valid_ii_name'] = raw_data['chain_ii_name'].apply(lambda x: is_valid_mhc(x))

    print("Apply filters")
    clean_data = raw_data.\
            query("reference_id.notna()").\
            query("receptor_type_x == 'alphabeta'").\
            query("chain_type.isin(@CDR3_CHAINS)").\
            query("species_name_receptor == species_name_chain == host_organism == mhc_source and species_name_receptor.isin(@HOST_SPECIES)").\
            query("cdr3_seq.notna()").\
            query("cdr3_seq.str.contains(@PROTEIN_CHECK)").\
            query("e_region_domain_flag == 'Exact Epitope'").\
            query("linear_peptide_seq.notna()").\
            query("linear_peptide_seq.str.contains(@PROTEIN_CHECK)").\
            query("linear_peptide_modification.isna()").\
            query("disc_region.isna()").\
            query("ant_type == 'Epitope'").\
            query("valid_i_name and valid_ii_name").\
            query("as_char_value.str.contains('Positive')")
    
    print("Get id")
    unique_id = clean_data["distinct_receptor_id"].unique()
    replacement = {k:str(uuid.uuid4()) for k in unique_id}
    clean_data["id"] = clean_data["distinct_receptor_id"].map(replacement)
    clean_data["database"] = database
    clean_data_selected = clean_data.filter(items = list(table_schema.keys()), axis = 1).rename(columns = table_schema)
    
    print("Reshape")
    pivoted_part = clean_data_selected.filter(["id", "chain", "sequence","V","D","J"],axis=1).sort_values(by='V', na_position='last')
    annotation_part = clean_data_selected.filter(["id", "epitope", "mhc_alpha","mhc_beta","mhc_class","host_species","epitope_species","epitope_source","database"],axis=1)
    pivoted_part_pivot = pivoted_part.pivot_table(index = "id",columns='chain', values=['sequence', 'V', 'D', 'J'], aggfunc = 'first').reset_index()
    wide = pd.merge(pivoted_part_pivot, annotation_part, on="id")

    print("Save")
    wide.to_csv(os.path.join(output, f"{database}_clean.csv"), sep = ";", index = False)

    
@main.command()
@input_option
@output_option
def vdjdb(input, output):
    raw_data = pd.read_csv(input, sep = ";", header = 0)

@main.command()
@input_option
@output_option
def mcpas(input, output):
    raw_data = pd.read_csv(input, sep = ";", header = 0)

@main.command()
@input_option
@output_option
def pird(input, output):
    raw_data = pd.read_csv(input, sep = ";", header = 0)
            
@main.command()
@input_option
@output_option
def mira(input, output):
    raw_data = pd.read_csv(input, sep = ";", header = 0)

if __name__ == "__main__":
    main()