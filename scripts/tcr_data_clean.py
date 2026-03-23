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

def fix_mhc_name(allele_str, chain = 'alpha'):
    try:
        # Предварительные проверки и подготовления
        if pd.isna(allele_str):
            return pd.NA     
        allele_first = allele_str.split(" ")[0]
        if allele_first == "B2M":
            return "B2M"
        # Получение конкретной аллели
        parsed_allele = mhcgnomes.parse(allele_first)
        if isinstance(parsed_allele, mhcgnomes.pair.Pair):
            if chain == 'alpha':
                allele = parsed_allele.alpha
            elif chain == 'beta':
                allele = parsed_allele.beta
            else:
                raise ValueError('Unknown chain')
        elif isinstance(parsed_allele, mhcgnomes.allele.Allele):
            allele = parsed_allele
        else:
            raise ValueError('Not allele')
        # Проверка
        if allele.gene.species.name == "Homo sapiens":
            if len(allele.allele_fields) < 2:
                raise ValueError('Too many allele fields. Need at least 2.')
            elif len(allele.allele_fields) == 2:
                return allele.to_string()
            else:
                return allele.restrict_allele_fields(2, drop_annotations=True, drop_mutations=True).to_string()
        elif allele.gene.species.name == "Mus musculus":
            return allele.to_string()
        else:
            raise ValueError('Wrong species. Only human and mouse are allowed.')
    except (mhcgnomes.errors.ParseError, TypeError, ValueError):
        return pd.NA



@click.group()
def main():
    pass


@main.command()
@input_option
@output_option
@click.option("-d","--database","database", required = True, type=click.Choice(['IEDB', 'CEDAR'], case_sensitive=True))
def mysql(input, output,database):

    def calculate_receptor_id(data: pd.DataFrame, column: str) -> pd.DataFrame:
        cleaned_data = data.copy(deep=True)
        unique_id = cleaned_data[column].unique()
        replacement = {k:str(uuid.uuid4()) for k in unique_id}
        cleaned_data["id"] = cleaned_data[column].map(replacement)
        return cleaned_data
    
    raw_data = pd.read_csv(input, sep = ";", header = 0)
    print(raw_data.shape)
    species_columns_to_clean = ["species_name_receptor","species_name_chain","mhc_source","host_organism"]
    print("Unificate species")
    for col in species_columns_to_clean:
        raw_data.loc[raw_data[col].str.contains("Homo sapiens"),col] = "human"
        raw_data.loc[raw_data[col].str.contains("Mus musculus"),col] = "mouse"
    print("Check mhc names")

    raw_data.loc[raw_data['chain_ii_name'].str.contains("Beta-2-microglobulin"),'chain_ii_name'] = "B2M"
    raw_data['valid_i_name'] = raw_data['chain_i_name'].apply(lambda x: fix_mhc_name(x, 'alpha'))
    raw_data['valid_ii_name'] = raw_data['chain_ii_name'].apply(lambda x: fix_mhc_name(x, 'beta'))

    table_schema = {
            "id": "id",
            "chain_type": "chain",
            "cdr3_seq": "sequence",
            "linear_peptide_seq": "epitope",
            "valid_i_name": "mhc_alpha",
            "valid_ii_name": "mhc_beta",
            "class": "mhc_class",
            "host_organism": "host_species",
            "epitope_source": "epitope_species",
            "source_antigen_accession": "epitope_source",
            "v_gene": "V",
            "d_gene": "D",
            "j_gene": "J",
            "database": "database"
        }


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
            query("valid_i_name.notna() and valid_ii_name.notna()").\
            query("as_char_value.str.contains('Positive')")        
    print(clean_data.shape)
    
    print("Get id")
    clean_data_with_id = calculate_receptor_id(clean_data,"curated_receptor_id")
    clean_data_with_id["database"] = database
    clean_data_selected = clean_data_with_id.filter(items = list(table_schema.keys()), axis = 1).rename(columns = table_schema)
        
    print("Reshape")
    pivoted_part = clean_data_selected.filter(["id", "chain", "sequence","V","D","J"],axis=1).sort_values(by='V', na_position='last')
    annotation_part = clean_data_selected.filter(["id", "epitope", "mhc_alpha","mhc_beta","mhc_class","host_species","epitope_species","epitope_source","database"],axis=1)
    pivoted_part_pivot = pivoted_part.pivot_table(index = "id",columns='chain', values=['sequence', 'V', 'D', 'J'], aggfunc = 'first').reset_index()
    pivoted_part_pivot.columns = pivoted_part_pivot.columns.to_flat_index().str.join('_')
    pivoted_part_pivot.rename(columns = {"id_":"id", "sequence_alpha":"cdr3_alpha", "sequence_beta":'cdr3_beta'},inplace=True)
    wide = pd.merge(pivoted_part_pivot, annotation_part, on="id")
    print("Save")
    wide.to_csv(os.path.join(output, f"{database}_clean.csv"), sep = ";", index = False)

    
@main.command()
@input_option
@output_option
def vdjdb(input, output):

    def calculate_receptor_id(tbl: pd.DataFrame, column_id: str):
        data = tbl.copy(deep=True)
        unique_id = tbl[column_id].unique()
        # receptor_id == 0 means that receptor isn`t paired
        # they need own id, not just replace
        index = np.argwhere(unique_id == 0)
        paired_id = np.delete(unique_id, index)
        replacement = {k:str(uuid.uuid4()) for k in paired_id}
        data["id"] = tbl[column_id].map(replacement)
        
        for i in data.index:
            if pd.isna(data.loc[i, "id"]):
                data.loc[i,"id"] = str(uuid.uuid4())

        return data
    
    raw_data = pd.read_csv(input, sep = ";", header = 0)
    print(raw_data.shape)
    print("Unificate species")
    raw_data.loc[raw_data["Species"].str.contains("HomoSapiens"),"Species"] = "human"
    raw_data.loc[raw_data["Species"].str.contains("MusMusculus"),"Species"] = "mouse"
    print("Unificate chains")
    raw_data.loc[raw_data["Gene"].str.contains("TRA"),"Gene"] = "alpha"
    raw_data.loc[raw_data["Gene"].str.contains("TRB"),"Gene"] = "beta"
    print("Check mhc names")
    raw_data['valid_i_name'] = raw_data['MHC A'].apply(lambda x: fix_mhc_name(x, 'alpha'))
    raw_data['valid_ii_name'] = raw_data['MHC B'].apply(lambda x: fix_mhc_name(x, 'beta'))
    raw_data["database"] = "VDJdb"
    raw_data['MHC class'] = raw_data['MHC class'].str.replace("MHC","")

    table_schema = {
            "id": "id",
            "Gene": "chain",
            "CDR3": "sequence",
            "Epitope": "epitope",
            "valid_i_name": "mhc_alpha",
            "valid_ii_name": "mhc_beta",
            "MHC class": "mhc_class",
            "Species": "host_species",
            "Epitope species": "epitope_species",
            "Epitope gene": "epitope_source",
            "V": "V",
            "D": "D",
            "J": "J",
            "database": "database"
        }


    print("Apply filters")
    clean_data = raw_data.\
            query("Reference.notna()").\
            query("Gene.isin(@CDR3_CHAINS)").\
            query("Species.isin(@HOST_SPECIES)").\
            query("CDR3.notna()").\
            query("CDR3.str.contains(@PROTEIN_CHECK)").\
            query("Epitope.notna()").\
            query("Epitope.str.contains(@PROTEIN_CHECK)").\
            query("valid_i_name.notna() and valid_ii_name.notna()")
    
    print(clean_data.shape)
    
    print("Get id")
    clean_data_with_id = calculate_receptor_id(clean_data,"receptor_id")
    clean_data_selected = clean_data_with_id.filter(items = list(table_schema.keys()), axis = 1).rename(columns = table_schema)
    print("Reshape")
    pivoted_part = clean_data_selected.filter(["id", "chain", "sequence","V","J"],axis=1).sort_values(by='V', na_position='last')
    annotation_part = clean_data_selected.filter(["id", "epitope", "mhc_alpha","mhc_beta","mhc_class","host_species","epitope_species","epitope_source","database"],axis=1)
    pivoted_part_pivot = pivoted_part.pivot_table(index = "id",columns='chain', values=['sequence', 'V', 'J'], aggfunc = 'first').reset_index()
    pivoted_part_pivot.columns = pivoted_part_pivot.columns.to_flat_index().str.join('_')
    pivoted_part_pivot.rename(columns = {"id_":"id", "sequence_alpha":"cdr3_alpha", "sequence_beta":'cdr3_beta'},inplace=True)
    wide = pd.merge(pivoted_part_pivot, annotation_part, on="id")
    wide["D_alpha"] = pd.NA
    wide["D_beta"] = pd.NA
    print("Save")
    wide.to_csv(os.path.join(output, "VDJdb_clean.csv"), sep = ";", index = False)

        

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