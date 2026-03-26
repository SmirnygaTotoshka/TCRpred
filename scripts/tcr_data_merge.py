#!/usr/bin/env python

import os
import click
import warnings
import mhcgnomes
import traceback
import numpy as np
import pandas as pd

output_option = click.option("-o","--output","output", required = True, type=str, help = "Output directory")
keep_option = click.option("-k","--keep","keep", is_flag=True, help = "Keep unneeded columns for PASS?")


@click.group()
def main():
    pass

@main.command()
@output_option
@click.argument('files', nargs=-1)
def simple(output, files):
    data = pd.concat([pd.read_csv(f, sep = ";", header = 0) for f in files],ignore_index = True).drop_duplicates(subset = ['cdr3_alpha', 'cdr3_beta', 'epitope', 'mhc_alpha', 'mhc_beta'])
    data['species_mhc_alpha'] = data['mhc_alpha'].apply(lambda x: mhcgnomes.parse(x).species.name)
    data['species_mhc_beta'] = data['mhc_beta'].apply(lambda x: mhcgnomes.parse(x).species.name)
    data_checked= data.query("species_mhc_alpha == species_mhc_beta").drop(columns = ['species_mhc_alpha','species_mhc_beta'])
    data_checked.to_csv(os.path.join(output, "TCR_info_simple_merge.csv"), sep = ";", index=False)
    print(data_checked.shape)
    
@main.command()
@keep_option
@output_option
@click.argument('files', nargs=-1)
def epitope(keep, output, files):
    if keep:
        use_cols = None
    else:
        use_cols = ['id','cdr3_alpha','cdr3_beta','epitope', 'host_species', 'epitope_species', 'epitope_source', 'database', 'mhc_alpha', 'mhc_beta']
    merged_data = pd.concat([pd.read_csv(f, sep = ";", header = 0, usecols = use_cols) for f in files],ignore_index = True).drop_duplicates(subset = ['cdr3_alpha', 'cdr3_beta', 'epitope', 'mhc_alpha', 'mhc_beta'])
    merged_data['species_mhc_alpha'] = merged_data['mhc_alpha'].apply(lambda x: mhcgnomes.parse(x).species.name)
    merged_data['species_mhc_beta'] = merged_data['mhc_beta'].apply(lambda x: mhcgnomes.parse(x).species.name)
    data_checked= merged_data.query("species_mhc_alpha == species_mhc_beta").drop(columns = ['species_mhc_alpha','species_mhc_beta'])
    
    alpha_human = data_checked.query("cdr3_alpha.notna() and host_species == 'human'").rename(columns = {'cdr3_alpha':'Structure', 'epitope': 'Activity'})
    alpha_human['chain'] = 'alpha'
    print(f'Alpha_human {alpha_human.shape}')
    alpha_human.to_csv(os.path.join(output, "cdr3alpha_epitope_human.csv"), sep = ";", index=False)
    
    beta_human = data_checked.query("cdr3_beta.notna() and host_species == 'human'").rename(columns = {'cdr3_beta':'Structure', 'epitope': 'Activity'})
    beta_human['chain'] = 'beta'
    print(f'Beta_human {beta_human.shape}')
    beta_human.to_csv(os.path.join(output, "cdr3beta_epitope_human.csv"), sep = ";", index=False)
    
    alpha_mouse = data_checked.query("cdr3_alpha.notna() and host_species == 'mouse'").rename(columns = {'cdr3_alpha':'Structure', 'epitope': 'Activity'})
    alpha_mouse['chain'] = 'alpha'
    print(f'Alpha_mouse {alpha_mouse.shape}')
    alpha_mouse.to_csv(os.path.join(output, "cdr3alpha_epitope_mouse.csv"), sep = ";", index=False)
    
    beta_mouse = data_checked.query("cdr3_beta.notna() and host_species == 'mouse'").rename(columns = {'cdr3_beta':'Structure', 'epitope': 'Activity'})
    beta_mouse['chain'] = 'beta'
    print(f'Beta_mouse {beta_mouse.shape}')
    beta_mouse.to_csv(os.path.join(output, "cdr3beta_epitope_mouse.csv"), sep = ";", index=False)

@main.command()
@keep_option
@output_option
@click.argument('files', nargs=-1)
def mhc(keep, output, files):
    if keep:
        use_cols = None
    else:
        use_cols = ['id','cdr3_alpha','cdr3_beta','epitope', 'host_species', 'epitope_species', 'epitope_source', 'database', 'mhc_alpha', 'mhc_beta']
    merged_data = pd.concat([pd.read_csv(f, sep = ";", header = 0, usecols = use_cols) for f in files],ignore_index = True).drop_duplicates(subset = ['cdr3_alpha', 'cdr3_beta', 'epitope', 'mhc_alpha', 'mhc_beta'])
    merged_data['species_mhc_alpha'] = merged_data['mhc_alpha'].apply(lambda x: mhcgnomes.parse(x).species.name)
    merged_data['species_mhc_beta'] = merged_data['mhc_beta'].apply(lambda x: mhcgnomes.parse(x).species.name)
    data_checked= merged_data.query("species_mhc_alpha == species_mhc_beta").drop(columns = ['species_mhc_alpha','species_mhc_beta'])
    merged_data['Activity'] = merged_data['mhc_alpha'] + "/" + merged_data['mhc_beta']
    
    alpha_human = data_checked.query("cdr3_alpha.notna() and host_species == 'human'").rename(columns = {'cdr3_alpha':'Structure'})
    alpha_human['chain'] = 'alpha'
    print(f'Alpha_human {alpha_human.shape}')
    alpha_human.to_csv(os.path.join(output, "cdr3alpha_mhc_human.csv"), sep = ";", index=False)
    
    beta_human = data_checked.query("cdr3_beta.notna() and host_species == 'human'").rename(columns = {'cdr3_beta':'Structure'})
    beta_human['chain'] = 'beta'
    print(f'Beta_human {beta_human.shape}')
    beta_human.to_csv(os.path.join(output, "cdr3beta_mhc_human.csv"), sep = ";", index=False)
    
    alpha_mouse = data_checked.query("cdr3_alpha.notna() and host_species == 'mouse'").rename(columns = {'cdr3_alpha':'Structure'})
    alpha_mouse['chain'] = 'alpha'
    print(f'Alpha_mouse {alpha_mouse.shape}')
    alpha_mouse.to_csv(os.path.join(output, "cdr3alpha_mhc_mouse.csv"), sep = ";", index=False)
    
    beta_mouse = data_checked.query("cdr3_beta.notna() and host_species == 'mouse'").rename(columns = {'cdr3_beta':'Structure'})
    beta_mouse['chain'] = 'beta'
    print(f'Beta_mouse {beta_mouse.shape}')
    beta_mouse.to_csv(os.path.join(output, "cdr3beta_mhc_mouse.csv"), sep = ";", index=False)    
    
if __name__ == "__main__":
    main()