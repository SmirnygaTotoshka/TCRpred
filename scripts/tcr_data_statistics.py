#!/usr/bin/env python

import sqlite3
import click
import pandas as pd
import numpy as np
from parallelbar import progress_map
from functools import reduce

class TCRDatasetStatistics:
    """
    В конструкторе в переменные записывается статистика для дальнейшего перевода и сохранения в JSON. JSON в дальнейшем будет использоваться в дашбоарде.
    """
    
    @staticmethod
    def process_sequence(args):
        """
        функция для параллельного расчета матрицы встречаемости аминокислотных остатков в определенных позициях
        """
        seq, target_max = args
        zeros = np.zeros(shape = (20, target_max))
        matrix = pd.DataFrame(zeros, index = list("ACDEFGHIKLMNPQRSTVWY"), columns = list(range(1, target_max+1)))
        for index, char in enumerate(seq):
            matrix.loc[char, index+1] += 1
        return matrix

    @staticmethod
    def get_mhc_class(allele):
        """
        функция определяющая класс mhc аллели
        """
        if "H2-" in allele:
            if "H2-I" in allele:
                return "II"
            else:
                return "I"
        elif "HLA-" in allele:
            if "D" in allele:
                return "II"
            else:
                return "I"
        else:
            raise ValueError(f"Unknown allele {allele}")
    
    
    def __init__(self, dataset, contain_mhc=False, contain_epitope = True, n_cpu = 10, process_timeout=300):   
        self.num_paired_receptors = self.get_num_paired_receptors(dataset)
        self.sources = dataset["Database"].value_counts(dropna=False)
        self.chains = dataset["Chain"].value_counts(dropna=False)
        self.host_species = dataset["Species"].value_counts(dropna=False)
        self.num_na = dataset.isna().sum()
        self.unique_receptors = len(dataset["ReceptorID"].unique())
        self.unique_structures = len(dataset["Structure"].unique())
        self.unique_activities = len(dataset["Activity"].unique())
        self.host_chain_stat = dataset.loc[:,["Species","Chain"]].value_counts(dropna=False).reset_index()
        self.sample_size_by_activity = dataset["Activity"].value_counts(dropna=False)
        
        self.structure_logo = self.calc_logo(dataset, "Structure", n_cpu, process_timeout)
        
        if contain_epitope:
            self.activity_logo = self.calc_logo(dataset, column = "Activity")
        else:
            self.activity_logo = None
        
        if contain_mhc:
            self.mhc_stat = self.get_stat_by_mhc_class(dataset, column_mhc = "Activity")
        else:
            self.mhc_stat = None

    def calc_logo(self, dataset, column, n_cpu = 10, process_timeout=300):
        target_max = dataset[column].str.len().max()
        sequences = dataset[column].tolist() 
        args = [(s, target_max) for s in sequences]
        mm = progress_map(TCRDatasetStatistics.process_sequence, args, n_cpu = n_cpu, process_timeout=process_timeout, return_failed_tasks = False)
        sum_matrix = reduce(lambda x, y: x + y, mm) # TODO suspicious
        return sum_matrix / len(dataset.index)
    
    def get_stat_by_mhc_class(self, dataset, column_mhc, n_cpu = 10, process_timeout=300):
        classes = progress_map(TCRDatasetStatistics.get_mhc_class, dataset[column_mhc].tolist(), n_cpu = n_cpu, process_timeout=process_timeout, return_failed_tasks = False)
        return pd.DataFrame(classes, columns=['MHC'])["MHC"].value_counts(dropna=False)
    
    def get_num_paired_receptors(self, dataset):
        chains_by_id = dataset[["ReceptorID","Chain"]].value_counts().reset_index()[["ReceptorID","Chain"]]
        counts = chains_by_id["ReceptorID"].value_counts()
        return len(counts[counts == 2])
    
    def to_sqlite(self, db:str) -> None:
        with sqlite3.connect(db) as conn:
            cursor = conn.cursor()

            conn.execute("BEGIN")

            cursor.execute("""
                CREATE TABLE SingleStat (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL,
                    value REAL NOT NULL
                );
            """)

            self.structure_logo.to_sql("StructureLogo", conn)
            if self.activity_logo is not None:
                self.activity_logo.to_sql("ActivityLogo", conn)
            

            cursor.executemany("""
            INSERT INTO SingleStat (name, value) VALUES (?, ?)
            """, [("Количество парных рецепторов",self.num_paired_receptors), 
                ("Количество уникальных рецепторов",self.unique_receptors),
                ("Количество уникальных структур",self.unique_structures),
                ("Количество уникальных активностей",self.unique_activities)
                ]
            )  

            self.sources.to_sql("DataSources", conn)
            self.chains.to_sql("TCRChains", conn)
            self.host_species.to_sql("Hosts", conn)
            self.num_na.to_sql("NumNA", conn)
            self.host_chain_stat.to_sql("HostChain", conn, index = False)
            self.sample_size_by_activity.to_sql("SampleSizes", conn)
            
            if self.mhc_stat is not None:
                self.mhc_stat.to_sql("MHCClasses", conn)
    
@click.command
@click.option("-i","--input","input", required = True, type=str,help = "Path to input csv file. Expect csv with ';' as separator")
@click.option("-o","--output","output", required = True, type=str,help = "Output csv with ';' as separator")
@click.option("-m","--mhc","contain_mhc", required = False, default = False,  is_flag=True, help = "TCR-MHC input dataset")
@click.option("-e","--epitope","contain_epitope", required = False, default = False, is_flag=True, help = "TCR-Epi input dataset")
@click.option("-n","--threads","threads", required = False, default = 10, type=int,help = "Number of CPU")
@click.option("-t","--timeout","timeout", required = False, default = 300, type=int,help = "Multithreads operation timeout") 
def main(input, output, contain_mhc, contain_epitope, threads, timeout):
    dataset = pd.read_csv(input, sep = ";", header = 0)
    tcr_stat = TCRDatasetStatistics(dataset, contain_mhc, contain_epitope, threads, timeout)
    tcr_stat.to_sqlite(output)
    print("Stat SUCCESS")


if __name__ == "__main__":
    main()