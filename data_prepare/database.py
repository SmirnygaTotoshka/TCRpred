import os
import mhcgnomes
import pandas as pd
from abc import ABC, abstractmethod

class Database(ABC):

    def __init__(self, database = None):
        super().__init__()
        self._final_columns_epitope = {}
        self._final_columns_mhc = {}
        self._database = database
        self._tcr_epitope_filters = {}
        self._tcr_mhc_filters = {}
        self._seq_pattern = f"^[ACDEFGHIKLMNPQRSTVWY]+$" # Use in queries
        self._included_species = ['human','mouse'] # Use in queries

    @abstractmethod
    def acquire(self):
        raise NotImplementedError("Should be used certain acquiring for the database")
    
    @abstractmethod
    def clean_tcr_epitope(self, raw_data: pd.DataFrame):
        raise NotImplementedError("Should be used certain cleaner for the database")
    
    @abstractmethod
    def clean_tcr_mhc(self, raw_data: pd.DataFrame):
        raise NotImplementedError("Should be used certain cleaner for the database")
    
    @abstractmethod
    def get_latest_update_date(self):
        raise NotImplementedError("Should be used certain cleaner for the database")
    
    def get_mhc_class(self, allele: str) -> str:
        a = mhcgnomes.parse(allele)
        return "I" if a.is_class1 else "II"

    def fix_allele(self, broken_allele: str) -> str:
        allele = mhcgnomes.parse(broken_allele)
        if allele.gene.species == "Homo sapiens":
            if len(allele.allele_fields) < 2:
                raise mhcgnomes.ParseError("Don`t need serotypes")
            elif len(allele.allele_fields) == 2:
                return allele.to_string()
            else:
                return allele.restrict_allele_fields(2, drop_annotations=True, drop_mutations=True).to_string()
        elif allele.gene.species == "Mus musculus":
            return allele.to_string()
        else:
            raise mhcgnomes.ParseError('Unknown species')
        

    def allele_to_compact_name(self, allele: str) -> str:
        return mhcgnomes.parse(allele).to_string()
    
    def save(self, data: pd.DataFrame, output: str) -> None:
        data.to_csv(output,sep = ";",index = False)