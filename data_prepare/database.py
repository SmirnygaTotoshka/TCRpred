import os
import mhcgnomes
from abc import ABC, abstractmethod

class Database(ABC):

    def __init__(self, database = None):
        super().__init__()
        self._final_columns_epitope = ["ReceptorID", 
                                      "Database", 
                                      "Chain",
                                      "Species", 
                                      "Structure",
                                      "Activity",
                                      "EpitopeProtein", 
                                      "EpitopeOrganism", 
                                      "V",
                                      "D",
                                      "J"]
        self._final_columns_mhc = ["ReceptorID", 
                                      "Database", 
                                      "Chain",
                                      "Species", 
                                      "Structure",
                                      "Activity",
                                      "V",
                                      "D",
                                      "J"]
        self._database = database

    @abstractmethod
    def acquire(self):
        raise NotImplementedError("Should be used certain acquiring for the database")
    
    @abstractmethod
    def clean_tcr_epitope(self):
        raise NotImplementedError("Should be used certain cleaner for the database")
    
    @abstractmethod
    def clean_tcr_mhc(self):
        raise NotImplementedError("Should be used certain cleaner for the database")
    
    @abstractmethod
    def get_latest_update_date(self):
        raise NotImplementedError("Should be used certain cleaner for the database")
    
    def get_mhc_class(self, allele: str) -> str:
        a = mhcgnomes.parse(allele)
        return "I" if a.is_class1 else "II"

    def fix_allele(self, broken_allele: str) -> str:
        return mhcgnomes.parse(broken_allele).to_string()

    def allele_to_compact_name(self, allele: str) -> str:
        return mhcgnomes.parse(allele).to_string()
    
    def save(self, data: pd.DataFrame, output: str) -> None:
        data.to_csv(output,sep = ";",index = False)