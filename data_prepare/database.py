import uuid
import mhcgnomes
import pandas as pd
from abc import ABC, abstractmethod
from tools import get_mhc_class

class Database(ABC):

    def __init__(self, database = None):
        super().__init__()
        self._tcr_epitope_schema = {}
        self._tcr_mhc_schema = {}
        self._database = database
        self._filter_rules = {}
        self._seq_pattern = f"^[ACDEFGHIKLMNPQRSTVWY]+$" # Use in queries
        self._included_species = ['human','mouse'] # Use in queries
        self._species_cols = []
        self._species_names = {
            "Homo sapiens":"human",
            "Mus musculus":"mouse"
        }
        self._receptor_col = ""
        self._mhc_cols = []
        self._mhc_class_col = "MHCclass"
        self._mhc_fix_col = "MHCchain"
        self._dup_cols = []

    @abstractmethod
    def acquire(self):
        raise NotImplementedError("Should be used certain acquiring for the database")
    
    @abstractmethod
    def clean(self, raw_data: pd.DataFrame, dataset: str):
        assert dataset == "tcr-epitope" or dataset == "tcr-mhc", f"Incorrect dataset {dataset}"
        for k in self._tcr_epitope_schema.keys():
            if k != self._mhc_fix_col: # _mhc_fix_col introduced only in cleaning
                assert k in raw_data.columns, f"Incorrect column names {k}. The raw data has been downloaded using Database.acquire method?"
        data_with_filters = self._preprocess(raw_data)
        dedup_data = data_with_filters.drop_duplicates(subset=self._dup_cols, ignore_index=True)
        
        return (
            dedup_data,
            self._apply_filters(dedup_data, dataset)
        )
    
    
    def _preprocess(self, data:pd.DataFrame) -> pd.DataFrame:
        fixed_species_names = self._fix_species(data, self._species_cols, self._species_names)
        fixed_receptor_id = self._calculate_receptor_id(fixed_species_names, self._receptor_col)
        mhc_reshaped = pd.melt(fixed_receptor_id, 
                              id_vars = fixed_receptor_id.columns[~fixed_receptor_id.columns.isin(self._mhc_cols)],               
                              value_vars = self._mhc_cols,
                              var_name = "mhc_chain_type",
                              value_name = self._mhc_fix_col,
                              ignore_index = True)
        fixed_mhc_allele = self._fix_mhc_allele(mhc_reshaped)
        calculated_filters = self._calculate_filters(fixed_mhc_allele)
        return calculated_filters

    
    @abstractmethod
    def get_latest_update_date(self):
        raise NotImplementedError("Should be used certain cleaner for the database")
    

    def _fix_mhc_allele(self, data: pd.DataFrame) -> pd.DataFrame:
        fixed_data = data.copy(deep=True)
        fixed_data["success_parse"] = False
        for i in fixed_data.index:
            allele = mhcgnomes.parse(fixed_data.loc[i, self._mhc_fix_col])
            if allele.gene.species == "Homo sapiens":
                if len(allele.allele_fields) < 2:
                    fixed_data.loc[i, "success_parse"] = False
                elif len(allele.allele_fields) == 2:
                    fixed_data.loc[i, self._mhc_fix_col] = allele.to_string()
                    fixed_data.loc[i,  self._mhc_class_col] = get_mhc_class(fixed_data.loc[i, self._mhc_fix_col])
                    fixed_data.loc[i, "success_parse"] = True
                else:
                    fixed_data.loc[i, self._mhc_fix_col] = allele.restrict_allele_fields(2, drop_annotations=True, drop_mutations=True).to_string()
                    fixed_data.loc[i,  self._mhc_class_col] = get_mhc_class(fixed_data.loc[i, self._mhc_fix_col])
                    fixed_data.loc[i, "success_parse"] = True
            elif allele.gene.species == "Mus musculus":
                fixed_data.loc[i, self._mhc_fix_col] = allele.to_string()
                fixed_data.loc[i,  self._mhc_class_col] = get_mhc_class(fixed_data.loc[i, self._mhc_fix_col])
                fixed_data.loc[i, "success_parse"] = True
            else:
                fixed_data.loc[i, "success_parse"] = False

        self._filter_rules["success_parse"] = "success_parse == True"
        return fixed_data
        

    def allele_to_compact_name(self, allele: str) -> str:
        return mhcgnomes.parse(allele).to_string()
    
    def save(self, data: pd.DataFrame, output: str) -> None:
        data.to_csv(output,sep = ";",index = False)

    def _fix_species(self, data: pd.DataFrame, columns: list[str], values: dict[str,str]) -> pd.DataFrame:
        cleaned_data = data.copy(deep=True)
        for col in columns:
            for k, v in values.items():
                assert v in self._included_species, f"{v} doesn`t include"
                cleaned_data.loc[cleaned_data[col].str.contains(k),col] = v
        return cleaned_data

    def _calculate_filters(self, data: pd.DataFrame) -> pd.DataFrame:
        cleaned_data = data.copy(deep=True)
        for f, e in self._filter_rules.items():
            cleaned_data[f] = cleaned_data.eval(e)
        return cleaned_data
    
    def _calculate_receptor_id(self, data: pd.DataFrame, column: str) -> pd.DataFrame:
        cleaned_data = data.copy(deep=True)
        unique_id = cleaned_data[column].unique()
        replacement = {k:str(uuid.uuid4()) for k in unique_id}
        cleaned_data[column] = cleaned_data[column].map(replacement)
        return cleaned_data
    
    def _apply_filters(self, data:pd.DataFrame, dataset: str) -> pd.DataFrame:
        final_columns = self._tcr_epitope_schema if dataset == "tcr-epitope" else self._tcr_mhc_schema
        cleaned_data = data.copy(deep=True)
        filter_names = list(self._filter_rules.keys())
        cleaned_data["PASSED"] = cleaned_data.loc[:,filter_names].all()
        selected_columns = list(final_columns.keys())
        return cleaned_data.query("PASSED").loc[:,selected_columns].rename(columns = final_columns)