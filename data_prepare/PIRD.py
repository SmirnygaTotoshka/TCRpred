import pandas as pd
import traceback
import time
import uuid
import requests
import tempfile
from datetime import datetime
from database import Database

class McPAS(Database):

    def __init__(self, database=None):
        super().__init__("PIRD")
        self.__url = 'https://ftp.cngb.org/pub/SciRAID/PIRD/TBAdb/TBAdb.xlsx'
        self._mhc_fix_col = "MHC"
        self._tcr_epi_schema = {
                                  "ReceptorID":"ReceptorID",
                                  "Database": "Database",
                                  "chain": "Chain",
                                  "Species":"Species",
                                  "cdr3_seq": "Structure",
                                  "Antigen.sequence":"Activity",
                                  "Antigen":"EpitopeProtein",
                                  "Disease.name":"EpitopeOrganism",
                                  "V":"V",
                                  "D":"D",
                                  "J":"J"}
        
        self._tcr_mhc_schema = {
                                  "ReceptorID":"ReceptorID",
                                  "Database": "Database",
                                  "chain": "Chain",
                                  "Species":"Species",
                                  "cdr3_seq": "Structure",
                                  "HLA":"Activity",
                                  "V":"V",
                                  "D":"D",
                                  "J":"J"}

        self._filter_rules = {
            "has_reference": "Reference.notna()",
            "ab_chains": "Gene == 'TRA' or Gene == 'TRB'",
            "included_species": "Species.isin(@self._included_species)",
            "has_cdr3": "CDR3.notna()",
            "valid_cdr3": "CDR3.str.contains(@self._seq_pattern)",
            "has_epitope": "Epitope.notna()",
            "valid_epitope": "Epitope.str.contains(@seq_pattern)",
        } #TODO

        self._species_cols = ["Species"]
        self._species_names = {"Human":"human",
                         "Mouse":"mouse"
                        }
        self._receptor_col = "receptor_id"


    def acquire(self):
        ready = None
        try:         
            with tempfile.NamedTemporaryFile() as t:
            # Send a GET request to download the file
                response = requests.get(self.__url)
                print("Download...")
                response.raise_for_status()
                with open(t.name, "wb") as file:
                    file.write(response.content)
                print(f"Downloaded {t.name}")
                time.sleep(30)  # Увеличьте время ожидания при необходимости
                print("Convert...")
                ready = pd.read_excel(t.name, sheet_name = "TCR-AB")
        except requests.HTTPError as e:
            traceback.print_exc()
            print(f"Something with HTTP connection {e}")    
        finally:
            return ready
    
    def clean(self, raw_data: pd.DataFrame, dataset: str):
        self._dup_cols = ["chain","cdr3_seq", "Epitope.peptide","Species"] if dataset == "tcr-epitope" else ["chain", "cdr3_seq", self._mhc_fix_col,"Species"]
        for k in self._tcr_epitope_schema.keys():
            assert k in raw_data.columns, "Incorrect column names. The raw data has been downloaded using McPAS.acquire method?"
        data_with_filters = self._preprocess(raw_data)
        dedup_data = data_with_filters.drop_duplicates(subset=self._dup_cols, ignore_index=True)
        
        return (
            dedup_data,
            self._apply_filters(dedup_data, dataset)
        )

    def _preprocess(self, data: pd.DataFrame) -> pd.DataFrame:
        mcpas_cdr3_reshaped = pd.melt(data, 
                              id_vars = data.columns[2:],                  
                              value_vars = ["CDR3.alpha.aa","CDR3.beta.aa"],
                              var_name = "chain",
                              value_name = "cdr3_seq",
                              ignore_index = True)
        fixed_species = self._fix_species(mcpas_cdr3_reshaped, self._species_cols, self._species_names)

        fixed_species.loc[fixed_species["chain"] == 'CDR3.alpha.aa','chain'] = 'alpha'
        fixed_species.loc[fixed_species["chain"] == 'CDR3.beta.aa','chain'] = 'beta'
        #TODO need test
        fixed_species.loc[fixed_species["chain"] == 'alpha','V'] = fixed_species.loc[:,"TRAV"]
        fixed_species.loc[fixed_species["chain"] == 'alpha','J'] = fixed_species.loc[:,"TRAJ"]
        fixed_species.loc[fixed_species["chain"] == 'beta','V'] = fixed_species.loc[fixed_species["chain"] == 'beta',"TRAV"]
        fixed_species.loc[fixed_species["chain"] == 'beta','D'] = fixed_species.loc[fixed_species["chain"] == 'beta',"TRBD"]
        fixed_species.loc[fixed_species["chain"] == 'beta','J'] = fixed_species.loc[fixed_species["chain"] == 'beta',"TRAJ"]

        fixed_mhc_alleles = self._fix_mhc_allele(fixed_species)
        calculated_filters = self._calculate_filters(fixed_mhc_alleles)
        return calculated_filters

    def _calculate_receptor_id(self, data: pd.DataFrame, column: str) -> pd.DataFrame:
        tmp = data.copy(deep=True)
        tmp[column] = [str(uuid.uuid4()) for _ in tmp.index]
        return tmp
        
    
    def get_latest_update_date(self):
        return datetime.strptime("24.07.2024","%d.%m.%Y")
