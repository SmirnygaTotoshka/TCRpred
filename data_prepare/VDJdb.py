import pandas as pd
import traceback
import uuid
import mhcgnomes
import requests
import tempfile
from datetime import datetime
from database import Database
from tools import get_chrome_driver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import WebDriverException

class VDJdb(Database):

    def __init__(self):
        super().__init__("VDJdb")
        self._final_columns_epitope = {
                                "receptor_id":"ReceptorID",  
                                "Database": "Database",
                                "Gene": "Chain",
                                "Species":"Species",
                                "CDR3": "Structure",
                                "Epitope":"Activity",
                                "Epitope gene":"EpitopeProtein",
                                "Epitope species":"EpitopeOrganism",
                                "V":"V",
                                "D":"D",
                                "J":"J"}
        
        self._final_columns_mhc = {
                                "receptor_id":"ReceptorID", 
                                "Database": "Database",
                                "Gene": "Chain",
                                "Species":"Species",
                                "CDR3": "Structure",
                                "mhc_chain":"Activity",
                                "V":"V",
                                "D":"D",
                                "J":"J"}

        self._tcr_epitope_filters = {
            "has_reference": "Reference.notna()",
            "ab_chains": "Gene == 'TRA' or Gene == 'TRB'",
            "included_species": "Species.isin(@self._included_species)",
            "has_cdr3": "CDR3.notna()",
            "valid_cdr3": "CDR3.str.contains(@self._seq_pattern)",
            "has_epitope": "Epitope.notna()",
            "valid_epitope": "Epitope.str.contains(@seq_pattern)",
        }
        self._tcr_mhc_filters = {
            "concordant_species": "not (mhc_chain.str.contains('HLA') and Species == 'mouse') or not ((mhc_chain.str.contains('H2') and Species == 'human')",
            "canonical_mhc": "`MHC class` == 'I' or `MHC class` == 'II'",
            "not_b2m": "mhc_chain != 'B2M'",
        }  


    def acquire(self):
        ready = None
        data_url = "https://vdjdb.com/api/database/search"
        meta_url = "https://vdjdb.com/api/database/meta"
        header = {"Content-Type": "application/json"}
        data = {"filters":[] }

        try:
            session = requests.Session()

            print("Get data...")
            data_result = session.post(data_url, json = data, headers = header)
            print("Get meta...")
            meta_result = session.get(meta_url)

            if data_result.status_code != requests.codes['ok'] or meta_result.status_code == requests.codes['ok']:
                data_result = session.post(data_url, json = data, headers = header, proxies={})
                meta_result = session.get(meta_url, proxies={})

            data_result.raise_for_status()
            meta_result.raise_for_status()

            print("Convert...")
            data_json = data_result.json()
            meta_json = meta_result.json()
            table = [d['entries'] for d in data_json["rows"]]
            columns = [c['title'] for c in meta_json['metadata']['columns']]
            ready = pd.DataFrame(table, columns = columns)
            ready["receptor_id"] = [d['metadata']['pairedID'] for d in data_json["rows"]]
        
        except requests.HTTPError as e:
            traceback.print_exc()
            print(f"Something with HTTP connection {e}")
        finally:
            print("Finished")
            return ready
    
    def clean_tcr_epitope(self, raw_data: pd.DataFrame):
        """
        For VDJdb use these filters:

        TCR-EPI:
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta
            3) epitopes with valid sequence.
            4) cdr3 with valid sequence
            5) receptor, host organism is human or mouse

        TCR-MHC:
            1) MHC first or second class
            3) Chains are known until certain protein or more info

        """
        cleaned_data = raw_data.copy(deep=True)
        # fix species. Experiment show that in all rows between all these cols values are identical
        cleaned_data.loc[cleaned_data["Species"] == "HomoSapiens","Species"] = "human"
        cleaned_data.loc[cleaned_data["Species"] == "MusMusculus","Species"] = "mouse"
        cleaned_data.loc[cleaned_data["Gene"] == "TRA","Gene"] = "alpha"
        cleaned_data.loc[cleaned_data["Gene"] == "TRB","Gene"] = "beta"

        for filter, eval in self._tcr_epitope_filters.items():
            cleaned_data[filter] = cleaned_data.eval(eval)        
        
        unique_id = cleaned_data["receptor_id"].unique()
        replacement = {k:str(uuid.uuid4()) for k in unique_id}
        cleaned_data["receptor_id"] = cleaned_data["receptor_id"].map(replacement)

        filter_names = list(self._tcr_epitope_filters.keys())
        cleaned_data["PASSED"] = cleaned_data.loc[:,filter_names].all()
        selected_columns = list(self._final_columns_epitope.keys())
        return (
            cleaned_data.query("PASSED").loc[:,selected_columns].rename(columns = self._final_columns_epitope),
            cleaned_data
        )

    
    def clean_tcr_mhc(self, raw_data: pd.DataFrame):
        """
        For VDJdb use these filters:

        TCR-EPI:
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta
            3) epitopes with valid sequence.
            4) cdr3 with valid sequence
            5) receptor, host organism is human or mouse
            +
            1) MHC first or second class
            3) Chains are known until certain protein or more info

        Nothing to return. Assign values to _tcr_epi and _tcr_mhc
        """
        _, data_with_filters = self.clean_tcr_epitope(raw_data)
        
        # melt all mhc alleles to one column
        mhc_value_columns = ["MHC A","MHC B"]
        iedb_mhc_reshaped = pd.melt(data_with_filters, 
                              id_vars = data_with_filters.columns[~data_with_filters.columns.isin(mhc_value_columns)],                  
                              value_vars = mhc_value_columns,
                              var_name = "mhc_chain_type",
                              value_name = "mhc_chain",
                              ignore_index = True)
        
        for i in iedb_mhc_reshaped.index:
            try:
                iedb_mhc_reshaped.loc[i, "mhc_chain"] = self.fix_allele(iedb_mhc_reshaped.loc[i, "mhc_chain"])
                iedb_mhc_reshaped.loc[i, "MHC class"] = self.get_mhc_class(iedb_mhc_reshaped.loc[i, "mhc_chain"])
                iedb_mhc_reshaped.loc[i, "success_parse"] = True
            except mhcgnomes.ParseError:
                iedb_mhc_reshaped.loc[i, "success_parse"] = False

        for filter, eval in self._tcr_epitope_filters.items():
            iedb_mhc_reshaped[filter] = iedb_mhc_reshaped.eval(eval)

        filter_names = list(self._tcr_epitope_filters.keys()).extend(list(self._tcr_mhc_filters.keys())).append("success_parse")
        iedb_mhc_reshaped["PASSED"] = iedb_mhc_reshaped.loc[:,filter_names].all()
        selected_columns = list(self._final_columns_mhc.keys())
        return (
            iedb_mhc_reshaped.query("PASSED").loc[:,selected_columns].rename(columns = self._final_columns_mhc),
            iedb_mhc_reshaped
        )
    
    def get_latest_update_date(self):
        vdjdb_url = 'https://vdjdb.cdr3.net/overview'
        vdjdb_last_update_xpath = "/html/body/application/div/overview/div/div/div/div/pre/code"
        last_update_date = None
        driver = None
        try:
            with tempfile.TemporaryDirectory() as tmpdirname:
                driver = get_chrome_driver(tmpdirname)
                driver.get(vdjdb_url)
                last_update_string = WebDriverWait(driver, 40).until(
                    EC.presence_of_element_located((By.XPATH, vdjdb_last_update_xpath))
                ).text
                clean_str = last_update_string.replace("Last updated on ", "")
                last_update_date = datetime.strptime(clean_str, "%d %B, %Y")
        except WebDriverException:
            last_update_date = datetime.strptime("28.12.2026","dd.mm.yyyy")
        finally:
            if driver is not None:
                driver.quit()
            return last_update_date
        