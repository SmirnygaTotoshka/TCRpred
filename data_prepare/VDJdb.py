import pandas as pd
import traceback
import requests
import tempfile
from datetime import datetime
from database import Database
from tools import get_chrome_driver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import WebDriverException

class VDJdb(Database):

    def __init__(self):
        super().__init__("VDJdb")
        self._tcr_epitope_schema = {
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
        
        self._tcr_mhc_schema = {
                                "receptor_id":"ReceptorID", 
                                "Database": "Database",
                                "Gene": "Chain",
                                "Species":"Species",
                                "CDR3": "Structure",
                                self._mhc_fix_col:"Activity",
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
        self._tcr_mhc_filters = {
            "concordant_species": "not (mhc_chain.str.contains('HLA') and Species == 'mouse') or not ((mhc_chain.str.contains('H2') and Species == 'human')",
            "canonical_mhc": "`MHC class` == 'I' or `MHC class` == 'II'",
            "not_b2m": "mhc_chain != 'B2M'",
        }  

        self._species_cols = ["Species"]
        self._species_names = {"HomoSapiens":"human",
                         "MusMusculus":"mouse"
                        }
        self._receptor_col = "receptor_id"
        self._mhc_value_columns = ["MHC A","MHC B"]


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
    
    def clean(self, raw_data: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame,pd.DataFrame]:     
        raw_data.loc[raw_data["Gene"] == "TRA","Gene"] = "alpha"
        raw_data.loc[raw_data["Gene"] == "TRB","Gene"] = "beta"

        self._dup_cols = ["Gene", "CDR3", "Epitope","Species"] if dataset == "tcr-epitope" else ["Gene", "CDR3", self._mhc_fix_col, "Species"]
        return super().clean(raw_data, dataset)

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
        