import pandas as pd
import traceback
import time
import uuid
import requests
import tempfile
from datetime import datetime
from database import Database
from tools import get_chrome_driver
from io import StringIO
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import WebDriverException

class McPAS(Database):

    def __init__(self, database=None):
        super().__init__("McPAS")
        self.__url = 'https://friedmanlab.weizmann.ac.il/McPAS-TCR/'
        self._mhc_fix_col = "MHC"
        self._tcr_epi_schema = {
                              "ReceptorID":"ReceptorID",
                              "Database": "Database",
                              "chain": "Chain",
                              "Species":"Species",
                              "cdr3_seq": "Structure",
                              "Epitope.peptide":"Activity",
                              "Antigen.protein":"EpitopeProtein",
                              "Pathology":"EpitopeOrganism",
                              "V":"V",
                              "D":"D",
                              "J":"J"}
        
        self._tcr_mhc_schema = {
                              "ReceptorID":"ReceptorID",
                              "Database": "Database",
                              "chain": "Chain",
                              "Species":"Species",
                              "cdr3_seq": "Structure",
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

        self._species_cols = ["Species"]
        self._species_names = {"Human":"human",
                         "Mouse":"mouse"
                        }
        self._receptor_col = "receptor_id"


    def acquire(self):
        download_button_path = "/html/body/div[1]/div/section/div/div/div[1]/div/div/div/div/div/div[2]/a"
        driver = None
        ready = None
        try:
            with tempfile.TemporaryDirectory() as tmpdirname:
                driver = get_chrome_driver(tmpdirname)
                driver.get(self.__url)
                time.sleep(10) # Нужно дать время, чтобы сервер сгенерировал сессионную ссылку
                link = driver.find_element(By.XPATH, download_button_path).get_attribute("href")
                print("Download...")
                data_result = requests.get(link)
                data_result.raise_for_status()
                print("Convert...")
                result = StringIO(data_result.text)
                ready = pd.read_csv(result, sep = ",")
        except requests.HTTPError as e:
            traceback.print_exc()
            print(f"Something with HTTP connection {e}")
        except WebDriverException as e:
            traceback.print_exc()
            print(f"Something with Selenium {e}")
        finally:
            if driver is not None:
                driver.quit()
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
        mcpas_last_update_path = "/html/body/div[1]/div/section/div/div/div[1]/div/div/div/div/div/div[2]/div[3]/p"
        driver = None
        last_update_date = None
        try:
            with tempfile.TemporaryDirectory() as tmpdirname:
                driver = get_chrome_driver(tmpdirname)
                driver.get(self.__url)
                last_update_string = WebDriverWait(driver, 40).until(
                    EC.presence_of_element_located((By.XPATH, mcpas_last_update_path))
                ).text
                clean_str = last_update_string.replace("The database was last updated on: ", "")
                last_update_date = datetime.strptime(clean_str, "%B %d, %Y")
        except WebDriverException as e:
            last_update_date = datetime.strptime("10.09.2022","%d.%m.%Y")
        finally:
            if driver is not None:
                driver.quit()
            return last_update_date