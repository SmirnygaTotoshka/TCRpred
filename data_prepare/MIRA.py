import uuid
import zipfile
import tempfile
import requests
import traceback
import pandas as pd
from database import Database
from datetime import datetime
from os.path import join
class MIRA(Database):

    def __init__(self, database=None):
        super().__init__('MIRA')
        self.__url = 'https://adaptivepublic.blob.core.windows.net/publishedproject-supplements/covid-2020/ImmuneCODE-MIRA-Release002.1.zip'
        self._tcr_epi_schema = {
                                  "TCR BioIdentity":"ReceptorID",
                                  "Database": "Database",
                                  "Chain": "Chain",
                                  "Species":"Species",
                                  "CDR3": "Structure",
                                  "epitope":"Activity",
                                  "ORF":"EpitopeProtein",
                                  "EpitopeOrganism":"EpitopeOrganism",
                                  "V":"V",
                                  "D":"D",
                                  "J":"J"}
        
        self._tcr_mhc_schema = {
                                  "TCR BioIdentity":"ReceptorID",
                                  "Database": "Database",
                                  "Chain": "Chain",
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

        self._species_cols = ["Species"]
        self._species_names = {"Human":"human",
                         "Mouse":"mouse"
                        }
        self._receptor_col = "receptor_id"
        

    def acquire(self):
        ready = None
        try:
            with tempfile.TemporaryDirectory() as tmpdirname:
                # Send a GET request to download the file
                response = requests.get(self.__url)
                print("Download...")
                response.raise_for_status()

                archive = join(tmpdirname, 'MIRA.zip')
                with open(archive, "wb") as file:
                    file.write(response.content)
                print(f"Downloaded {archive}")
                
                with zipfile.ZipFile(archive, 'r') as zip_ref:
                    zip_ref.extractall(tmpdirname)
                
                home = join(tmpdirname,"ImmuneCODE-MIRA-Release002.1")
                print("Convert...")
            
                #Subject metadata
                subject_metadata = pd.read_csv(join(home,"subject-metadata.csv"),encoding = "cp1251")
                hla_columns = subject_metadata.columns[8:]
                subject_metadata_melt = pd.melt(subject_metadata, 
                                        id_vars = ["Experiment","Subject","Cell Type"], 
                                        value_vars = hla_columns,
                                        var_name = "allele_type",
                                        value_name = "allele",
                                        ignore_index = True)
                subject_metadata_ready = subject_metadata_melt[subject_metadata_melt['allele'].notnull()]
                # Таблицы с данным невозможно присоединить по номеру участника, так как в них только информация об эксперименте. 
                # Поэтому участники, которые задействованы в нескольких экспериментах надо удалить для однозначного соответствия.
                times_participate = subject_metadata[['Subject', 'Experiment']].groupby('Subject').size().sort_values()
                multiple_times = list(times_participate.index[times_participate != 1])
                subject_metadata_ready = subject_metadata_ready.query("~Subject.isin(@multiple_times)")
                
                # Minigene data
                
                minigene_detail = pd.read_csv(join(home,"minigene-detail.csv"),encoding = "cp1251")
                minigene_detail[["CDR3","V","J"]] = minigene_detail["TCR BioIdentity"].str.split('\\+',expand=True)
                ids = {k:uuid.uuid4() for k in minigene_detail["TCR BioIdentity"]}
                minigene_detail["TCR BioIdentity"] = minigene_detail["TCR BioIdentity"].map(ids)
                
                # Peptide class I
                
                peptide_detail_i = pd.read_csv(join(home,"peptide-detail-ci.csv"),encoding = "cp1251")
                peptide_detail_i[["CDR3","V","J"]] = peptide_detail_i["TCR BioIdentity"].str.split('\\+',expand=True)
                ids = {k:uuid.uuid4() for k in peptide_detail_i["TCR BioIdentity"]}
                peptide_detail_i["TCR BioIdentity"] = peptide_detail_i["TCR BioIdentity"].map(ids)
                fn = lambda x: pd.Series([i for i in x.split(',')])
                epi = peptide_detail_i['Amino Acids'].apply(fn)
                peptide_detail_i_epi = pd.concat([peptide_detail_i,epi], axis=1)
                peptide_detail_i_melt = pd.melt(peptide_detail_i_epi, 
                                                id_vars = ["TCR BioIdentity","Experiment",'CDR3','V','J','ORF Coverage'], 
                                                value_vars = epi.columns,
                                                var_name = "epi_seq",
                                                value_name = "epitope",
                                                ignore_index = True)
                peptide_detail_i_ready = peptide_detail_i_melt[(peptide_detail_i_melt['epitope'].notnull()) & (peptide_detail_i_melt['CDR3'].notnull())]
                
                # Peptide class II
                
                peptide_detail_ii = pd.read_csv(join(home,"peptide-detail-cii.csv"),encoding = "cp1251")
                peptide_detail_ii[["CDR3","V","J"]] = peptide_detail_ii["TCR BioIdentity"].str.split('\\+',expand=True)
                ids = {k:uuid.uuid4() for k in peptide_detail_ii["TCR BioIdentity"]}
                peptide_detail_ii["TCR BioIdentity"] = peptide_detail_ii["TCR BioIdentity"].map(ids)
                epi = peptide_detail_ii['Amino Acids'].apply(fn)
                peptide_detail_ii_epi = pd.concat([peptide_detail_ii,epi], axis=1)
                peptide_detail_ii_melt = pd.melt(peptide_detail_ii_epi, 
                                        id_vars = ["TCR BioIdentity","Experiment",'CDR3','V','J','ORF Coverage'], 
                                        value_vars = epi.columns,
                                        var_name = "epi_seq",
                                        value_name = "epitope",
                                        ignore_index = True)
                peptide_detail_ii_ready = peptide_detail_ii_melt[(peptide_detail_ii_melt['epitope'].notnull()) & (peptide_detail_ii_melt['CDR3'].notnull())]
                
                #Join I class
                metadata_I = subject_metadata_ready.query("allele_type == 'I'")
                first = pd.merge(peptide_detail_i_ready, metadata_I, on = "Experiment").\
                drop(columns = ["epi_seq","Subject"]).\
                rename(columns = {"ORF Coverage":"ORF"})

                first["type"] = "peptide"

                #Join II class
                
                metadata_II = subject_metadata_ready.query("allele_type == 'II'")
                second = pd.merge(peptide_detail_ii_ready, metadata_II, on = "Experiment").\
                drop(columns = ["epi_seq","Subject"]).\
                rename(columns = {"ORF Coverage":"ORF"})

                second["type"] = "peptide"
                
                # Join minigene
                
                minigene = pd.merge(minigene_detail, metadata_II, on = "Experiment").\
                loc[:,["TCR BioIdentity","Experiment","CDR3","V","J","ORF","Amino Acid","Cell Type"]].\
                rename(columns={"Amino Acid":"epitope"})
                
                minigene["type"] = "minigene"      

                #Rbind all together
                ready = pd.concat([first,second, minigene],ignore_index = True)


        except requests.HTTPError as e:
            traceback.print_exc()
            print(f"Something with Internet {e}")

        finally:
            return ready
    
    def clean(self, raw_data: pd.DataFrame, dataset: str):
        self._dup_cols = ["Chain", "CDR3", "epitope","Species"] if dataset == "tcr-epitope" else ["Gene", "CDR3", self._mhc_fix_col, "Species"]
        return super().clean(raw_data, dataset)

    
    def get_latest_update_date(self):
        return datetime.strptime("04.08.2020","dd.mm.yyyy")