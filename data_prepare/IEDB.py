import pandas as pd
import traceback
from datetime import datetime
from collections import namedtuple
from mysql.connector import connect, Error as MysqlError
from database import Database

class IEDB(Database):

    def __init__(self, server: str, user: str, password:str, database:str = "IEDB"):
        assert database == "IEDB" or database == "CEDAR", f"Неподходящий объект для БД {database}"
        super().__init__(database=database)
        self.__server = server
        self.__user = user
        self.__password = password


    def acquire(self) -> pd.DataFrame | None:
        '''
        Download tables from IEDB and CEDAR connected with TCR and bind all into one raw big table by primary keys.
        The merge scheme see experiments/3_download_iedb_cedar.ipynb
        WARNING! Large tables store in RAM. Execute only cluster with slurm and huge RAM
        '''
        ready = None
        try:
            table_names = ["tcell","curated_epitope","object","epitope_object","epitope", "tcell_receptor","curated_receptor","distinct_receptor","distinct_chain","mhc_allele_restriction"]
            tables = {}
        
            #DOWNLOAD PART
            chunk_size = 100000
            with connect(host = self.__server,
                         user = self.__user,
                         password = self.__password,
                         database = self.database) as connection:
                for name in table_names:
                    print(name)
                    tbl = pd.DataFrame()
                    for chunk in pd.read_sql(f"SELECT * FROM {name};", con = connection, chunksize = chunk_size):
                        tbl = pd.concat([tbl, chunk], ignore_index = False)
                    tables[name] = tbl
                tables["organism_names"] = pd.read_sql("SELECT organism_id, name_txt FROM organism_names WHERE name_class = 'scientific name';", con = connection)
                
            #MERGE PART
            print("Merging...Can take a long time")
            #TCR INFO
            print("TCR information...")
            
            distinct_receptor_melt = pd.melt(tables["distinct_receptor"], 
                                            id_vars = ["distinct_receptor_id","receptor_type","species"], 
                                            value_vars = ["distinct_chain1_id","distinct_chain2_id"],
                                            var_name = "chain",
                                            value_name = "chain_id",
                                            ignore_index = True)
            distinct_receptor_melt_filt = distinct_receptor_melt[distinct_receptor_melt["chain_id"].notna()]
            chain_with_receptor = pd.merge(distinct_receptor_melt_filt, tables["distinct_chain"], left_on = "chain_id",right_on = "distinct_chain_id",suffixes = ("_receptor","_chain"))
            merged_curated_receptor = pd.merge(chain_with_receptor, tables["curated_receptor"], on = "distinct_receptor_id")
            tcr_information = pd.merge(merged_curated_receptor,  tables["tcell_receptor"], on = "curated_receptor_id")
            
            tcr_info_species_receptor_names = pd.merge(tcr_information, 
                                                    tables["organism_names"], 
                                                    left_on = "species_receptor", 
                                                    right_on = "organism_id").\
            drop(columns = ["organism_id"]).\
            rename(columns = {"name_txt":"species_name_receptor"})
            
            tcr_information_ready = pd.merge(tcr_info_species_receptor_names, 
                                            tables["organism_names"], 
                                            left_on = "species_chain", 
                                            right_on = "organism_id").\
            drop(columns = ["organism_id"]).\
            rename(columns = {"name_txt":"species_name_chain"})
            
            #Epitope INFO
            print("Epitope information...")
            
            epi_epi_object = pd.merge(tables["epitope"],tables["epitope_object"], on = "epitope_id")
            merged_epi_obj = pd.merge(tables["object"],epi_epi_object , on = "object_id")
            epitope_information = pd.merge(merged_epi_obj,tables["curated_epitope"], left_on = "object_id", right_on = "e_object_id")
            
            epitope_information_ready = pd.merge(epitope_information, 
                                                tables["organism_names"], 
                                                left_on = "source_organism_org_id", 
                                                right_on = "organism_id").\
            drop(columns = ["organism_id_x","organism_id_y","organism2_id","source_organism_org_id"]).\
            rename(columns = {"name_txt":"epitope_source"})

            #MHC INFO
            print("MHC information...")
            mhc_ready = pd.merge(tables["mhc_allele_restriction"], 
                                tables["organism_names"], 
                                left_on = "organism_ncbi_tax_id", 
                                right_on = "organism_id").\
            drop(columns = ["organism_id"]).\
            rename(columns = {"name_txt":"mhc_source"})
            
            #MERGING TOGETHER
            print("Merging...")
            
            tcell_with_names = pd.merge(tables['tcell'],
                                        tables["organism_names"],
                                        left_on = "h_organism_id", 
                                        right_on = "organism_id").\
            drop(columns = ["organism_id"]).\
            rename(columns = {"name_txt":"host_organism"})
            
            tcell_with_mhc = pd.merge(tcell_with_names,
                                    mhc_ready,
                                    on = "mhc_allele_restriction_id",
                                    suffixes = ("_tcell","_mhc"))
            
            tcell_tcr_mhc = pd.merge(tcell_with_mhc, 
                                    tcr_information_ready,
                                    on = "tcell_id")

            ready = pd.merge(tcell_tcr_mhc,
                            epitope_information_ready,
                            on = "curated_epitope_id")       
        except MysqlError as e:
            traceback.print_exc()
            print(f"Something with MySQL server connection {e}")
        finally:
            print("Finished")
            return ready
        
    def get_latest_update_date(self) -> namedtuple:
        Dates = namedtuple("Dates", [upload_on_server, last_update])
        if self._database == 'IEDB':
            return Dates(
                datetime.strptime("25.01.2026", "dd.mm.yyyy"),
                datetime.strptime("18.01.2026", "dd.mm.yyyy")
            )
        elif self._database == 'CEDAR':
            return Dates(
                datetime.strptime("25.01.2026", "dd.mm.yyyy"),
                datetime.strptime("06.01.2026", "dd.mm.yyyy")
            )
    
    def clean_tcr_epitope(self) -> pd.DataFrame:
        return super().clean_tcr_epitope()
    
    def clean_tcr_mhc(self) -> pd.DataFrame:
        return super().clean_tcr_mhc()