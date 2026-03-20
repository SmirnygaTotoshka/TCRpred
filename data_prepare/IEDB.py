import pandas as pd
import traceback
from datetime import datetime
from mysql.connector import connect, Error as MysqlError
from database import Database

class IEDB(Database):

    def __init__(self, server: str, user: str, password:str, database:str = "IEDB"):
        assert database == "IEDB" or database == "CEDAR", f"Неподходящий объект для БД {database}"
        super().__init__(database=database)
        self.__server = server
        self.__user = user
        self.__password = password
        self._tcr_epitope_schema = {
                                "curated_receptor_id":"ReceptorID",
                                "Database": "Database",
                                "chain_type": "Chain",
                                "host_organism":"Species",
                                "cdr3_seq": "Structure",
                                "linear_peptide_seq":"Activity",
                                "source_antigen_accession":"EpitopeProtein",
                                "epitope_source":"EpitopeOrganism",
                                "v_gene":"V",
                                "d_gene":"D",
                                "j_gene":"J"}
        
        self._tcr_mhc_schema = {
                                "curated_receptor_id":"ReceptorID",
                                "Database": "Database",
                                "chain_type": "Chain",
                                "host_organism":"Species",
                                "cdr3_seq": "Structure",
                                self._mhc_fix_col:"Activity",
                                "v_gene":"V",
                                "d_gene":"D",
                                "j_gene":"J"}
        
        self._filter_rules = {
            "has_reference": "reference_id.notna()",
            "ab_receptor": "receptor_type_x == 'alphabeta'",
            "ab_chains": "chain_type == 'alpha' or chain_type == 'beta'",
            "included_species": "species_name_receptor == species_name_chain and species_name_receptor.isin(@self._included_species)",
            "has_cdr3": "cdr3_seq.notna()",
            "valid_cdr3": "cdr3_seq.str.contains(@self._seq_pattern)",
            "is_epitope": "e_region_domain_flag == 'Exact Epitope' and ant_type == 'Epitope'",
            "has_epitope": "linear_peptide_seq.notna()",
            "valid_epitope": "linear_peptide_seq.str.contains(@self._seq_pattern)",
            "no_modif_epitope": "linear_peptide_modification.isna()",
            "is_linear": "disc_region.isna()",
            "is_positive": "as_char_value.str.contains('Positive')",
            "is_included_host": "host_organism == mhc_source and host_organism.isin(@self._included_species)",
            "is_complete_mhc": "restriction_level == 'complete molecule'",
            "without_mutations": "chain_i_mutation.isna() and chain_ii_mutation.isna()",
            "canonical_mhc": "`class` == 'I' or `class` == 'II'",
            "not_b2m": "@self._mhc_fix_col != 'Beta-2-microglobulin'"
        }#TODO

        self._species_cols = ["species_name_receptor","species_name_chain","mhc_source","host_organism"]
        self._species_names = {"Homo sapiens":"human",
                         "Mus musculus":"mouse"
                        }
        self._receptor_col = "curated_receptor_id"
        self._mhc_cols = ["chain_i_name","chain_ii_name"]

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
                         database = self._database) as connection:
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
            ready["Database"] = self._database
        except MysqlError as e:
            traceback.print_exc()
            print(f"Something with MySQL server connection {e}")
        finally:
            print("Finished")
            return ready
        
    def get_latest_update_date(self) -> tuple[datetime,datetime]:
        if self._database == 'IEDB':
            return (
                datetime.strptime("25.01.2026", "%d.%m.%Y"),
                datetime.strptime("18.01.2026", "%d.%m.%Y")
            )
        elif self._database == 'CEDAR':
            return (
                datetime.strptime("25.01.2026", "%d.%m.%Y"),
                datetime.strptime("06.01.2026", "%d.%m.%Y")
            )
        else:
            raise ValueError(f"Incorrect database {self._database}")
    
    def clean(self, raw_data: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame,pd.DataFrame]:
        self._dup_cols = ["chain_type", "cdr3_seq", "linear_peptide_seq","host_organism"] if dataset == "tcr-epitope" else ["chain_type","host_organism", "cdr3_seq", self._mhc_fix_col]
        return super().clean(raw_data, dataset)
