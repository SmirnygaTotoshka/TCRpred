import pandas as pd
import traceback
import uuid
import mhcgnomes
from datetime import datetime
from mysql.connector import connect, Error as MysqlError
from database import Database

class IEDB(Database):
#TODO common rules for tcr-epi and tcr-mhc

    def __init__(self, server: str, user: str, password:str, database:str = "IEDB"):
        assert database == "IEDB" or database == "CEDAR", f"Неподходящий объект для БД {database}"
        super().__init__(database=database)
        self.__server = server
        self.__user = user
        self.__password = password
        self._final_columns_epitope = {
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
        
        self._final_columns_mhc = {
                                "curated_receptor_id":"ReceptorID",
                                "Database": "Database",
                                "chain_type": "Chain",
                                "host_organism":"Species",
                                "cdr3_seq": "Structure",
                                "mhc_chain":"Activity",
                                "v_gene":"V",
                                "d_gene":"D",
                                "j_gene":"J"}
        
        self._tcr_epitope_filters = {
            "has_reference": "reference_id.notna()",
            "ab_receptor": "receptor_type_x == 'alphabeta'",
            "ab_chains": "chain_type == 'alpha' or chain_type == 'beta'",
            "included_species": "species_name_receptor == species_name_chain and species_name_receptor.isin(@self._included_species)",
            "has_cdr3": "cdr3_seq.notna()",
            "valid_cdr3": "cdr3_seq.str.contains(@seq_pattern)",
            "is_epitope": "e_region_domain_flag == 'Exact Epitope' and ant_type == 'Epitope'",
            "has_epitope": "linear_peptide_seq.notna()",
            "valid_epitope": "linear_peptide_seq.str.contains(@seq_pattern)",
            "no_modif_epitope": "linear_peptide_modification.isna()",
            "is_linear": "disc_region.isna()",
            "is_positive": "as_char_value.str.contains('Positive')"
        }
        self._tcr_mhc_filters = {
            "is_included_host": "host_organism == mhc_source and host_organism.isin(@self._included_species)",
            "is_complete_mhc": "restriction_level == 'complete molecule'",
            "without_mutations": "chain_i_mutation.isna() and chain_ii_mutation.isna()",
            "canonical_mhc": "`class` == 'I' or `class` == 'II'",
            "not_b2m": "mhc_chain != 'Beta-2-microglobulin'",
        }      

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
        except MysqlError as e:
            traceback.print_exc()
            print(f"Something with MySQL server connection {e}")
        finally:
            print("Finished")
            return ready
        
    def get_latest_update_date(self) -> tuple[datetime,datetime]:
        if self._database == 'IEDB':
            return (
                datetime.strptime("25.01.2026", "dd.mm.yyyy"),
                datetime.strptime("18.01.2026", "dd.mm.yyyy")
            )
        elif self._database == 'CEDAR':
            return (
                datetime.strptime("25.01.2026", "dd.mm.yyyy"),
                datetime.strptime("06.01.2026", "dd.mm.yyyy")
            )
        else:
            raise ValueError(f"Incorrect database {self._database}")
    
    def clean_tcr_epitope(self, raw_data: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame]:
        """
        For IEDB and CEDAR use these filters:

        TCR-EPI:
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta
            3) full exact epitopes with valid sequence. it must be linear without modification and discontinious regions
            4) cdr3 with valid sequence
            5) receptor, host organism is human or mouse
            6) only positive results. negative sample is too small.

        Return tuple (filtered_data, total_data_with_filters_sresult)
        """
        cleaned_data = raw_data.copy(deep=True)
        species_columns_to_clean = ["species_name_receptor","species_name_chain","mhc_source","host_organism"]
        for col in species_columns_to_clean:
            cleaned_data.loc[cleaned_data[col].str.contains("Homo sapiens"),col] = "human"
            cleaned_data.loc[cleaned_data[col].str.contains("Mus musculus"),col] = "mouse"

        for filter, eval in self._tcr_epitope_filters.items():
            cleaned_data[filter] = cleaned_data.eval(eval)        
        
        unique_id = cleaned_data["curated_receptor_id"].unique()
        replacement = {k:str(uuid.uuid4()) for k in unique_id}
        cleaned_data["curated_receptor_id"] = cleaned_data["curated_receptor_id"].map(replacement)

        filter_names = list(self._tcr_epitope_filters.keys())
        cleaned_data["PASSED"] = cleaned_data.loc[:,filter_names].all()
        selected_columns = list(self._final_columns_epitope.keys())
        return (
            cleaned_data.query("PASSED").loc[:,selected_columns].rename(columns = self._final_columns_epitope),
            cleaned_data
        )
    
    def clean_tcr_mhc(self, raw_data: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame]:
        """
        For IEDB and CEDAR use these filters:

        TCR-MHC:
        same as tcr-epi
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta
            3) full exact epitopes with valid sequence. it must be linear without modification and discontinious regions
            4) cdr3 with valid sequence
            5) receptor, host organism is human or mouse
            6) only positive results. negative sample is too small.
        + additional rules
            7) MHC first or second class
            8) MHC without mutations
            9) Chains are known until certain protein or more info

        """
        _, data_with_filters = self.clean_tcr_epitope(raw_data)
        
        # melt all mhc alleles to one column
        mhc_value_columns = ["chain_i_name","chain_ii_name"]
        iedb_mhc_reshaped = pd.melt(data_with_filters, 
                              id_vars = data_with_filters.columns[~data_with_filters.columns.isin(mhc_value_columns)],                  
                              value_vars = mhc_value_columns,
                              var_name = "mhc_chain_type",
                              value_name = "mhc_chain",
                              ignore_index = True)
        
        for i in iedb_mhc_reshaped.index:
            try:
                iedb_mhc_reshaped.loc[i, "mhc_chain"] = self.fix_allele(iedb_mhc_reshaped.loc[i, "mhc_chain"])
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
