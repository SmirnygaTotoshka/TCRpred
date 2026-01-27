import click
import os
import uuid
import pandas as pd
import numpy as np
import warnings
from abc import ABC, abstractmethod

warnings.filterwarnings("ignore")

class Cleaner(ABC):
    """
    All inherited classes directed to clean datasets for models cdr3->epitope cdr3->mhc. They don`t clean
    side fields (V,J genes for example). 
    """
    def __init__(self, input, output):
        self._raw_data = pd.read_csv(input, sep = ";", header = 0)
        self._output = output
        self._database = None

        self._tcr_epi = None
        self._tcr_mhc = None

        protein_alphabet = "ACDEFGHIKLMNPQRSTVWY"
        self._seq_pattern = f"^[{protein_alphabet}]+$" # Use in queries
        self._included_species = ['human','mouse'] # Use in queries
        
        self._tcr_epi_schema = {}#TODO
        self._tcr_mhc_schema = {}#TODO
        
    @staticmethod
    def trim(data: pd.Series):
        replaced_non_breaking_spaces = data.replace("","")
        
    @abstractmethod
    def clean(self):
        raise NotImplementedError("Should be used certain cleaner for the database") 
        
    @abstractmethod
    def _unificate(self):
        self._tcr_epi = self._convertToSchema(self._tcr_epi, self._tcr_epi_schema)
        self._tcr_mhc = self._convertToSchema(self._tcr_mhc, self._tcr_mhc_schema)

    def _convertToSchema(self, tbl: pd.DataFrame, schema: dict) -> pd.DataFrame:
        old_col_names = list(schema.keys())
        return tbl.loc[:,old_col_names].\
    rename(columns = schema).\
    apply(lambda x: x.str.strip())
    
    def _modificate_receptor_id(self, tbl: pd.DataFrame, column_id: str):
        unique_id = tbl[column_id].unique()
        replacement = {k:str(uuid.uuid4()) for k in unique_id}
        tbl[column_id] = tbl[column_id].map(replacement)
    
    def save(self):
        if self._tcr_epi is None or self._tcr_mhc is None:
            raise ValueError("Cleaned datasets is not defined. Please call Cleaner.clean() first.")
        else:
            self._tcr_epi.to_csv(os.path.join(self._output,f"{self._database}_tcr_epi_cleaned.csv"), sep = ";", index = False, header = True)
            self._tcr_mhc.to_csv(os.path.join(self._output,f"{self._database}_tcr_mhc_cleaned.csv"), sep = ";", index = False, header = True)
            
    def trimHLA(self, allele_name: str):
        num_points = allele_name.count(":")
        if num_points == 0:
            #Мы хотим конкретный белковый продукт,
            #если запись только до серотипа,то нам она не нужна
            return None
        elif num_points == 1:
            return allele_name
        else:
            return self.trimHLA(allele_name[0:allele_name.rfind(":")])
            
    

class IEDBCleaner(Cleaner):

    def __init__(self, input, output, database):
        super().__init__(input, output)
        self._database = database
        self._raw_data["Database"] = self._database
        self._tcr_epi_schema = {
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
                                "mhc_chain":"Activity",
                                "v_gene":"V",
                                "d_gene":"D",
                                "j_gene":"J"}

    def clean(self):
        """
        For IEDB and CEDAR use these filters:

        TCR-EPI:
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta
            3) full exact epitopes with valid sequence. it must be linear without modification and discontinious regions
            4) cdr3 with valid sequence
            5) receptor, host organism is human or mouse
            6) only positive results. negative sample is too small.

        TCR-MHC:
            1) MHC first or second class
            2) MHC without mutations
            3) Chains are known until certain protein or more info

        Nothing to return. Assign values to _tcr_epi and _tcr_mhc
        """
        # fix species. Experiment show that in all rows between all these cols values are identical
        species_columns_to_clean = ["species_name_receptor","species_name_chain","mhc_source","host_organism"]
        for col in species_columns_to_clean:
            self._raw_data.loc[self._raw_data[col].str.contains("Homo sapiens"),col] = "human"
            self._raw_data.loc[self._raw_data[col].str.contains("Mus musculus"),col] = "mouse"
        
        self._tcr_epi = self._raw_data.\
            query("reference_id.notna()").\
            query("receptor_type_x == 'alphabeta'").\
            query("chain_type == 'alpha' or chain_type == 'beta'").\
            query("species_name_receptor == species_name_chain and species_name_receptor.isin(@self._included_species)").\
            query("cdr3_seq.notna()").\
            query("cdr3_seq.str.contains(@self._seq_pattern)").\
            query("e_region_domain_flag == 'Exact Epitope'").\
            query("linear_peptide_seq.notna()").\
            query("linear_peptide_seq.str.contains(@self._seq_pattern)").\
            query("linear_peptide_modification.isna()").\
            query("disc_region.isna()").\
            query("as_char_value.str.contains('Positive')").\
            query("ant_type == 'Epitope'")         
        
        self._modificate_receptor_id(self._tcr_epi, "curated_receptor_id")
        
        #-------------------------------------------------------------------------------------
        # TCR-EPI STOP, NEXT STEP FOR TCR-MHC                                                # 
        #-------------------------------------------------------------------------------------

        # melt all mhc alleles to one column
        mhc_value_columns = ["chain_i_name","chain_ii_name"]
        iedb_mhc_reshaped = pd.melt(self._tcr_epi, 
                              id_vars = self._tcr_epi.columns[~self._tcr_epi.columns.isin(mhc_value_columns)],                  
                              value_vars = mhc_value_columns,
                              var_name = "mhc_chain_type",
                              value_name = "mhc_chain",
                              ignore_index = True)
        
        self._tcr_mhc = iedb_mhc_reshaped.\
            query("host_organism == mhc_source and host_organism.isin(@self._included_species)").\
            query("restriction_level == 'complete molecule'").\
            query("chain_i_mutation.isna() and chain_ii_mutation.isna()").\
            query("`class` == 'I' or `class` == 'II'").\
            query("mhc_chain != 'Beta-2-microglobulin'").\
            query("mhc_chain.str.count(':') >= 1 or mhc_source == 'mouse'")#mouse MHC are written very clearly
        
        #-------------------------------------------------------------------------------------
        # TCR-MHC STOP, Make table to the one from                                           # 
        #-------------------------------------------------------------------------------------
       
        self._unificate()
        
    #---------------------------------------------------------------------------------------------------
                             
    def _unificate(self):
        super()._unificate()
        #trim HLA
        for i in self._tcr_mhc.index:
            if self._tcr_mhc.loc[i, "Species"] == 'human':
                self._tcr_mhc.loc[i, "Activity"] = self.trimHLA(self._tcr_mhc.loc[i, "Activity"])


class VDJdbCleaner(Cleaner):

    def __init__(self, input, output):
        super().__init__(input, output)
        self._database = "VDJdb"
        self._raw_data["Database"] = self._database
        self._raw_data["D"] = pd.NA
        self._tcr_epi_schema = {
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
                              "mhc_chain":"Activity",
                              "V":"V",
                              "D":"D",
                              "J":"J"}

    def clean(self):
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

        Nothing to return. Assign values to _tcr_epi and _tcr_mhc
        """
        # fix species. Experiment show that in all rows between all these cols values are identical
        self._raw_data.loc[self._raw_data["Species"] == "HomoSapiens","Species"] = "human"
        self._raw_data.loc[self._raw_data["Species"] == "MusMusculus","Species"] = "mouse"

        self._tcr_epi = self._raw_data.\
            query("Reference.notna()").\
            query("Gene == 'TRA' or Gene == 'TRB'").\
            query("Species.isin(@self._included_species)").\
            query("CDR3.str.contains(@self._seq_pattern)").\
            query("Epitope.str.contains(@self._seq_pattern)")
        
        self._modificate_receptor_id(self._tcr_epi, "receptor_id")
        
        #-------------------------------------------------------------------------------------
        # TCR-EPI STOP, NEXT STEP FOR TCR-MHC                                                # 
        #-------------------------------------------------------------------------------------
        
        mhc_value_columns = ["MHC A","MHC B"]
        vdjdb_mhc_reshaped = pd.melt(self._tcr_epi, 
                              id_vars = self._tcr_epi.columns[~self._tcr_epi.columns.isin(mhc_value_columns)],                  
                              value_vars = mhc_value_columns,
                              var_name = "chain",
                              value_name = "mhc_chain",
                              ignore_index = True)

        self._tcr_mhc = vdjdb_mhc_reshaped.\
            query("`MHC class` == 'MHCI' or `MHC class` == 'MHCII'").\
            query("mhc_chain != 'B2M'").\
        query("mhc_chain.str.count(':') >= 1 or Species == 'mouse'").\
        query("~(mhc_chain.str.contains('HLA') and Species == 'mouse')") #discovered mouse rows with HLA
        
        #-------------------------------------------------------------------------------------
        # TCR-MHC STOP, Make table to the one from                                           # 
        #-------------------------------------------------------------------------------------
       
        self._unificate()
    #---------------------------------------------------------------------------------------------------
    def _unificate(self):
        
        # fix tcr chain type
        
        self._tcr_epi.loc[self._tcr_epi["Gene"] == "TRA","Gene"] = "alpha"
        self._tcr_epi.loc[self._tcr_epi["Gene"] == "TRB","Gene"] = "beta"
        
        self._tcr_mhc.loc[self._tcr_mhc["Gene"] == "TRA","Gene"] = "alpha"
        self._tcr_mhc.loc[self._tcr_mhc["Gene"] == "TRB","Gene"] = "beta"
        
        # fix MHC names
        
        for i in self._tcr_mhc.index:
            if 'H-2' in self._tcr_mhc.loc[i,'mhc_chain']:#Vectorize format not works, Have to make such of this.
                self._tcr_mhc.loc[i,'mhc_chain'] = self._tcr_mhc.loc[i,'mhc_chain'].replace("H-2","H2-").strip()
            if 'H2-Kb' in self._tcr_mhc.loc[i,'mhc_chain'] or 'H2-Kb' in self._tcr_mhc.loc[i,'mhc_chain'] or 'H2-KB' in self._tcr_mhc.loc[i,'mhc_chain']:
                self._tcr_mhc.loc[i,'mhc_chain'] = 'H2-Kb'
            if self._tcr_mhc.loc[i,'Species'] == "human" and not self._tcr_mhc.loc[i,'mhc_chain'].startswith('HLA-'):
                self._tcr_mhc.loc[i,'mhc_chain'] = f"HLA-{self._tcr_mhc.loc[i,'mhc_chain']}"
        self._tcr_mhc.loc[self._tcr_mhc['mhc_chain'].str.contains('I-Ab'),'mhc_chain'] = 'H2-IAb alpha'
        self._tcr_mhc.loc[self._tcr_mhc['mhc_chain'].str.contains('Aa'),'mhc_chain'] = 'H2-IAa alpha'
        self._tcr_mhc.loc[self._tcr_mhc['mhc_chain'].str.contains('Eb1'),'mhc_chain'] = 'H2-IEb1 alpha'
        
        self._tcr_mhc.loc[self._tcr_mhc['mhc_chain'].str.contains('I-Ab'),'mhc_chain'] = 'H2-IAb beta'
        self._tcr_mhc.loc[self._tcr_mhc['mhc_chain'].str.contains('Aa'),'mhc_chain'] = 'H2-IAa beta'
        self._tcr_mhc.loc[self._tcr_mhc['mhc_chain'].str.contains('Eb1'),'mhc_chain'] = 'H2-IEb1 beta'
        
        super()._unificate()
        
        #trim HLA
        for i in self._tcr_mhc.index:
            if self._tcr_mhc.loc[i, "Species"] == 'human':
                self._tcr_mhc.loc[i, "Activity"] = self.trimHLA(self._tcr_mhc.loc[i, "Activity"])
    #---------------------------------------------------------------------------------------------------
    
    def _modificate_receptor_id(self, tbl: pd.DataFrame, column_id: str):
        unique_id = tbl[column_id].unique()
        # receptor_id == 0 means that receptor isn`t paired
        # they need own id, not just replace
        index = np.argwhere(unique_id == 0)
        paired_id = np.delete(unique_id, index)
        replacement = {k:str(uuid.uuid4()) for k in paired_id}
        tbl[column_id] = tbl[column_id].map(replacement)
        
        for i in tbl.index:
            if pd.isna(tbl.loc[i, column_id]):
                tbl.loc[i,column_id] = str(uuid.uuid4())
        
        
class McPASCleaner(Cleaner):

    def __init__(self, input, output):
        super().__init__(input, output)
        self._database = "McPAS"
        self._raw_data["Database"] = self._database
        self._raw_data["ReceptorID"] = [str(uuid.uuid4()) for _ in self._raw_data.index]
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
                              "MHC":"Activity",
                              "V":"V",
                              "D":"D",
                              "J":"J"}


    def clean(self):
        """
        For McPAS-TCR use these filters:

        TCR-EPI:
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta
            3) epitopes with valid sequence.
            4) cdr3 with valid sequence
            5) receptor, host organism is human or mouse

        TCR-MHC:
            1) MHC first or second class
            3) Chains are known until certain protein or more info

        Nothing to return. Assign values to _tcr_epi and _tcr_mhc
        """
        
        mcpas_cdr3_reshaped = pd.melt(self._raw_data, 
                              id_vars = self._raw_data.columns[2:],                  
                              value_vars = ["CDR3.alpha.aa","CDR3.beta.aa"],
                              var_name = "chain",
                              value_name = "cdr3_seq",
                              ignore_index = True)
        
        # fix species. Experiment show that in all rows between all these cols values are identical
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["Species"] == 'Human','Species'] = 'human'
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["Species"] == 'Mouse','Species'] = 'mouse'

        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'CDR3.alpha.aa','chain'] = 'alpha'
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'CDR3.beta.aa','chain'] = 'beta'
        #TODO need test
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'alpha','V'] = mcpas_cdr3_reshaped.loc[:,"TRAV"]
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'alpha','J'] = mcpas_cdr3_reshaped.loc[:,"TRAJ"]
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'beta','V'] = mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'beta',"TRAV"]
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'beta','D'] = mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'beta',"TRBD"]
        mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'beta','J'] = mcpas_cdr3_reshaped.loc[mcpas_cdr3_reshaped["chain"] == 'beta',"TRAJ"]

        self._tcr_epi = mcpas_cdr3_reshaped.\
            query("`PubMed.ID`.notna()").\
            query("chain == 'alpha' or chain == 'beta'").\
            query("Species.isin(@self._included_species)").\
            query("cdr3_seq.notna()").\
            query("cdr3_seq.str.contains(@self._seq_pattern)").\
            query("`Epitope.peptide`.notna()").\
            query("`Epitope.peptide`.str.contains(@self._seq_pattern)")
        
        #-------------------------------------------------------------------------------------
        # TCR-EPI STOP, NEXT STEP FOR TCR-MHC                                                # 
        #-------------------------------------------------------------------------------------

        self._tcr_mhc = self._tcr_epi.\
        query("MHC.notna()").\
            query("MHC.str.count(':') >= 1 or Species == 'mouse'").\
        query("MHC != 'DR3*02:02'")
        
        #-------------------------------------------------------------------------------------
        # TCR-MHC STOP, Make table to the one from                                           # 
        #-------------------------------------------------------------------------------------
       
        self._unificate()
    #---------------------------------------------------------------------------------------------------
    
    def _unificate(self):
        
        # fix HLA names
        
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "HLA-A2:01","MHC"] = "HLA-A*02:01"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "DRB1*15:03","MHC"] = "HLA-DRB1*15:03"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "DRB1*04:01","MHC"] = "HLA-DRB1*04:01"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "DRB1*04:01","MHC"] = "HLA-DRB1*04:01"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "DRB1*04:01","MHC"] = "HLA-DRB1*04:01"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "HLA-Cw* 16:01","MHC"] = "HLA-C*16:01"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "HLA-A*2:01","MHC"] = "HLA-A*02:01"
        
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "H2-db","MHC"] = "H2-Db"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "H2-kb","MHC"] = "H2-Kb"
        self._tcr_mhc.loc[self._tcr_mhc["MHC"] == "H2-Db","MHC"] = "H2-Db"
        # fix mouse H2 names
        for i in self._tcr_mhc.index:
            if 'H-2' in self._tcr_mhc.loc[i,"MHC"]:
                self._tcr_mhc.loc[i,"MHC"] = self._tcr_mhc.loc[i,"MHC"].replace("H-2","H2-")
            if 'H2-db' in self._tcr_mhc.loc[i,"MHC"] or 'H2-Db' in self._tcr_mhc.loc[i,"MHC"]:
                self._tcr_mhc.loc[i,"MHC"] = "H2-Db"
            if 'H2-kb' in self._tcr_mhc.loc[i,"MHC"]:
                self._tcr_mhc.loc[i,"MHC"] = "H2-Kb"
        included_mouse_mhc = ["H2-Kb","H2-Db","H2-Kd"]
        self._tcr_mhc = self._tcr_mhc.query("MHC in @included_mouse_mhc or Species =='human'")
        
        super()._unificate()
        
        #trim HLA
        for i in self._tcr_mhc.index:
            if self._tcr_mhc.loc[i, "Species"] == 'human':
                self._tcr_mhc.loc[i, "Activity"] = self.trimHLA(self._tcr_mhc.loc[i, "Activity"])
                

class PIRDCleaner(Cleaner):

    def __init__(self, input, output):
        super().__init__(input, output)
        self._database = "PIRD"
        self._raw_data["Database"] = self._database
        self._raw_data["ReceptorID"] = [str(uuid.uuid4()) for _ in self._raw_data.index]
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


    def clean(self):
        """
        For PIRD use these filters:

        TCR-EPI:
            1) Has reference
            2) receptor - alphabeta and chain is alpha or beta 
            3) epitopes with valid sequence and not null.
            4) cdr3 with valid sequence and not null
            5) receptor, host organism is human or mouse

        TCR-MHC:
            1) MHC first or second class
            3) Chains are known until certain protein or more info

        Nothing to return. Assign values to _tcr_epi and _tcr_mhc
        """
        
        receptor_columns = ["CDR3.alpha.aa","CDR3.beta.aa"]
        pird_cdr3_reshaped = pd.melt(self._raw_data, 
                              id_vars = self._raw_data.columns[~self._raw_data.columns.isin(receptor_columns)],                  
                              value_vars = receptor_columns,
                              var_name = "chain",
                              value_name = "cdr3_seq",
                              ignore_index = True)
        
        pird_cdr3_reshaped.loc[pird_cdr3_reshaped["Species"] == 'Homo Sapiens','Species'] = 'human'
        pird_cdr3_reshaped.loc[pird_cdr3_reshaped["Species"] == 'Mus Musculus','Species'] = 'mouse'

        pird_cdr3_reshaped.loc[pird_cdr3_reshaped["chain"] == 'CDR3.alpha.aa','chain'] = 'alpha'
        pird_cdr3_reshaped.loc[pird_cdr3_reshaped["chain"] == 'CDR3.beta.aa','chain'] = 'beta'
        
        # For union V,D,J columns from different chains
        alpha = pird_cdr3_reshaped.query("chain =='alpha'")
        beta = pird_cdr3_reshaped.query("chain =='beta'")
        alpha1 = alpha.rename(columns = {"Valpha":"V","Jalpha":"J"}).drop(columns = ["Vbeta","Jbeta","Dbeta"])
        alpha1["D"] = pd.NA

        beta1 = alpha.rename(columns = {"Vbeta":"V","Jbeta":"J","Dbeta":"D"}).drop(columns = ["Valpha","Jalpha"])
        pird_raw_ready = pd.concat([alpha1,beta1],ignore_index = True)        

        self._tcr_epi = pird_raw_ready.query("`Pubmed.id`.notna()").\
            query("cdr3_seq.notna()").\
            query("`Antigen.sequence`.notna()").\
            query("chain == 'alpha' or chain == 'beta'").\
            query("Species.isin(@self._included_species)").\
            query("cdr3_seq.str.contains(@self._seq_pattern)").\
            query("`Antigen.sequence`.str.contains(@self._seq_pattern)")
                
        #-------------------------------------------------------------------------------------
        # TCR-EPI STOP, NEXT STEP FOR TCR-MHC                                                # 
        #-------------------------------------------------------------------------------------

        self._tcr_mhc = self._tcr_epi.query("HLA.notna()").query("HLA.str.count(':') >= 1")
        self._tcr_mhc["HLA"] = "HLA-" + self._tcr_mhc["HLA"]

        
        #-------------------------------------------------------------------------------------
        # TCR-MHC STOP, Make table to the one from                                           # 
        #-------------------------------------------------------------------------------------
       
        self._unificate()
    #---------------------------------------------------------------------------------------------------
    def _unificate(self):
        super()._unificate()       
        self._tcr_epi["V"] = self._tcr_epi["V"].replace({"-":pd.NA})
        self._tcr_epi["J"] = self._tcr_epi["J"].replace({"-":pd.NA})
        self._tcr_mhc["V"] = self._tcr_mhc["V"].replace({"-":pd.NA})
        self._tcr_mhc["J"] = self._tcr_mhc["J"].replace({"-":pd.NA})

class MIRACleaner(Cleaner):

    def __init__(self, input, output):
        super().__init__(input, output)
        self._database = "MIRA"
        self._raw_data["Database"] = self._database
        self._raw_data['Chain'] = "beta"
        self._raw_data['Species'] = "human"
        self._raw_data['EpitopeOrganism'] = "SARS-CoV2"
        self._raw_data['D'] = pd.NA
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
                                  "allele":"Activity",
                                  "V":"V",
                                  "D":"D",
                                  "J":"J"}


    def clean(self):
        """
        For MIRA use these filters:

        TCR-EPI:
            1) epitopes with valid sequence and not null.
            2) cdr3 with valid sequence and not null

        TCR-MHC:
            1) Chains are known until certain protein or more info

        Nothing to return. Assign values to _tcr_epi and _tcr_mhc
        """       

        self._tcr_epi = self._raw_data.\
            query("CDR3.str.contains(@self._seq_pattern)").\
            query("epitope.str.contains(@self._seq_pattern)")
                
        #-------------------------------------------------------------------------------------
        # TCR-EPI STOP, NEXT STEP FOR TCR-MHC                                                # 
        #-------------------------------------------------------------------------------------

        self._tcr_mhc = self._tcr_epi.\
        query("allele.notna()").\
        query("allele.str.count(':') >= 1")
        self._tcr_mhc["allele"] = "HLA-" + self._tcr_mhc["allele"]

        
        #-------------------------------------------------------------------------------------
        # TCR-MHC STOP, Make table to the one from                                           # 
        #-------------------------------------------------------------------------------------
       
        self._unificate()
    #---------------------------------------------------------------------------------------------------
    def _unificate(self):
        super()._unificate()       


def get_cleaner(database, input, output):
    if database == 'IEDB' or database == 'CEDAR':
        return IEDBCleaner(input, output, database)
    elif database == 'VDJdb':
        return VDJdbCleaner(input, output)
    elif database == 'McPAS':
        return McPASCleaner(input, output)
    elif database == 'PIRD':
        return PIRDCleaner(input, output)
    elif database == 'MIRA':
        return MIRACleaner(input, output)
    else:
        return ValueError("Unknown database. Correct choices are ['IEDB', 'CEDAR', 'VDJdb','McPAS','PIRD','MIRA']")

@click.command
@click.option("-d","--database","database", required = True, type=click.Choice(['IEDB', 'CEDAR', 'VDJdb','McPAS','PIRD','MIRA']),help = "Database to clean")
@click.option("-i","--input","input", required = True, type=str,help = "Path to input csv file. Expect csv with ';' as separator")
@click.option("-o","--output","output", required = True, type=str,help = "Output csv with ';' as separator")
def main(database, input, output):
    cleaner = get_cleaner(database, input, output)
    cleaner.clean()
    cleaner.save()
    print(f"{database} SUCCESS!")

if __name__ == "__main__":
    main()