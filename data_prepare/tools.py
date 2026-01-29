from selenium import webdriver
import pandas as pd
import numpy as np
import mhcgnomes
from parallelbar import progress_map
from functools import reduce

def get_chrome_driver(output):
    """
    Get headless Google Chrome selenium driver with output dir as Downloads
    """
    print("Set Chrome options")
    options = webdriver.ChromeOptions()
    options.add_argument('--headless=new')
    options.add_argument('--no-sandbox')
    options.add_argument("--disable-extensions")
    options.add_argument('--remote-debugging-pipe')
    options.add_argument('--disable-dev-shm-usage')

    # Настройки для загрузки файлов
    prefs = {
        "download.default_directory": output,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True
    }
    options.add_experimental_option("prefs", prefs)
    print("Get Chrome driver...")
    driver = webdriver.Chrome(options=options)
    return driver

def get_mhc_class(allele):
    a = mhcgnomes.parse(allele)
    if a.is_class1:
        return "I"
    elif a.is_class2:
        return "II"
    else:
        return "Non-canonical"  


class TCRDatasetStatistics:
    """
    В конструкторе в переменные записывается статистика для дальнейшего перевода и сохранения в JSON. JSON в дальнейшем будет использоваться в дашбоарде.
    """
    
    @staticmethod
    def process_sequence(args):
        """
        функция для параллельного расчета матрицы встречаемости аминокислотных остатков в определенных позициях
        """
        seq, target_max = args
        zeros = np.zeros(shape = (20, target_max))
        matrix = pd.DataFrame(zeros, index = list("ACDEFGHIKLMNPQRSTVWY"), columns = list(range(1, target_max+1)))
        for index, char in enumerate(seq):
            matrix.loc[char, index+1] += 1
        return matrix
    
    def __init__(self, dataset, contain_mhc=False, contain_epitope = True, n_cpu = 10, process_timeout=300):   
        self.num_paired_receptors = self.get_num_paired_receptors(dataset)
        self.sources = dataset["Database"].value_counts(dropna=False)
        self.chains = dataset["Chain"].value_counts(dropna=False)
        self.host_species = dataset["Species"].value_counts(dropna=False)
        self.num_na = dataset.isna().sum()
        self.unique_receptors = len(dataset["ReceptorID"].unique())
        self.unique_structures = len(dataset["Structure"].unique())
        self.unique_activities = len(dataset["Activity"].unique())
        self.host_chain_stat = dataset.loc[:,["Species","Chain"]].value_counts(dropna=False).reset_index()
        self.sample_size_by_activity = dataset["Activity"].value_counts(dropna=False)
        
        self.structure_logo = self.calc_logo(dataset, "Structure", n_cpu, process_timeout)
        
        if contain_epitope:
            self.activity_logo = self.calc_logo(dataset, column = "Activity")
        else:
            self.activity_logo = None
        
        if contain_mhc:
            self.mhc_stat = self.get_stat_by_mhc_class(dataset, column_mhc = "Activity")
        else:
            self.mhc_stat = None

    def calc_logo(self, dataset, column, n_cpu = 10, process_timeout=300):
        target_max = dataset[column].str.len().max()
        sequences = dataset[column].tolist() 
        args = [(s, target_max) for s in sequences]
        mm = progress_map(TCRDatasetStatistics.process_sequence, args, n_cpu = n_cpu, process_timeout=process_timeout, return_failed_tasks = False)
        sum_matrix = reduce(lambda x, y: x + y, mm) # TODO suspicious
        return sum_matrix / len(dataset.index)
    
    def get_stat_by_mhc_class(self, dataset, column_mhc, n_cpu = 10, process_timeout=300):
        classes = progress_map(get_mhc_class, dataset[column_mhc].tolist(), n_cpu = n_cpu, process_timeout=process_timeout, return_failed_tasks = False)
        return pd.DataFrame(classes, columns=['MHC'])["MHC"].value_counts(dropna=False)
    
    def get_num_paired_receptors(self, dataset):
        chains_by_id = dataset[["ReceptorID","Chain"]].value_counts().reset_index()[["ReceptorID","Chain"]]
        counts = chains_by_id["ReceptorID"].value_counts()
        return len(counts[counts == 2])
    

    def save_statistics(self, cleaned_data: pd.DataFrame, output_dir: str) -> None:
            pass