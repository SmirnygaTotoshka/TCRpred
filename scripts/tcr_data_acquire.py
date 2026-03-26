#!/usr/bin/env python

import os
import time
import uuid
import click
import zipfile
import requests
import warnings
import traceback
import pandas as pd

from tqdm import tqdm
from io import StringIO
from shutil import rmtree
from datetime import datetime
from sqlalchemy import create_engine, URL
from sqlalchemy.exc import OperationalError

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException, WebDriverException


warnings.filterwarnings("ignore")
            
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

output_option = click.option("-o","--output","output", required = True, type=str, help = "Output directory")

@click.group()
def main():
    pass


@main.command()
@click.option("-s","--server","server", required = True, type=str,help = "Server IP address with local IEDB databases. Should be nextflow secret")
@click.option("-u","--user","user", required = True, type=str,help = "Database user with SELECT privileges in the selected database. Should be nextflow secret")
@click.option("-p","--password","password", required = True, type=str,help = "Database user password. MUST be nextflow secret")
@click.option("-d","--database","database", required = True, type=click.Choice(['IEDB', 'CEDAR'], case_sensitive=True))
@click.option("-c","--chunk_size","chunk_size", required = False, type=int, default = 100000)
@output_option
def mysql(server, user, password, database, output, chunk_size):
    '''
    Download tables from IEDB and CEDAR connected with TCR and bind all into one raw big table by primary keys.
    The merge scheme see experiments/3_download_iedb_cedar.ipynb
    WARNING! Large tables store in RAM. Execute only cluster with slurm and huge RAM
    '''
    date = datetime.today().strftime('%Y-%m-%d')
    try:
        table_names = ["tcell","curated_epitope","object","epitope_object","epitope", "tcell_receptor","curated_receptor","distinct_receptor","distinct_chain","mhc_allele_restriction"]
        tables = {}
        assert chunk_size >= 100 and chunk_size <= 100000000, "Chunk size should be between 100 and 100000000" 
        #DOWNLOAD PART
        url = URL.create(
            "mysql+mysqlconnector",
            username=user,
            password=password, 
            host=server,
            port = 3306,
            database=database
        )
        conn = create_engine(url)
        for name in tqdm(table_names):
            tbl = pd.DataFrame()
            for chunk in tqdm(pd.read_sql(f"SELECT * FROM {name};", con = conn.connect(), chunksize = chunk_size)):
                tbl = pd.concat([tbl, chunk], ignore_index = False)
            tables[name] = tbl
        tables["organism_names"] = pd.read_sql("SELECT organism_id, name_txt FROM organism_names WHERE name_class = 'scientific name';", con = conn.connect())
            
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
        print("Saving...")    
        ready.to_csv(os.path.join(output,f"{database}.csv"),sep = ";",index = False)
    except (OperationalError, pd.errors.DatabaseError) as e:
        traceback.print_exc()
        print(f"Something with MySQL server connection {e}")
    else:
        with open(os.path.join(output, "README"),"w") as readme:
            readme.write(f'''
                The database {database} has been downloaded from the server {server} at {date}.
                It was uploaded to the server 26.03.21 (yy-mm-dd).
                IEDB Last Update: March 08, 2026
                CEDAR Last Update: March 18, 2026
            ''')
    finally:
        print("Finished")

@main.command()
@output_option
def vdjdb(output):
    '''
    Download VDJdb using their server API
    '''
    last_update = ""
    date = datetime.today().strftime('%Y-%m-%d')
    driver = None
    try:
        # Текст, содержащий информацию о последнем обновлении БД
        vdjdb_url = 'https://vdjdb.com/overview'
        vdjdb_last_update_path = "/html/body/application/div/overview/div/div/div/div/pre/code"
        driver = get_chrome_driver(output)
        print(f"Driver = {driver}")
        driver.get(vdjdb_url)
        last_update = WebDriverWait(driver, 40).until(
            EC.presence_of_element_located((By.XPATH, vdjdb_last_update_path))
        ).text
        print(last_update)

        data_url = "https://vdjdb.com/api/database/search"
        meta_url = "https://vdjdb.com/api/database/meta"
        header = {"Content-Type": "application/json"}
        data = {"filters":[] }

        print("Get data...")
        data_result = requests.post(data_url, json = data, headers = header)
        print("Get meta...")
        meta_result = requests.get(meta_url)

        if data_result.status_code == requests.codes['ok'] and meta_result.status_code == requests.codes['ok']:
            print("Convert...")
            data_json = data_result.json()
            meta_json = meta_result.json()
            table = [d['entries'] for d in data_json["rows"]]
            columns = [c['title'] for c in meta_json['metadata']['columns']]
            vdjdb = pd.DataFrame(table, columns = columns)
            vdjdb["receptor_id"] = [d['metadata']['pairedID'] for d in data_json["rows"]]
            vdjdb.to_csv(os.path.join(output,f"VDJdb.csv"),sep = ";", index = False)
        else:
            if data_result.status_code != requests.codes['ok']:
                data_result.raise_for_status()
            else:
                meta_result.raise_for_status()
        
    except requests.HTTPError as e:
        traceback.print_exc()
        print(f"Something with HTTP connection {e}")
    except (TimeoutException, NoSuchElementException, WebDriverException) as e:
        traceback.print_exc()
        print(f"Something with Selenium {e}")
    else:
        with open(os.path.join(output, "README"),"w") as readme:
            readme.write(f'''
                The database VDJdb has been downloaded from the server at {date}.
                {last_update}.
            ''')
    finally:
        if driver is not None:
            driver.quit()
        print("Finished")

@main.command()
@output_option
def mcpas(output):
    '''
    Download McPAS-TCR using selenium web-scrabbing and requests library
    '''
    last_update = ""
    date = datetime.today().strftime('%Y-%m-%d')
    driver = None
    try:
        mcpas_url = 'https://friedmanlab.weizmann.ac.il/McPAS-TCR/'
        download_button = "/html/body/div[1]/div/section/div/div/div[1]/div/div/div/div/div/div[2]/a"
        # Текст, содержащий информацию о последнем обновлении БД
        mcpas_last_update_path = "/html/body/div[1]/div/section/div/div/div[1]/div/div/div/div/div/div[2]/div[3]/p"

        driver = get_chrome_driver(output)
        print(f"Driver = {driver}")
        driver.get(mcpas_url)
        last_update = WebDriverWait(driver, 40).until(
            EC.presence_of_element_located((By.XPATH, mcpas_last_update_path))
        ).text
        print(last_update)
        time.sleep(10) # Нужно дать время, чтобы сервер сгенерировал сессионную ссылку
        link = driver.find_element(By.XPATH,download_button).get_attribute("href")
        print("Download...")
        data_result = requests.get(link)

        if data_result.status_code == requests.codes['ok']:
            print("Convert...")
            result = StringIO(data_result.text)
            mcpas = pd.read_csv(result, sep = ",")
            mcpas.to_csv(os.path.join(output,f"McPAS.csv"),sep = ";", index = False)
        else:
            data_result.raise_for_status()
        
    except requests.HTTPError as e:
        traceback.print_exc()
        print(f"Something with HTTP connection {e}")
    except (TimeoutException, NoSuchElementException, WebDriverException) as e:
        traceback.print_exc()
        print(f"Something with Selenium {e}")
    else:
        with open(os.path.join(output, "README"),"w") as readme:
            readme.write(f'''
                The database McPAS-TCR has been downloaded from the server at {date}.
                {last_update}.
            ''')
    finally:
        if driver is not None:
            driver.quit()
        print("Finished")

@main.command()
@output_option
def pird(output):
    '''
    Download PIRD using selenium web-scrabbing and requests library
    '''
    date = datetime.today().strftime('%Y-%m-%d')
    driver = None
    try:
        
        tbadb_url = "https://ftp.cngb.org/pub/SciRAID/PIRD/TBAdb/TBAdb.xlsx" # TODO это прямая ссылка. Селениум здесь не нужен. 
        file_name = os.path.join(output, f"TBAdb.xlsx")
        # Send a GET request to download the file
        response = requests.get(tbadb_url)
        print("Download...")

        # Check if the request was successful
        if response.status_code == 200:
            # Write the content to a file
            with open(file_name, "wb") as file:
                file.write(response.content)
            print(f"Downloaded {file_name}")
        else:
            raise ValueError(f"Failed to download file: {response.status_code} - {response.text}")

        time.sleep(30)  # Увеличьте время ожидания при необходимости
        print("Convert...")
        tbadb = pd.read_excel(os.path.join(output, "TBAdb.xlsx"), sheet_name = "TCR-AB")
        tbadb.to_csv(os.path.join(output,f"PIRD.csv"),sep = ";", index = False)
        
    except (TimeoutException, NoSuchElementException, WebDriverException) as e:
        traceback.print_exc()
        print(f"Something with Selenium {e}")
    else:
        with open(os.path.join(output, "README"),"w") as readme:
            readme.write(f'''
                The database PIRD has been downloaded from the server at {date}.
                The database doesn`t contain any field about last update information.
            ''')
    finally:
        if driver is not None:
            driver.quit()
        print("Finished")

if __name__ == "__main__":
    main()