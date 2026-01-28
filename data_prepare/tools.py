from selenium import webdriver

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