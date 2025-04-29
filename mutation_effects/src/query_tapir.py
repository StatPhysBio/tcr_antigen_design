
'''
Generated with a lot of help from ChatGPT
'''

import numpy as np
import pandas as pd

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.keys import Keys
import time

# ---------- Configuration ----------
LOGIN_URL = "https://vcreate.io/tapir/signin"  # Update with actual login page URL
FORM_URL = "https://vcreate.io/tapir?page=18"    # Update with actual form URL

# Replace with your actual credentials
USERNAME = ""
PASSWORD = ""


data_batch = []

for tcr in range(5, 8):

    if tcr != 7:
        mhc = "HLA-A*02:01"
    else:
        mhc = "HLA-B*27:05"

    antigens = pd.read_csv(f'../mskcc/mskcc_tcr{tcr}_ec50_sat_mut_af3.csv')['sequence'].values
    for antigen in antigens:

        data_batch.append({
            "file_path": f"/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/mskcc/mskcc_tcr{tcr}_tapir.tsv",
            "file_format": "Vcreate",
            "antigen": antigen,
            "mhc": mhc,
            "model": "Tapir"
        })


# ---------- Start Browser ----------
driver = webdriver.Chrome()
wait = WebDriverWait(driver, 10)

# ---------- Log In ----------
print('logging in...')
driver.get(LOGIN_URL)

# Fill in login form (update selectors as needed)
wait.until(EC.presence_of_element_located((By.XPATH, '//input[@type="email"]'))).send_keys(USERNAME)
driver.find_element(By.XPATH, '//input[@type="password"]').send_keys(PASSWORD)
driver.find_element(By.XPATH, '//input[@type="submit"]').click()

print('logged in!')

# Wait for redirect / confirmation
time.sleep(3)


# ---------- Form Submission Loop ----------
for entry in data_batch:
    # driver.get(FORM_URL)

    # Upload file
    wait.until(EC.presence_of_element_located((By.XPATH, '//input[@type="file"]'))).send_keys(entry["file_path"])

    # Select file format
    file_format_select = Select(driver.find_element(By.XPATH, '//select[option[contains(text(), "Select a file type")]]'))
    file_format_select.select_by_value(entry["file_format"])

    # Antigen
    antigen_input = driver.find_element(By.XPATH, '//input[@type="text"]').send_keys(entry["antigen"])

    # MHC
    # Click the React Select input area
    mhc_input = wait.until(EC.element_to_be_clickable((By.XPATH, '//div[@class="react-select-virtualized css-b62m3t-container"]')))
    mhc_input.click()

    # Find the input field (appears only after clicking)
    search_box = wait.until(EC.presence_of_element_located((By.XPATH, '//div[contains(@class, "react-select")]//input[@type="text"]')))

    # Type the allele and press Enter
    search_box.send_keys(entry["mhc"])
    search_box.send_keys(Keys.ENTER)

    # Model
    model_select = Select(driver.find_element(By.XPATH, '//select[option[contains(text(), "Tapir")]]'))
    model_select.select_by_value(entry["model"])

    # Submit
    driver.find_element(By.XPATH, '//input[@type="submit"]').click()
    time.sleep(4)

driver.quit()

