#
# merge sen scores with region files
#

MODEL = "ir" # ir, rs, aad, uni

BASE_PATH = "..."

SEN_PATH = f"{BASE_PATH}/nuc1-norm-{MODEL}-xe24.csv"


CODED_PATH = f"coded-sen-{MODEL}.csv" 
CODED_NUC_PATH = CODED_PATH.replace(".csv", "-nuc.csv")

REGION_FAT_PATH = f"{BASE_PATH}/out-fat.csv"
REGION_EPI_PATH = f"{BASE_PATH}/out-epi.csv"
REGION_TDLU_PATH = f"{BASE_PATH}/out-tdlu.csv"

TEXT_COLS = [0,1] 
DATA_COLS = ['yng', 'sen']
ALL_COLS = ['key', 'y' ] + DATA_COLS    


REGION_TEXT_COLS = [0,1]
REGION_COLS = ['key', 'ishere']


CANCER_META_TEXT_COLS = [0,1,2,5,6,7,8,9,11]
CANCER_META_COLS = ['Barcode','Race','Hispanic','Age','BMI','Donation Year','Biopsy Side','Cancer Side','Self-Report Histology','Cancer Registry Histology','Year Diagnosed','Notes']

TEXT_COLS_MORPH = [0,1]
DATA_COLS_MORPH = ['x', 'y', 'w', 'h', 'min_wh1', 'min_wh2', 'perimeter', 'area', 'hull_perimeter', 'hull_area', 'imean', 'imed', 'i95', 'i99']
ALL_COLS_MORPH = ['key', 'idx'] + DATA_COLS_MORPH  

TISSUE_ORDER = ['epi', 'fat', 'other']


import analyze_utils as autils
import numpy as np
import pandas as pd
import os, sys, re


sys.path.append(os.path.abspath(".."))
from sampler import SampleManager

autils.prep_pub()


data_mp = autils.load_csv(SEN_PATH, TEXT_COLS)
pdf = pd.DataFrame(data_mp, columns=ALL_COLS)
def get_code(key):
    pos = key.find("/")
    key = key[0:pos]
    pos = key.find(" ")
    if pos > 0:
        key = key[0:pos]
    return key

pdf['code'] = pdf['key'].apply(get_code)

def rekey(key):
    key = key.replace(".tif", "")
    key = key.replace(".jpg", "")
    key = re.sub(r'_\d+(_\d+_\d+)', '\\1', key)
    return key

def tissue(row):
    fat = row['isfat']
    epi = row['isepi']
    tdlu = row['istdlu']
    
    if tdlu:
        return "tdlu"
    elif fat and epi:
        return "both"
    elif fat:
        return "fat"
    elif epi:
        return "epi"
    else:
        return "stroma"
    
data_mp = autils.load_csv(REGION_FAT_PATH, header=False, text_cols=REGION_TEXT_COLS, convert_str_None=True)
reg_fat_df = pd.DataFrame(data_mp, columns=REGION_COLS)

data_mp = autils.load_csv(REGION_EPI_PATH, header=False, text_cols=REGION_TEXT_COLS, convert_str_None=True)
reg_epi_df = pd.DataFrame(data_mp, columns=REGION_COLS)

data_mp = autils.load_csv(REGION_TDLU_PATH, header=False, text_cols=REGION_TEXT_COLS, convert_str_None=True)
reg_tdlu_df = pd.DataFrame(data_mp, columns=REGION_COLS)


reg_fat_df['key'] = reg_fat_df['key'].apply(rekey)
reg_epi_df['key'] = reg_epi_df['key'].apply(rekey)
reg_tdlu_df['key'] = reg_tdlu_df['key'].apply(rekey)

reg_fat_df['ishere'] = (reg_fat_df['ishere'] == 'True')
reg_epi_df['ishere'] = (reg_epi_df['ishere'] == 'True')
reg_tdlu_df['ishere'] = (reg_tdlu_df['ishere'] == 'True')


dup_keys = ['K107124/Series 1-0-_0_11x10_270_762', 'K106234/Series 1-0-_0_6x13_504_392','K107710/Series 1-0-_0_27x61_525_547','K108523/Series 1-0-_0_30x56_836_508']

for key in dup_keys:
    reg_fat_df = reg_fat_df.query(f"key!='{key}'")
    reg_epi_df = reg_epi_df.query(f"key!='{key}'")
    reg_tdlu_df = reg_tdlu_df.query(f"key!='{key}'")
    pdf = pdf.query(f"key!='{key}'")

pdf['isfat'] = pdf['key'].map( reg_fat_df.set_index("key")['ishere'])
pdf['isepi'] = pdf['key'].map( reg_epi_df.set_index("key")['ishere'])
pdf['istdlu'] = pdf['key'].map( reg_tdlu_df.set_index("key")['ishere'])
pdf['tissue'] = pdf.apply(tissue, axis=1)

pdfo = pdf.groupby(["code","tissue"]).mean()

pdfo.to_csv(CODED_PATH)

pdf[['code','tissue','key','sen']].to_csv(CODED_NUC_PATH)

