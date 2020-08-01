import pandas as pd
import os
from phosphosite_analysis_7_16 import write_df

def write_phos(phosphosites, uniprot_id, dir = './phosphosites'):
    path = dir + "/" + uniprot_id + "_phos.csv"
    if not os.path.exists(dir):
        os.mkdir(dir)
    if os.path.exists(path):
        os.remove(path)
    phosphosites.to_csv(path, header=False, index=False)

def dataset_in(file_name= "Lanz.xlsx", sheet_name = "Lanz et al. Dataset"):
    print("Reading Excel...")
    try:
        excel = pd.read_excel("Lanz.xlsx", sheet_name = sheet_name).loc[:37245]
        print("Excel Read")
    except:
        print("Excel file problem")
        raise SystemExit  
    return excel

def multiple_proteins(df):
    excluded_phosphosites = (df[df['Phosphosite'].str.contains(';')]).loc[:, ['Phosphosite', 'Uniprot_id'] ]
    print(excluded_phosphosites)
    write_df(excluded_phosphosites, "excluded_phosphosites.csv", dir = "./")
    return df[~df['Phosphosite'].str.contains(';')]

def solve():
    df = dataset_in()
    df = multiple_proteins(df)

    df['SiteLoc'] = df.apply(lambda entry: entry['Phosphosite'].split("_")[-1], axis=1)

    df_split = df.groupby(['Uniprot_id'])

    for i in df_split.groups.keys():
        
        phos = df_split.get_group(i)['SiteLoc']
        write_phos(phos, i)
    
    write_df(df.loc[:, ['Uniprot_id', 'Gene(s)', 'SystematicGeneName']].drop_duplicates(), 'Lanz_systematic_common_gene.csv', dir = '.')

   # Systematic & common gene name write_df(x.loc[:, ['Uniprot_id', '']])
    
solve()
    
'''
dir = "./phosphosite_analysis_results"
    path = dir + "/" + file_name
    if not os.path.exists(dir):
        os.mkdir(dir)
    if os.path.exists(path):
        os.remove(path)
    df.to_csv(path, index=False)
'''
