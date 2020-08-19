import pandas as pd
import os
from phosphosite_analysis import write_df, write_phos

def dataset_in(file_name = "lanz_dataset.xlsx", sheet_name = "Lanz et al. Dataset"):
    print("Reading Excel...")
    try:
        excel = pd.read_excel(io = file_name, sheet_name = sheet_name)
        print("Excel Read")
    except:
        print("Excel file problem")
        raise SystemExit  
    return excel

# Not needed anymore
def multiple_proteins(df):
    excluded_phosphosites = (df[df['Phosphosite'].str.contains(';')]).loc[:, ['Phosphosite', 'Uniprot_id'] ]
    print(excluded_phosphosites)
    write_df(excluded_phosphosites, "excluded_phosphosites.csv", dir = "./")
    return df[~df['Phosphosite'].str.contains(';')]

def common_uniprot_lanz(df):
    #df.loc[:, ['Gene(s)_split', 'Uniprot_id_split']] = 
    df['Gene(s)_split'] = df['Gene(s)'].str.split(';')
    df['Uniprot_id_split'] = df['Uniprot_id'].str.split(';')

    df = df.loc[:, ['Gene(s)_split', 'Uniprot_id_split']].apply(pd.Series.explode).reset_index(drop=True)
    return pd.Series(df.loc[:, 'Uniprot_id_split'].values,index=df.loc[:, 'Gene(s)_split']).to_dict()

def main():
    df = dataset_in()
    common_uniprot_lanz_dic = common_uniprot_lanz(df)
    print(common_uniprot_lanz_dic)
    x = (df['Phosphosite'].str.split(';').explode().reset_index(drop=True).str.rsplit('_', n = 1, expand = True))
    print(x)
    x.columns = ['UniProt', 'Residue Number']
    x['UniProt'] = x['UniProt'].replace(common_uniprot_lanz_dic)
    print(x)
    x.to_csv('lanz_dataset_analyzed.csv')
    

    
    #df['SiteLoc'] = df['Phosphosite'].str.split
    df_split = x.groupby(['UniProt'])

    for uniprot_id in df_split.groups.keys():
        
        phos = df_split.get_group(uniprot_id)['Residue Number'].unique()
        write_phos(phos, uniprot_id)
    
    #write_df(df.loc[:, ['Uniprot_id', 'Gene(s)', 'SystematicGeneName']].drop_duplicates(), 'Lanz_systematic_common_gene.csv', dir = '.')

   # Systematic & common gene name write_df(x.loc[:, ['Uniprot_id', '']])
    
if __name__ == '__main__':
    main()
    
'''
dir = "./phosphosite_analysis_results"
    path = dir + "/" + file_name
    if not os.path.exists(dir):
        os.mkdir(dir)
    if os.path.exists(path):
        os.remove(path)
    df.to_csv(path, index=False)
'''
