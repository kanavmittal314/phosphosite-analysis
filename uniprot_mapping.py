# Maps each UniProt ID to a common gene and systematic gene

import pandas as pd 

# UniProt, Systematic mapping through dictionary
def uniprot_systematic(df):
    return pd.Series(df.loc[:, 'Gene names  (ordered locus )'].values,index=df.loc[:, 'Entry']).fillna(value='').to_dict()

# UniProt, Common mapping through dictionary
def uniprot_common(df):
    return pd.Series(df.loc[:, 'Gene names  (primary )'].values,index=df.loc[:, 'Entry']).fillna(value='').to_dict()

# Reads UniProt table file and gets relevant mappings
def get_mappings(path = './uniprot_table.tab'):
    df = pd.read_csv(path, sep = '\t')
    return uniprot_systematic(df), uniprot_common(df)

def get_all_uniprot_ids(path = './uniprot_table.tab'):
    return pd.read_csv(path, sep = '\t').loc[:, 'Entry']

if __name__ == '__main__':
    print('hello')
    x, y = (get_mappings())
    print(len(x))
    print(len(y))
    print(type(y['P0CD97']))
    