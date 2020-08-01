import pandas as pd
import os
from phosphosite_analysis import write_phos

def get_entrez_uniprot_mapping():

    df = pd.read_csv('entrez_uniprot_mapping', sep = '\t', header = None, names=['Entrez', 'UniProt'])

    with_commas = df[df.iloc[:,0].str.contains(',')].to_numpy()
    commas_mapping = []

    for i in with_commas:
        for j in i[0].split(','):
            commas_mapping.append((j, i[1]))
    df_with_commas = pd.DataFrame(commas_mapping, columns = ['Entrez', 'UniProt'])

    df_nocommas = df[~df.iloc[:,0].str.contains(',')]
    df_all = pd.concat([df_nocommas, df_with_commas]).reset_index(drop = True)

    entrez_uniprot_mapping = {} # Each entrez ID mapped to corresponding UniProt ID

    grouped = df_all.groupby('Entrez')

    for i in grouped.groups:
        entrez_uniprot_mapping[int(i)] = grouped.get_group(i).loc[:, 'UniProt'].to_list()
    
    entrez_uniprot_mapping[850488] = ['P43547']
    entrez_uniprot_mapping[854523] = ['Q12501']
    return entrez_uniprot_mapping


entrez_uniprot_mapping = get_entrez_uniprot_mapping()

df2 = pd.read_csv('phosphogrid_info.txt', sep = '\t')

lst = (list(df2.loc[:,['Entrez Gene ID', 'Position']].to_records(index=False)))
lst2 = ([(entrez_uniprot_mapping[elem[0]], elem[1]) for elem in lst])


final_dict = {}
for i in lst2:
    for j in i[0]:
        if j in final_dict:
            final_dict[j].append(i[1])
        else:
            final_dict[j] = [i[1]]

for uniprot in final_dict:
    write_phos(pd.Series(ef[uniprot]).drop_duplicates(), uniprot, dir = "./phosphosites_phosphogrid")

#print(df2['UniProt'])

