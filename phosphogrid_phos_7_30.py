import pandas as pd
import os
#from lanz_phos_7_17 import write_phos

def write_phos(phosphosites, uniprot_id, dir = './phosphosites'):
    path = dir + "/" + uniprot_id + "_phos.csv"
    if not os.path.exists(dir):
        os.mkdir(dir)
    if os.path.exists(path):
        os.remove(path)
    phosphosites.to_csv(path, header=False, index=False)
df = pd.read_csv('entrez_uniprot_mapping', sep = '\t', header = None, names=['Entrez', 'UniProt'])

df_commas = df[df.iloc[:,0].str.contains(',')]
ab = []
numpy_commas = df_commas.to_numpy()
for i in numpy_commas:
    for j in i[0].split(','):
        ab.append((j, i[1]))
dfx = pd.DataFrame(ab, columns = ['Entrez', 'UniProt'])

df_nocommas = df[~df.iloc[:,0].str.contains(',')]
df_all = (pd.concat([df_nocommas, dfx]).reset_index(drop = True))

cd = {}
print(df_all)
grouped = (df_all.groupby('Entrez'))
print(grouped.get_group('855816'))

for i in grouped.groups:
    cd[int(i)] = grouped.get_group(i).loc[:, 'UniProt'].to_list()
print(cd)
print(850488 in cd or 854523 in cd)
cd[850488] = ['P43547']
cd[854523] = ['Q12501']
print(len(cd))



df2 = pd.read_csv('phosphogrid_info.txt', sep = '\t')
print(df2)
lst = (list(df2.loc[:,['Entrez Gene ID', 'Position']].to_records(index=False)))
lst2 = ([(cd[elem[0]], elem[1]) for elem in lst])
print(lst2)
print(len(lst2))
print('here')

ef = {}
for i in lst2:
    for j in i[0]:
        if j in ef:
            ef[j].append(i[1])
        else:
            ef[j] = [i[1]]
print(ef)
print(ef['P43547'])

for uniprot in ef:
    write_phos(pd.Series(ef[uniprot]).drop_duplicates(), uniprot, dir = "./phosphosites_phosphogrid")
raise SystemExit
print(len(lst2))
df2['UniProt'] = (df2.loc[:, 'Entrez Gene ID'].apply(lambda elem: cd[elem]))
print(df2.columns)
print(df2)
#print(df2['UniProt'])

