def func():
	import pandas as pd
	df = pd.read_excel('lanz_dataset.xlsx', sheet_name = 'Lanz et al. Dataset')
	print(df[df['Gene(s)'].apply(lambda elem: type(elem) != str)])
