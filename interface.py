## analyzes lanz_interface
import pandas as pd
import itertools
def analyze_interface():
    df = pd.read_excel('lanz_interface.xlsx', sheet_name = "INSIDER_Phosphosites")
    df['P1_Interface_residues_numbers'] = df['P1_Interface_residues'].apply(convert_list_to_numbers)
    df['P2_Interface_residues_numbers'] = df['P2_Interface_residues'].apply(convert_list_to_numbers)
    return df

#pd.read_csv()

def foo(df):
    df_split = df.groupby(df.columns[0])
    i = 0
    dic = {}
    for uniprot_id in df_split.groups:
        group = (df_split.get_group(uniprot_id))
        lst = list(set(itertools.chain.from_iterable(group['Interface_residues_numbers'].to_list())))
        dic[uniprot_id] = lst
    return dic


# take loc of P1, P1_interface_phosphosites

# converts the list given in interface spreadsheet to a list of single numbers
def convert_list_to_numbers(lst_str):
    elems = lst_str[1:-1]
    if len(elems) == 0:
        return []
    nums = elems.split(',')
    num_list = []
    for num in nums:
        if '-' in num:
            start_end_lst = num.split('-')
            start = int(start_end_lst[0])
            end = int(start_end_lst[1])
            num_list.extend(list(range(start, end + 1)))
        else:
            num_list.append(int(num))
    return num_list

def get_interface():
    x = analyze_interface()
    print(x.head())
    print(x.tail())
    p1 = x.loc[:, ["P1", "P1_Interface_residues_numbers"]]
    p1.columns = ['Protein', 'Interface_residues_numbers']
    p2 = x.loc[:, ["P2", "P2_Interface_residues_numbers"]]
    p2.columns = ['Protein', 'Interface_residues_numbers']
    y = pd.concat([p1, p2])
    print(y.head())
    print(y.tail())
    return foo(y)

if __name__ == "__main__":
    print(get_interface())
    raise SystemExit
    print(analyze_interface())
    dic = (foo(analyze_interface().loc[:, ["P1", "P1_Interface_residues_numbers"]]))
    print(dic)
    #print(dic[])
    #print(convert_list_to_numbers('[19-21,23-24,50,60,63,67,81,84-85,90,95-97,100,103-104,132,143,146-147,150,160,162-163,167,169-173,175-176,179-180,183,186-187,190-191,195,198-200,204-205,208-209,217,220,233-234,236-237,240,258,262-263,265-266,268-271,277,280-281,283-284]'))