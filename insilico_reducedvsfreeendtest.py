import pandas as pd

def key(df):
    return list(zip(df['Hex'], df['HexNAc'], df.get('NeuAc',0), df.get('NeuGc',0), df.get('KDN',0), df.get('Fuc',0)))

A = pd.read_csv("zf_sPerMeNG_intestine_insilico_20250824.csv.csv").assign(t=lambda d: key(d)).drop_duplicates('t').set_index('t') #reduced
B = pd.read_csv("zf_sPerMeNG_intestine_insilico_20250825.csv.csv").assign(t=lambda d: key(d)).drop_duplicates('t').set_index('t') #non reduced

mass_col = 'Mass'  # or your actual mass column name
both = A[[mass_col]].join(B[[mass_col]], lsuffix='_r', rsuffix='_nr', how='inner')
both['delta_Da'] = both[f'{mass_col}_r'] - both[f'{mass_col}_nr']
print(both['delta_Da'].describe())
print(both['delta_Da'].round(4).value_counts().head())