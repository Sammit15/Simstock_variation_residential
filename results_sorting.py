import pandas as pd
import numpy as np

df = pd.read_csv('combined_results.csv')

res = ['Cooling', 'Equipment', 'Lights']

for item in res:

    df1 = df.filter(regex=item)

    new_list = list()
    for index, row in df1.transpose().iterrows():
        cleaned = [x for x in list(row) if ~np.isnan(x)]
        new_list.append([index, cleaned[0]])

    df2 = pd.DataFrame(new_list, columns = ['Floor', item])
    df2.to_csv('{}_floors.csv'.format(item))

    df2 = df2.groupby(df2.Floor.str[:13]).sum().reset_index()
    df2.to_csv('{}_polygons.csv'.format(item))
