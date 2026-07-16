"""Post processing for the Class 25 Learning Activity in REB, The Course"""

import pandas as pd
import os
from great_tables import GT

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# create, show, and save a results table
# read the data
df = pd.read_csv('results/results.csv')
# drop the last column of units
df = df.drop('units', axis=1)
# rename the third column
df = df.rename(columns = {'lower_limit': '95% CI'})
# add units to the first column and format as html
df['Parameter'] = ['V<sub>max</sub> <br> (mmol L<sup>-1</sup> min<sup>-1</sup>)'
        , 'K<sub>m</sub> <br> (mmol L<sup>-1</sup>)', 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
results_tbl.show()
results_tbl.gtsave('png/results_table.png', zoom=4.0)
