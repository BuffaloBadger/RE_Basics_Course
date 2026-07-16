"""Post processing for the Class 28 Practice Assignment from REB, The Course"""

import pandas as pd
import os
from great_tables import GT, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# read the experimental data file
df = pd.read_csv('results.csv')
# drop the last column of units
df = df.drop('units', axis=1)
# rename the third column
df = df.rename(columns = {'lower_limit': '95% Confidence Interval'})
# add units to the first column and format as html
df['Parameter'] = ['k<sub>0</sub> <br> (L mol<sup>-1</sup> min<sup>-1</sup> atm<sup>-2</sup>)'
        , 'E <br> (kJ mol<sup>-1</sup>)', 'α<sub>A</sub>', 'α<sub>B</sub>', 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .fmt_scientific(columns=[1,2,3], rows=[0], decimals=2, exp_style='x10n')
    .cols_merge(columns=[2,3], rows=[0,1,2,3], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[4])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[4], missing_text=' ')
)
results_tbl.show()
results_tbl.gtsave('png/results_table.png', zoom=4.0)