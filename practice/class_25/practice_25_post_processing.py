"""Pre- and Post-Processing for the Class 25 Learning Activity in 
REB, The Course"""

# import libraries
import pandas as pd
import os
from great_tables import GT, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# read the data from the csv file
df = pd.read_csv('practice_25_data.csv')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .tab_stub(rowname_col='Experiment')
    .tab_stubhead(label='Experiment')
    .tab_spanner(
        label='Adjusted Inputs',
        columns = [1,2,3]
    )
    .tab_spanner(
        label= 'Response',
        columns = [4]
    )
    .cols_label(T=html('T (°C)'),
                AtoB=html('Initial A to B'),
                t=html('t<sub>meas</sub><br>(min)'),
                yZ=html('y<sub>Z</sub>')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png')

# create, show, and save a results table
# read the data
df = pd.read_csv('results.csv')
# drop the last column of units
df = df.drop('units', axis=1)
# rename the third column
df = df.rename(columns = {'lower_limit': '95% CI'})
# add units to the first column and format as html
df['Parameter'] = ['k<sub>0</sub> <br> (L mol<sup>-1</sup> min<sup>-1</sup>)'
        , 'E <br> (kJ mol<sup>-1</sup>)', 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .fmt_scientific(columns=[1,2,3], rows=[0], exp_style='x10n')
    .fmt_number(columns=[1,2,3], rows=[1], decimals=1)
    .fmt_number(columns=[1],rows=[2], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
#results_tbl.show()
results_tbl.gtsave('png/results_table.png', zoom=4.0)