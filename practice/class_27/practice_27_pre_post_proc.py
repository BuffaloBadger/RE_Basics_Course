"""Post processing for the Class 27 Practice Assignment in REB, The Course"""

import pandas as pd
import os
from great_tables import GT, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# read the experimental data file
df = pd.read_csv('practice_27_data.csv')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .cols_label(T=html('T<br>(°C)'),
                CA0=html('C<sub>A,0</sub><br>(M)'),
                CB0=html('C<sub>B,0</sub><br>(M)'),
                t_meas=html('t<sub>meas</sub><br>(min)'),
                CA_meas=html('C<sub>A,meas</sub><br>(M)')
    )
    .fmt_number(columns=[0,1], decimals=0)
    .cols_align(
        align='center',
        columns=[0,1,2,3,4,5]
    )
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)

# create, show, and save a model plot results table
# read the data
df = pd.read_csv('model_plot_results.csv')
# build the table
model_tbl = (GT(df)
    .fmt_number(columns=[0], decimals=0)
    .fmt_number(columns=[1,2,3], n_sigfig=3)
    .fmt_number(columns=[4], decimals=3)
    .cols_label(
        T=html('T <br> (°C)'),
        k=html('k <br> (L mol<sup>-1</sup> min<sup>-1</sup>)'),
        k_lower=html('95% Confidence Interval'),
        R_sq=html('R<sup>2</sup>')
    )
    .cols_merge(columns=[2,3], pattern='[{0}, {1}]')
    .cols_align(align='center', columns=[0,1,2,3,4])
)
#model_tbl.show()
model_tbl.gtsave('png/model_plot_results_table.png', zoom=4.0)

# create, show, and save a results table
# read the data
df = pd.read_csv('Arrhenius_results.csv')
# drop the last column of units
df = df.drop('units', axis=1)
# rename the third column
df = df.rename(columns = {'lower_limit': '95% Confidence Interval'})
# add units to the first column and format as html
df['Parameter'] = ['k<sub>0</sub> <br> (L mol<sup>-1</sup> min<sup>-1</sup>)'
        , 'E <br> (kcal mol<sup>-1</sup>)', 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .fmt_scientific(columns=[1,2,3], rows=[0], decimals=2, exp_style='x10n')
    .fmt_number(columns=[1,2,3], rows=[1], decimals=1)
    .fmt_number(columns=[1], rows=[2], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
#results_tbl.show()
results_tbl.gtsave('png/results_table.png', zoom=4.0)