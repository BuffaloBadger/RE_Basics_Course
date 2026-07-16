"""Post processing for the Class 26 Practice Assignment in REB, The Course"""

import pandas as pd
import os
from great_tables import GT, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# read the experimental data file
df = pd.read_csv('practice_26_data.csv')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .cols_label(VFR_A_in=html("""
            <span style="position:relative">
            V
            <span style="position:absolute; top:-0.51em; left:0.21em;">
            ·
            </span>
            </span><sub>A,in</sub>
            <br>
            (cm<sup>3</sup> s<sup>-1</sup>)
            """),
        VFR_Y_in=html("""
            <span style="position:relative">
            V
            <span style="position:absolute; top:-0.51em; left:0.21em;">
            ·
            </span>
            </span><sub>Y,in</sub>
            <br>
            (cm<sup>3</sup> s<sup>-1</sup>)
            """),
        VFR_Z_in=html("""
            <span style="position:relative">
            V
            <span style="position:absolute; top:-0.51em; left:0.21em;">
            ·
            </span>
            </span><sub>Z,in</sub>
            <br>
            (cm<sup>3</sup> s<sup>-1</sup>)
            """),
        T=html('T<br>(°C)'),
        Y_to_A=html('α<br>(mol Y / mol A)')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)


# create, show, and save a results table
# read the data
df = pd.read_csv('results.csv')
# drop the last column of units
df = df.drop('units', axis=1)
# rename the third column
df = df.rename(columns = {'lower_limit': '95% CI'})
# add units to the first column and format as html
df['Parameter'] = ['k<sub>0</sub> <br> (L mol<sup>-1</sup> min<sup>-1</sup>)'
        , 'E <br> (kcal mol<sup>-1</sup>)', 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .fmt_scientific(columns=[1,2,3], rows=[0], decimals=2, exp_style='x10n')
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
#results_tbl.show()
results_tbl.gtsave('png/results_table.png', zoom=4.0)

