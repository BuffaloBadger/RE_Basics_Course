"""Pre- and PostPprocessing for the Class 28 Learning Activity in REB, The 
Course"""

import pandas as pd
import os
from great_tables import GT, html, loc

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# read the experimental data file
df = pd.read_csv('activity_28_data.csv')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .cols_label(yA=html('y<sub>A</sub>'),
                yB=html('y<sub>B</sub>'),
                yY=html('y<sub>Y</sub>'),
                yZ=html('y<sub>Z</sub>'),
                PA=html('P<sub>A,out</sub>')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
    .tab_footnote('atm.',
        locations=loc.column_labels(columns='PA'))
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)

# create, show, and save a results table for guess 1
# read the data
df = pd.read_csv('guess_1/results.csv')
# drop the last column of units
df = df.drop('Units', axis=1)
# rename the third column
df = df.rename(columns = {'lower limit': '95% CI'})
# add units to the first column and format as html
df['Parameter'] = [
    'k'
    , 'K<sub>A</sub>'
    , 'K<sub>B</sub>'
    , 'K<sub>Y</sub>'
    , 'K<sub>Z</sub>'
    , 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .opt_footnote_marks(marks=['*'])
    .tab_footnote(
        footnote=html('mol g<sup>-1</sup> min<sup>-1</sup> atm<sup>-2</sup>'),
        locations=loc.body(columns=[0], rows=[0])
    )
    .tab_footnote(
        footnote=html('atm<sup>-1</sup>'),
        locations=loc.body(columns=[0], rows=[1,2,3,4])
    )
    .fmt_scientific(columns=[1], rows=[0,2,3], decimals=2, exp_style='x10n')
    .cols_merge(columns=[2,3], rows=[0,1,2,3,4], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[5])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[5], missing_text=' ')
    .text_replace(pattern='inf', replacement='∞')
)
#results_tbl.show()
results_tbl.gtsave('guess_1/results_table.png', zoom=4.0)

# repeat for guess 2
df = pd.read_csv('results.csv')
df = df.drop('Units', axis=1)
df = df.rename(columns = {'lower limit': '95% CI'})
df['Parameter'] = [
    'k'
    , 'K<sub>A</sub>'
    , 'K<sub>B</sub>'
    , 'K<sub>Y</sub>'
    , 'K<sub>Z</sub>'
    , 'R<sup>2</sup>']
results_tbl = (GT(df)
    .opt_footnote_marks(marks=['*'])
    .tab_footnote(
        footnote=html('mol g<sup>-1</sup> min<sup>-1</sup> atm<sup>-2</sup>'),
        locations=loc.body(columns=[0], rows=[0])
    )
    .tab_footnote(
        footnote=html('atm<sup>-1</sup>'),
        locations=loc.body(columns=[0], rows=[1,2,3,4])
    )
    .fmt_scientific(columns=[1], rows=[0,2,3], decimals=2, exp_style='x10n')
    .cols_merge(columns=[2,3], rows=[0,1,2,3,4], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[5])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[5], missing_text=' ')
    .text_replace(pattern='inf', replacement='∞')
)
#results_tbl.show()
results_tbl.gtsave('png/results_table.png', zoom=4.0)

# repeat for follow up
df = pd.read_csv('results_b.csv')
df = df.drop('Units', axis=1)
df = df.rename(columns = {'lower limit': '95% CI'})
df['Parameter'] = ['k', 'K<sub>B</sub>', 'R<sup>2</sup>']
results_tbl = (GT(df)
    .opt_footnote_marks(marks=['*'])
    .tab_footnote(
        footnote=html('mol g<sup>-1</sup> min<sup>-1</sup> atm<sup>-2</sup>'),
        locations=loc.body(columns=[0], rows=[0])
    )
    .tab_footnote(
        footnote=html('atm<sup>-1</sup>'),
        locations=loc.body(columns=[0], rows=[1])
    )
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[5])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[3], missing_text=' ')
)
results_tbl.show()
results_tbl.gtsave('png/results_b_table.png', zoom=4.0)

