"""Pre- and PostPprocessing for the Class 29 Practice Assignment from REB, The 
Course"""

import pandas as pd
import os
from great_tables import GT, html, loc

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# generate a table showing the first 6 data points.
# read the data
df = pd.read_csv('practice_29_data.csv')
table_df = df.head(6)
data_tbl = (GT(table_df)
    .tab_footnote('Time in minutes, concentration in moles per liter')
    .tab_spanner(label='Experiment 1', columns=[0,1])
    .tab_spanner(label='Experiment 2', columns=[2,3])
    .tab_spanner(label='Experiment 3', columns=[4,5])
    .cols_label(
        t_meas_exp1=html('t<sub>meas</sub>'),
        CY_meas_exp1=html('C<sub>Y,meas</sub>'),
        t_meas_exp2=html('t<sub>meas</sub>'),
        CY_meas_exp2=html('C<sub>Y,meas</sub>'),
        t_meas_exp3=html('t<sub>meas</sub>'),
        CY_meas_exp3=html('C<sub>Y,meas</sub>')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4,5]
    )
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png')

# create, show, and save a model plots results table
# read the data
df = pd.read_csv('model_plot_results.csv')
# build the table
results_tbl = (GT(df)
    .opt_footnote_marks(marks=['*'])
    .tab_footnote(
        footnote=html('°C'),
        locations=loc.column_labels(columns=[0])
    )
    .tab_footnote(
        footnote=html('L mol<sup>-1</sup> min<sup>-1</sup>'),
        locations=loc.column_labels(columns=[1,2])
    )
    .cols_label(k_lower=html('95% Confidence Interval')
        , Rsq=html('R<sup>2</sup>'))
    .fmt_number(columns=[0], decimals=0)
    .fmt_scientific(columns=[1,2,3], decimals=2, exp_style='x10n')
    .fmt_number(columns=[4], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1,2], pattern='[{0}, {1}]')
    .cols_align(align='center', columns=[0,1,2,3,4])
)
results_tbl.show()
results_tbl.gtsave('png/model_results_table.png', zoom=4.0)

# create, show, and save an Arrhenius results table
# read the data
df = pd.read_csv('arrhenius_results.csv')
# drop the last column of units
df = df.drop('units', axis=1)
# rename the third column
df = df.rename(columns = {'lower_limit': '95% CI'})
# add units to the first column and format as html
df['Parameter'] = ['k<sub>0</sub>', 'E', 'R<sup>2</sup>']
# build the table
results_tbl = (GT(df)
    .opt_footnote_marks(marks=['*'])
    .tab_footnote(html('L mol<sup>-1</sup> min<sup>-1</sup>')
        , locations=loc.body(columns=[0], rows=[0]))
    .tab_footnote(html('kcal mol<sup>-1</sup>')
        , locations=loc.body(columns=[0], rows=[1]))
    .fmt_scientific(columns =[1,2,3], rows=[0], decimals=2)
    .fmt_number(columns=[1,2,3], rows=[1,2], n_sigfig=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2])
    .cols_align(align='center', columns=[0,1,2,3])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
results_tbl.show()
results_tbl.gtsave('png/arrhenius_results_table.png', zoom=4.0)


