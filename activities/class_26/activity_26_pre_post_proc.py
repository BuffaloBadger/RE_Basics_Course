"""Post processing for the Class 26 Learning Activity in REB, The Course"""

import pandas as pd
import os
from great_tables import GT, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# read the experimental data file
df = pd.read_csv('activity_26_data.csv')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .cols_label(T=html('T<br>(°C)'),
                CAin=html('C<sub>A,in</sub><br>(M)'),
                CBin=html('C<sub>B,in</sub><br>(M)'),
                tau=html('τ<br>(min)'),
                CA_meas=html('C<sub>A,meas</sub><br>(M)')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)


# create, show, and save a results table for the fitting function approach
# read the data
df = pd.read_csv('results/fitting_function_results.csv')
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

# create, show, and save a results table for the model plots
#read the data
df = pd.read_csv('results/model_plot_results.csv')
# build the table
model_tbl = (GT(df)
    .cols_label(T=html('T (°C)'),
            k=html('k (L mol<sup>-1</sup> min<sup>-1</sup>)'),
            R_squared=html('R<sup>2</sup>'))
    .cols_merge(columns=[1,2,3], pattern='{0}, 95% CI [{0}, {1}]')
    .fmt_number(columns=[1,2,3,4], n_sigfig=3)
    .cols_align(align='center')
)
#model_tbl.show()
model_tbl.gtsave('png/model_plots_results_table.png', zoom=4.0)

# create, show, and save a results table for the linear model approach
# read the data
df = pd.read_csv('results/linear_model_results.csv')
# add units to the first column and format as html
df['Parameter'] = ['k<sub>0</sub> <br> (L mol<sup>-1</sup> min<sup>-1</sup>)'
        , 'E <br> (kcal mol<sup>-1</sup>)', 'R<sup>2</sup>']
# build the table
model_approach_tbl = (GT(df)
    .fmt_scientific(columns=[1,2,3], rows=[0], decimals=2, exp_style='x10n')
    .cols_hide(columns=[2,3,4])
    .cols_align(align='center', columns=[0,1])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
#model_approach_tbl.show()
model_approach_tbl.gtsave('png/linear_results_table.png', zoom=4.0)