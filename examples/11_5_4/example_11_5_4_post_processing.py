"""Post processing for Example 11.5.4 from REB, The Book"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from great_tables import GT, md, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# table with the first six rows of the experimental data
table_df = pd.read_csv('example_11_5_4_data.csv')
table_df = table_df.head(6)
data_tbl = (GT(table_df)
    .cols_label(T=html('T <br> (°C)'),
                PA0=html('P<sub>A,0</sub> <br> (atm)'),
                PB0=html('P<sub>B,0</sub> <br> (atm)'),
                t_meas=html('t<sub>meas</sub> <br> (min)'),
                fA=html('f<sub>A</sub>')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
)
data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)