"""Post processing for Example 11.5.5 from REB, The Book"""

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
table_df = pd.read_csv('example_11_5_5_data.csv')
table_df = table_df.head(6)
data_tbl = (GT(table_df)
    .cols_label(yA_in=html('y<sub>A,in</sub>'),
                yB_in=html('y<sub>B,in</sub>'),
                yY_in=html('y<sub>Y,in</sub>'),
                yZ_in=html('y<sub>Z,in</sub>'),
                T=html('T <br> (°C)'),
                PA_out=html('P<sub>A,out</sub> <br> (atm)')
    )
    .fmt_number(columns=[4], decimals=0)
    .cols_align(
        align='center',
        columns=[0,1,2,3,4,5]
    )
)
data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)