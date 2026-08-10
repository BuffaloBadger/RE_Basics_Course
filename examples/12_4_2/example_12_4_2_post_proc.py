'''Post processing for Example 12.4.1 from REB, The Book'''

import numpy as np
import pandas as pd
from great_tables import GT, md, html

# results table
table_df = pd.read_csv('results.csv')
data_tbl = (GT(table_df)
    .cols_label(stream=html('Stream<br>Number'),
                conversion=html('Conversion<br>(%)'),
                temperature=html('Temperature<br>(°C)')
    )
    .fmt_number(columns=[1,2], n_sigfig=3)
    .cols_align(
        align='center',
        columns=[0,1,2]
    )
)
data_tbl.show()
data_tbl.gtsave('results_table.png', zoom=4.0)