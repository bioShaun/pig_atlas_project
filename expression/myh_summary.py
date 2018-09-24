import pandas as pd
import fire
import re
from pathlib import Path
import sys


def myh_summary(myh_table, outfile):
    outfile = Path(outfile)
    if outfile.is_file():
        sys.exit('myh_summary finished.')
    myh_df = pd.read_table(myh_table)

    def get_sample_id(sample_name):
        s_pattern = re.compile('.*([1,2,3]).*')
        try:
            sample_id = s_pattern.match(sample_name).groups()[0]
        except AttributeError:
            sample_id = 0
        return sample_id

    myh_df.loc[:, 'sample_id'] = myh_df.gene_name.map(get_sample_id)
    myh_summary = myh_df.groupby(['Part', 'sample_id']).agg(
        {'MYH7_ratio': 'mean', 'MYH2_ratio': 'mean',
         'MYH4_ratio': 'mean', 'MYH7_MYH2': 'mean'})
    myh_summary.to_csv(outfile, sep='\t')


if __name__ == '__main__':
    fire.Fire(myh_summary)
