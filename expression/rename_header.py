import fire
import pandas as pd


def rename_header(exp_table, name_map, rename_table):
    exp_df = pd.read_table(exp_table, index_col=0)
    sample_map = pd.read_table(name_map, header=None,
                               names=['new_name'],
                               index_col=0)
    rename_exp_df = exp_df.rename(
        columns=sample_map.loc[exp_df.columns].new_name)
    rename_exp_df.to_csv(rename_table, sep='\t')


if __name__ == '__main__':
    fire.Fire()
