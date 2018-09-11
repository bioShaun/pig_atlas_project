import pandas as pd
import fire


def expression_filter(exp_table, cutoff, outfile=None, method='max'):
    exp_df = pd.read_table(exp_table, index_col=0)
    if method == 'max':
        exp_df = exp_df[exp_df.max(1) > cutoff]
    if outfile is None:
        return exp_df
    else:
        exp_df.to_csv(outfile, sep='\t',
                      float_format='%.5f')


def exp_by_group(exp_table, group_inf, outfile=None, method='mean'):
    exp_df = pd.read_table(exp_table, index_col=0)
    group_df = pd.read_table(group_inf, index_col=1,
                             header=None,
                             names=['group_id'])
    group_exp_df = pd.merge(exp_df.T, group_df,
                            left_index=True,
                            right_index=True)
    group_mean_exp_df = group_exp_df.groupby('group_id').mean().T
    group_mean_exp_df.index.name = 'Gene_id'
    if outfile is None:
        return group_mean_exp_df
    else:
        group_mean_exp_df.to_csv(outfile, sep='\t',
                                 float_format='%.5f')


if __name__ == '__main__':
    fire.Fire()
