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


def exp_by_portion(exp_table, group_inf, exp_cutoff,
                   num_cutoff, por_cutoff, outfile=None):
    exp_df = pd.read_table(exp_table, index_col=0)
    exp_pass_df = exp_df > exp_cutoff
    group_df = pd.read_table(group_inf, index_col=1,
                             header=None,
                             names=['group_id'])
    group_exp_df = pd.merge(exp_pass_df.T, group_df,
                            left_index=True,
                            right_index=True)
    exp_pass_df = group_exp_df.groupby('group_id').sum().T
    num_pass_df = exp_pass_df > num_cutoff
    grp_por_df = num_pass_df.sum(1) / num_pass_df.shape[1]
    passed_circs = list(grp_por_df[grp_por_df > por_cutoff].index)
    passed_exp_df = exp_df.loc[passed_circs]
    passed_exp_df.to_csv(outfile, sep='\t', float_format='%.5f')


if __name__ == '__main__':
    fire.Fire()
