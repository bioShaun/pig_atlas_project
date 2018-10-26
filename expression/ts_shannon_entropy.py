import pandas as pd
import numpy as np
import click
import os


TS_CUTOFF = 2
UB_CUTOFF = 3.7


def shannon_entropy(series):
    p_norm = series[series > 0] / series.sum()
    return -sum(np.log2(p_norm) * p_norm)


@click.command()
@click.option(
    '-m',
    '--matrix',
    help='expression tpm/fpkm matrix.',
    required=True,
    type=click.Path(dir_okay=False)
)
@click.option(
    '-g',
    '--group',
    help='group information of each sample.',
    required=True,
    type=click.Path(dir_okay=False)
)
@click.option(
    '-c',
    '--gene_classify',
    help='gene classification',
    required=True,
    type=click.Path(dir_okay=False)
)
@click.option(
    '-o',
    '--out_dir',
    help='output directory.',
    default=os.getcwd(),
    type=click.Path(file_okay=False)
)
@click.option(
    '--exp_cut',
    default=0.1,
    type=click.FLOAT,
)
def main(matrix, group, gene_classify, out_dir, exp_cut):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    exp_df = pd.read_table(matrix, index_col=0)
    exp_df = exp_df[exp_df.T.max() > exp_cut]
    group_df = pd.read_table(group, index_col=1, header=None)
    gene_df = pd.read_table(gene_classify, index_col=0)
    group_df.columns = ['tissue']

    group_exp_df = pd.merge(exp_df.T, group_df,
                            left_index=True, right_index=True)
    group_max_df = group_exp_df.groupby('tissue').max().T
    passed_genes = group_max_df[group_max_df.T.max() > exp_cut].index
    group_mean_exp_df = group_exp_df.groupby('tissue').mean().T
    group_mean_exp_df = group_mean_exp_df.loc[passed_genes]
    group_mean_exp_df = group_mean_exp_df.round(3)
    group_mean_exp_df.index.name = 'Gene_id'
    group_mean_exp_file = os.path.join(out_dir, 'tissue.mean.exp.txt')
    group_mean_exp_df.to_csv(group_mean_exp_file, sep='\t')
    log_group_mean_exp_df = np.log2(group_mean_exp_df + 1)
    shannon_entropy_df = log_group_mean_exp_df.T.apply(shannon_entropy)
    max_exp_tissue = group_mean_exp_df.T.idxmax()
    ts_df = pd.concat([shannon_entropy_df, max_exp_tissue], axis=1)
    ts_df.columns = ['shannon_entropy', 'tissue']
    ts_df.index.name = 'Gene_id'
    ts_df.loc[:, 'ts'] = [
        each < TS_CUTOFF for each in ts_df.loc[:, 'shannon_entropy']]
    ts_by_gene_df = pd.merge(ts_df, gene_df,
                             left_index=True, right_index=True, how='left')
    all_out = os.path.join(out_dir, 'shannon_entropy.score.txt')
    ts_out = os.path.join(out_dir, 'tissue_specific.genes.txt')
    ub_out = os.path.join(out_dir, 'ubiquitous.genes.txt')
    ts_genes = ts_by_gene_df[ts_by_gene_df.shannon_entropy < TS_CUTOFF]
    ub_genes = ts_by_gene_df[ts_by_gene_df.shannon_entropy > UB_CUTOFF]
    ts_genes.to_csv(ts_out, sep='\t')
    ub_genes.to_csv(ub_out, sep='\t')
    ts_num_out = os.path.join(out_dir, 'tissue_specific.number.txt')
    ts_num_summary_out = os.path.join(
        out_dir, 'tissue_specific.number.summary.txt')
    ts_gene_num = ts_by_gene_df.groupby(['tissue', 'gene_biotype'])['ts'].sum()
    ts_gene_num_df = ts_gene_num.reset_index()
    ts_gene_num_df.to_csv(ts_num_out, sep='\t', index=False)
    ts_by_gene_df.to_csv(all_out, sep='\t')
    gene_num = ts_by_gene_df.gene_biotype.value_counts()
    ts_num = ts_by_gene_df.groupby('gene_biotype')['ts'].sum()
    ts_over_all_summary = pd.concat([gene_num, ts_num], axis=1)
    ts_over_all_summary.columns = ['detected_genes', 'ts_genes']
    ts_over_all_summary.loc[:, 'ts_portion'] = ts_over_all_summary.ts_genes / \
        ts_over_all_summary.detected_genes
    ts_over_all_summary.index.name = 'Gene_type'
    ts_over_all_summary.to_csv(ts_num_summary_out, sep='\t',
                               float_format='%.3f')
    ts_type_df = ts_by_gene_df.copy()
    ts_type_df.loc[:, 'exp_type'] = 'other'
    ts_type_df.loc[ub_genes.index, 'exp_type'] = 'Ubiquitous'
    ts_type_df.loc[ts_genes.index, 'exp_type'] = 'Tissue specific'
    ts_type_df = ts_type_df.loc[:, 'exp_type']
    ts_type_file = os.path.join(out_dir, 'gene.exp_type.txt')
    ts_type_df.to_csv(ts_type_file, sep='\t')
    # exp_num_dict = dict()
    # for each_tissue in group_mean_exp_df.columns:
    #     exp_num_dict.setdefault('tissue', []).append(each_tissue)
    #     exp_num_dict.setdefault('expressed_number', []).append(
    #         sum(group_mean_exp_df.loc[:, each_tissue] > 0.1)
    #     )
    # exp_num_df = pd.DataFrame(exp_num_dict)


if __name__ == '__main__':
    main()
