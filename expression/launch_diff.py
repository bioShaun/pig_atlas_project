from __future__ import print_function
import click
import pandas as pd
from pandas import DataFrame
import itertools
import os
import sys
from pathlib import PurePath, Path
from collections import Counter


def check_dir(diff_dir):
    diff_dir = Path(diff_dir)
    if diff_dir.exists() and list(diff_dir.iterdir()):
        return True
    else:
        return False


def get_compare_names(sample_inf, method='pairwise',
                      contrasts=None, replicates=False):
    if contrasts is None:
        group_sample_df = pd.read_table(
            sample_inf, header=None, index_col=0)
        compare_list = list(itertools.combinations(
            group_sample_df.index.unique(), 2))
        tmp_compare_list = compare_list[:]
        if replicates:
            for each_pair in tmp_compare_list:
                shape1 = len(group_sample_df.loc[each_pair[0]].shape)
                shape2 = len(group_sample_df.loc[each_pair[1]].shape)
                if shape1 < 2 or shape2 < 2:
                    compare_list.remove(each_pair)
        if method == 'pairwise':
            compare_name_list = ['{0}_vs_{1}'.format(
                each_compare[0], each_compare[1])
                for each_compare in compare_list]
        elif method == 'preferential':
            unique_group_list = set([item for sublist in compare_list
                                     for item in sublist])
            compare_name_list = ['{g}_vs_other'.format(g=each)
                                 for each in unique_group_list]
        else:
            pass
    else:
        contrasts_df = pd.read_table(
            contrasts, header=None)
        compare_name_list = ['{0}_vs_{1}'.format(contrasts_df.loc[i, 0],
                                                 contrasts_df.loc[i, 1])
                             for i in contrasts_df.index]
    return compare_name_list


def get_file_lines(table, gene_type_df, gene_type, sep='\t'):
    if not os.path.exists(table):
        return []
    table_df = pd.read_table(table, sep=sep, index_col=0)
    merged_table_df = pd.merge(table_df, gene_type_df,
                               left_index=True, right_index=True,
                               how='left')
    selected_table = merged_table_df[merged_table_df.gene_biotype == gene_type]
    return list(selected_table.index)


def cal_gene_type_m_fc(table_df, gene_type_df, gene_type):
    merged_table_df = table_df.merge(gene_type_df,
                                     left_index=True, right_index=True,
                                     how='left')
    selected_table = merged_table_df[merged_table_df.gene_biotype == gene_type]
    return selected_table.logFC.abs().quantile(0.5)


def get_fc_df(table, gene_type_df, gene_type):
    if not os.path.exists(table):
        return None
    df = pd.read_table(table, index_col=0)
    merged_table_df = df.merge(gene_type_df,
                               left_index=True, right_index=True,
                               how='left')
    selected_table = merged_table_df[merged_table_df.gene_biotype == gene_type]
    selected_table.index.name = 'gene_id'
    return selected_table.loc[:, ['gene_biotype', 'logFC']]


def cal_fc_stats(table, gene_type_df, gene_type):
    if not os.path.exists(table):
        return 0
    table_df = pd.read_table(table, index_col=0)
    return cal_gene_type_m_fc(table_df, gene_type_df, gene_type)


def cal_all_fc(table, gene_type_df, gene_type):
    if not os.path.exists(table):
        sys.exit('Diff table not found. [{df}]'.format(df=table))
    table_df = pd.read_table(table, index_col=0)
    up_df = table_df[table_df.logFC > 0]
    down_df = table_df[table_df.logFC < 0]
    up_fc = cal_gene_type_m_fc(up_df, gene_type_df, gene_type)
    down_fc = cal_gene_type_m_fc(down_df, gene_type_df, gene_type)
    return up_fc, down_fc


@click.command()
@click.option(
    '-s',
    '--sample_inf',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
    help='sample information with sample names and group names.'
)
@click.option(
    '-c',
    '--counts',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
    help='gene counts table.'
)
@click.option(
    '-t',
    '--tpm_table',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
    help='table with tpm value.'
)
@click.option(
    '-o',
    '--out_dir',
    type=click.Path(file_okay=False),
    required=True,
    help='output directory.'
)
@click.option(
    '-q',
    '--qvalue',
    type=click.FLOAT,
    default=0.05
)
@click.option(
    '-f',
    '--fc',
    type=click.FLOAT,
    default=1
)
@click.option(
    '-c',
    '--gene_class',
    type=click.Path(dir_okay=False, exists=True),
    required=True,
    help='gene classify table.'
)
@click.option(
    '-m',
    '--method',
    type=click.Choice(['pairwise', 'preferential']),
    default='pairwise',
    help='diff analysisi method.'
)
@click.option(
    '--lib_size',
    type=click.Path(dir_okay=False),
    default='',
    help='user provided sample library size.'
)
def main(sample_inf, counts, tpm_table,
         out_dir, qvalue, fc, gene_class,
         method, lib_size):
    out_dir = PurePath(out_dir)
    compare_list = get_compare_names(sample_inf, method=method)
    for each_compare in compare_list:
        each_compare_out = os.path.join(out_dir, each_compare)
        if not check_dir(each_compare_out):
            sys.exit('{cmp} result not exists!'.format(cmp=each_compare))

    sample_df = pd.read_table(sample_inf, index_col=1, header=None)
    sample_df.columns = ['tissue']
    tissues = sample_df.tissue.unique()
    diff_num_df = pd.DataFrame([[None for each in tissues]
                                for each in tissues])
    diff_num_df.columns = tissues
    diff_num_df.index = tissues

    gene_type_df = pd.read_table(gene_class, index_col=0)
    for each_type in gene_type_df.loc[:, 'gene_biotype'].unique():
        each_diff_num_df = diff_num_df.copy()
        each_fc_val_df = diff_num_df.copy()
        each_all_fc_val_df = diff_num_df.copy()
        each_diff_genes_dict = dict()
        each_diff_genes_app_dict = Counter()
        fc_dfs = []
        diff_fc_dfs = []
        for each_compare in compare_list:
            each_compare_out = os.path.join(out_dir, each_compare)
            each_compare_pair = each_compare.split('_vs_')
            each_up_gene_file = os.path.join(each_compare_out,
                                             '{c}.{p}-UP.edgeR.DE_results.txt'.format(
                                                 c=each_compare,
                                                 p=each_compare_pair[0]
                                             ))
            each_down_gene_file = os.path.join(each_compare_out,
                                               '{c}.{p}-UP.edgeR.DE_results.txt'.format(
                                                   c=each_compare,
                                                   p=each_compare_pair[1]
                                               ))
            all_de_file = os.path.join(each_compare_out,
                                       '{c}.edgeR.DE_results.txt'.format(c=each_compare))
            fc_dfs.append(get_fc_df(all_de_file, gene_type_df, each_type))
            diff_fc_dfs.append(
                get_fc_df(each_up_gene_file, gene_type_df, each_type))
            diff_fc_dfs.append(
                get_fc_df(each_down_gene_file, gene_type_df, each_type))
            each_up_genes = get_file_lines(
                each_up_gene_file, gene_type_df, each_type)
            each_up_fc = cal_fc_stats(
                each_up_gene_file, gene_type_df, each_type)
            each_down_genes = get_file_lines(
                each_down_gene_file, gene_type_df, each_type)
            each_down_fc = cal_fc_stats(
                each_down_gene_file, gene_type_df, each_type)
            each_all_up_fc, each_all_down_fc = cal_all_fc(all_de_file,
                                                          gene_type_df,
                                                          each_type)
            each_diff_num_df.loc[
                each_compare_pair[1], each_compare_pair[0]] = len(
                    each_up_genes)
            each_diff_num_df.loc[
                each_compare_pair[0], each_compare_pair[1]] = len(
                    each_down_genes)
            each_fc_val_df.loc[
                each_compare_pair[1], each_compare_pair[0]
            ] = each_up_fc
            each_fc_val_df.loc[
                each_compare_pair[0], each_compare_pair[1]
            ] = each_down_fc
            each_all_fc_val_df.loc[
                each_compare_pair[1], each_compare_pair[0]
            ] = each_all_up_fc
            each_all_fc_val_df.loc[
                each_compare_pair[0], each_compare_pair[1]
            ] = each_all_down_fc
            all_diff_genes = each_up_genes + each_down_genes
            each_diff_genes_dict.setdefault(
                each_compare_pair[0], []).extend(all_diff_genes)
            each_diff_genes_dict.setdefault(
                each_compare_pair[1], []).extend(all_diff_genes)
            each_diff_genes_app_dict.update(all_diff_genes)

        def merge_fc_df(fc_dfs):
            fc_dfs = [each for each in fc_dfs
                      if each is not None]
            fc_df = pd.concat(fc_dfs)
            fc_df.loc[:, 'logFC'] = fc_df.logFC.abs()
            fc_df = fc_df.groupby(['gene_id', 'gene_biotype'])[
                'logFC'].median()
            return DataFrame(fc_df)

        fc_df = merge_fc_df(fc_dfs)
        fc_file = out_dir / '{t}.fc.median.table.txt'.format(t=each_type)
        fc_df.to_csv(fc_file, sep='\t')
        diff_fc_df = merge_fc_df(diff_fc_dfs)
        diff_fc_file = out_dir / \
            '{t}.diff.fc.median.table.txt'.format(t=each_type)
        diff_fc_df.to_csv(diff_fc_file, sep='\t')

        each_diff_num_df = each_diff_num_df.dropna(how='all')
        each_diff_num_df = each_diff_num_df.loc[:, each_diff_num_df.index]
        for each_index in each_diff_num_df.index:
            each_diff_num_df.loc[each_index, each_index] = 0
        diff_matrix_file = out_dir / '{t}.diff.matrix.txt'.format(t=each_type)
        each_diff_num_df.to_csv(diff_matrix_file, sep='\t')
        each_diff_summary_dict = dict()
        for each_tissue in each_diff_genes_dict:
            each_t_diff_num = len(set(each_diff_genes_dict[each_tissue]))
            each_diff_summary_dict.setdefault(
                'tissue', []).append(each_tissue)
            each_diff_summary_dict.setdefault(
                'diff_num', []).append(each_t_diff_num)
        each_diff_summary_df = pd.DataFrame(each_diff_summary_dict)
        each_diff_summary_df.loc[:, 'gene_biotype'] = each_type
        each_diff_summary_file = out_dir / \
            '{t}.diff.num.txt'.format(t=each_type)
        each_diff_summary_df.to_csv(each_diff_summary_file, sep='\t',
                                    index=False)
        each_df_fc_file = out_dir / '{t}.diff.fc.txt'.format(t=each_type)
        each_fc_val_df.to_csv(each_df_fc_file, sep='\t')
        each_all_fc_file = out_dir / '{t}.all.fc.txt'.format(t=each_type)
        each_all_fc_val_df.to_csv(each_all_fc_file, sep='\t')
        each_diff_genes_app_file = out_dir / \
            '{t}.diff.gene.app.txt'.format(t=each_type)
        each_diff_genes_app_df = pd.DataFrame(
            list(each_diff_genes_app_dict.items()),
            columns=['Gene_ID',
                     'Diff_Pair_Number'])
        each_diff_genes_app_df.to_csv(each_diff_genes_app_file, sep='\t',
                                      index=False)


if __name__ == '__main__':
    main()
