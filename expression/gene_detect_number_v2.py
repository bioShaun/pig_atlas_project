from collections import OrderedDict
import pandas as pd
import fire
from pathlib import Path


def detected_genes(overall_df, cutoff):
    overall_cut = overall_df > cutoff
    return overall_cut[overall_cut].index


def check_col(exp_df, colname):
    if colname not in exp_df.columns:
        print(f'{colname} not in expression talbe.')
        return False
    else:
        return True


def exp_number_summary(gene_exp, group_inf, cutoff, outfile,
                       group_file1=None, group_file2=None):

    gene_exp_df = pd.read_table(gene_exp, index_col=0)
    tissue_df = pd.read_table(group_inf, index_col=1, header=None)
    tissue_df.columns = ['group_id']

    tissue_gene_exp_df = pd.merge(
        gene_exp_df.T, tissue_df,
        left_index=True, right_index=True,
        how='left')
    tissue_max_df = tissue_gene_exp_df.groupby(['group_id']).max().T

    exp_number_stats_dict = OrderedDict()

    if group_file1 is not None:
        group_df1 = pd.read_table(group_file1,
                                  header=None,
                                  index_col=0)
        all_groups = group_df1.index
    else:
        all_groups = tissue_max_df.columns

    for each_t in all_groups:
        if not check_col(tissue_max_df, each_t):
            continue
        each_tissue_df = tissue_max_df.loc[:, each_t]
        gene_index = detected_genes(each_tissue_df,
                                    cutoff)
        if group_file2 is None:
            col1, col2 = 'group_id', 'group_expressed_genes'
            exp_number_stats_dict.setdefault(
                col1, []).append(each_t)
            exp_number_stats_dict.setdefault(
                col2, []).append(len(gene_index))
        else:
            col1, col2 = 'group1_id', 'group1_expressed_genes'
            group_df2 = pd.read_table(group_file2,
                                      header=None,
                                      index_col=0)
            for each_t2 in group_df2.index:
                exp_number_stats_dict.setdefault(
                    col1, []).append(each_t)
                exp_number_stats_dict.setdefault(
                    col2, []).append(len(gene_index))
                if not check_col(tissue_max_df, each_t2):
                    continue
                each_tissue_df2 = tissue_max_df.loc[:, each_t2]
                gene_index2 = detected_genes(each_tissue_df2,
                                             cutoff)
                exp_number_stats_dict.setdefault(
                    'group2_id', []).append(each_t2)
                exp_number_stats_dict.setdefault(
                    'group2_expressed_genes', []).append(len(gene_index2))
                exp_number_stats_dict.setdefault(
                    'group1_unique_genes', []).append(
                        len(gene_index.difference(gene_index2))
                )
                exp_number_stats_dict.setdefault(
                    'group2_unique_genes', []).append(
                        len(gene_index2.difference(gene_index))
                )
                exp_number_stats_dict.setdefault(
                    'group1_group2_intersection', []).append(
                        len(gene_index.intersection(gene_index2))
                )
    exp_number_stats_df = pd.DataFrame(exp_number_stats_dict)
    exp_number_stats_df.loc[:, 'cutoff'] = cutoff
    # out_dir = Path(out_dir)
    # if not out_dir.exists():
    #     out_dir.mkdir()
    # exp_number_stats_file = out_dir / 'expressed_gene_number_by_group.txt'
    exp_number_stats_df.to_csv(
        outfile,
        sep='\t',
        index=False)


if __name__ == '__main__':
    fire.Fire(exp_number_summary)
