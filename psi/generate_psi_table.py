import pandas as pd
import fire
from pathlib import Path
import re


def cal_psi(psi_df, sample_id):
    psi_df.loc[:, 'exc_reads'] = psi_df.loc[:, 'exon.count_prev_and_next']
    psi_df.loc[:, 'inc_reads'] = (
        psi_df.loc[:, 'exon.count_prev_and_curr'] +
        psi_df.loc[:, 'exon.count_curr_and_next'])
    psi_df.loc[:, 'sample_id'] = sample_id

    def get_exon_id(chrom, start, end):
        return f'{chrom}_{start}_{end}'

    psi_df.loc[:, 'exon_cur'] = list(map(get_exon_id,
                                         psi_df.loc[:, '#chrom'],
                                         psi_df.loc[:, 'exon.start'],
                                         psi_df.loc[:, 'exon.end']))

    psi_df.loc[:, 'tr_exon_idx'] = psi_df.loc[:, 'transcript_id'] + \
        '_' + psi_df.loc[:, 'exon.index5_3']

    psi_df = psi_df.set_index('tr_exon_idx')

    def get_pre(exon_idx, tr_id):
        nonlocal psi_df
        exon_idx, total = [int(each) for each in exon_idx.split('/')]
        if exon_idx == 1:
            return 'none'
        else:
            pre_idx = f'{tr_id}_{exon_idx-1}/{total}'
            return psi_df.loc[pre_idx, 'exon_cur']

    def get_next(exon_idx, tr_id):
        nonlocal psi_df
        exon_idx, total = [int(each) for each in exon_idx.split('/')]
        if exon_idx == total:
            return 'none'
        else:
            pre_idx = f'{tr_id}_{exon_idx+1}/{total}'
            return psi_df.loc[pre_idx, 'exon_cur']

    psi_df.loc[:, 'exon_pre'] = list(map(get_pre,
                                         psi_df.loc[:, 'exon.index5_3'],
                                         psi_df.loc[:, 'transcript_id']))

    psi_df.loc[:, 'exon_next'] = list(map(get_next,
                                          psi_df.loc[:, 'exon.index5_3'],
                                          psi_df.loc[:, 'transcript_id']))

    psi_df.loc[:, 'exon_idx'] = 'p_' + psi_df.exon_pre + '_' + \
                                'c_' + psi_df.exon_cur + '_' + \
                                'n_' + psi_df.exon_next

    psi_df = psi_df.reset_index()
    psi_df = psi_df[psi_df.exon_pre != 'none']
    psi_df = psi_df[psi_df.exon_next != 'none']

    return psi_df.loc[:, ['sample_id', 'gene_id',
                          'exon_idx', 'inc_reads',
                          'exc_reads']].drop_duplicates()


def load_psi_data(psi_dir, sample_file, out_prefix):
    psi_dir = Path(psi_dir)
    psi_matrix_file = psi_dir / f'{out_prefix}.merged.psi.table.txt'
    if psi_matrix_file.exists():
        psi_matrix_df = pd.read_table(psi_matrix_file)
    else:
        sample_df = pd.read_table(sample_file, header=None,
                                  names=['group_id', 'sample_id'])
        psi_tables = []
        for each in sample_df.sample_id:
            each_table = psi_dir / f'{each}.psi.txt'
            psi_tables.append(each_table)
        psi_dfs = [pd.read_table(each) for each in psi_tables]
        psi_matrix_dfs = list(map(cal_psi, psi_dfs, sample_df.sample_id))
        psi_matrix_df = pd.concat(psi_matrix_dfs)
        psi_matrix_df.to_csv(psi_matrix_file, sep='\t', index=False)
    return psi_matrix_df


def generate_psi_matrix(psi_dir, sample_file, out_prefix,
                        gene_type, reads_cutoff=10, na_cutoff=0.1):
    psi_matrix_df = load_psi_data(psi_dir, sample_file, out_prefix)
    gene_type_df = pd.read_table(gene_type)
    psi_matrix_df = psi_matrix_df.merge(gene_type_df)

    def ret_current_exon(exon_idx):
        cur_exon = re.search('_c_(\S+)_n_', exon_idx).groups()[0]
        return cur_exon

    def extract_psi_matrix(psi_matrix_df, gene_type,
                           psi_dir, method='mean'):
        nonlocal out_prefix
        nonlocal reads_cutoff
        nonlocal na_cutoff
        psi_matrix_df = psi_matrix_df.loc[
            :, ['sample_id', 'exon_idx',
                'inc_reads', 'exc_reads']].drop_duplicates()
        psi_matrix_df = psi_matrix_df.groupby(
            ['sample_id', 'exon_idx']).agg('mean')
        psi_matrix_df = psi_matrix_df.reset_index()
        if method == 'mean':
            psi_matrix_df.loc[:, 'exon_idx'] = psi_matrix_df.exon_idx.map(
                ret_current_exon)
            psi_matrix_df = psi_matrix_df.groupby(
                ['sample_id', 'exon_idx']).agg('sum')
            psi_matrix_df = psi_matrix_df.reset_index()
        psi_matrix_df.loc[:, 'total_reads'] = psi_matrix_df.inc_reads + \
            psi_matrix_df.exc_reads
        psi_matrix_df = psi_matrix_df[
            psi_matrix_df.total_reads >= reads_cutoff]
        psi_matrix_df.loc[:, 'ave_inc_reads'] = psi_matrix_df.inc_reads / 2
        psi_matrix_df.loc[:, 'psi'] = psi_matrix_df.ave_inc_reads / \
            (psi_matrix_df.ave_inc_reads + psi_matrix_df.exc_reads)
        psi_only_df = psi_matrix_df.loc[:, ['sample_id', 'exon_idx', 'psi']]
        psi_only_df = psi_only_df.set_index(['sample_id', 'exon_idx'])
        psi_only_df = psi_only_df.unstack(['sample_id'])
        psi_only_df.columns = psi_only_df.columns.droplevel()
        psi_only_df.index.name = psi_only_df.columns.name
        psi_na_por = psi_only_df.isna().sum(1) / len(psi_only_df.columns)
        psi_na_pass = psi_na_por[psi_na_por < na_cutoff].index
        psi_only_df = psi_only_df.loc[psi_na_pass]
        psi_matrix_out = Path(psi_dir) / f'{gene_type}.{out_prefix}.{method}.cutoff{reads_cutoff}.psi.matrix.txt'
        psi_only_df.to_csv(psi_matrix_out, sep='\t', na_rep='NA',
                           float_format='%.3f')

    for each_type, group in psi_matrix_df.groupby('gene_biotype'):
        extract_psi_matrix(group, each_type, psi_dir)


if __name__ == '__main__':
    fire.Fire(generate_psi_matrix)
