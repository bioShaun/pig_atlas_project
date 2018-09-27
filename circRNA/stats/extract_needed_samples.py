import fire
import pandas as pd


def extract_samples_from_circ_summary(circ_summary, sample_inf, outfile):
    circ_df = pd.read_table(circ_summary, index_col=0)
    sample_df = pd.read_table(sample_inf, index_col=1,
                              header=None,
                              names=['group_id'])
    out_colunms = list(circ_df.columns[:14])
    out_colunms.extend(list(sample_df.index))
    host_gene_col = [f'{each}(host gene)' for
                     each in sample_df.index]
    out_colunms.extend(host_gene_col)
    out_colunms.extend(list(circ_df.columns[-3:]))
    out_circ_df = circ_df.loc[:, out_colunms]
    out_circ_df.to_csv(outfile, sep='\t')


if __name__ == '__main__':
    fire.Fire(extract_samples_from_circ_summary)
