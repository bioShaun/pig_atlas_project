import fire
import pandas as pd
from pathlib import Path


CFG_TEMP = '''
Annotated_3UTR=/project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/apa/DaPars/analysis_genes.extracted_3UTR.bed
Group1_Tophat_aligned_Wig={case_files}
Group2_Tophat_aligned_Wig={control_files}
Output_directory={out_dir}
Output_result_file={out_name}

# Parameters
Num_least_in_group1=1
Num_least_in_group2=1
Coverage_cutoff=0
FDR_cutoff=0.05
PDUI_cutoff=0.2
Fold_change_cutoff=0.59
'''


def generate_dapars_cfg(case_groups, control_groups,
                        genome_cov_dir, cfg_dir, result_dir,
                        by_group=False):
    case_df = pd.read_table(case_groups, header=None,
                            names=['sample_id'],
                            index_col=0)
    control_df = pd.read_table(control_groups, header=None,
                               names=['sample_id'],
                               index_col=0)
    control_genome_files = [str(Path(genome_cov_dir) / f'{each}.bedgraph')
                            for each in control_df.sample_id]
    control_files = ','.join(control_genome_files)
    case_file_dict = {}
    for each_group in case_df.index:
        each_group_samples = case_df.loc[each_group]
        if each_group_samples.shape == 1:
            each_group_samples_obj = [case_df.loc[each_group].sample_id]
        else:
            each_group_samples_obj = case_df.loc[each_group].sample_id
        if by_group:
            case_file_dict[each_group] = each_group_samples_obj
        else:
            for each_sample in each_group_samples_obj:
                case_file_dict[each_sample] = [each_sample]
    for each_sbj in case_file_dict:
        each_group_samples_obj = case_file_dict[each_sbj]
        case_genome_files = [str(Path(genome_cov_dir) / f'{each}.bedgraph')
                             for each in each_group_samples_obj]
        case_files = ','.join(case_genome_files)
        result_dir = Path(result_dir)
        each_group_result = result_dir / each_sbj
        each_group_result = each_group_result.resolve()
        each_group_cfg = CFG_TEMP.format(case_files=case_files,
                                         control_files=control_files,
                                         out_dir=each_group_result,
                                         out_name=each_sbj)
        cfg_dir = Path(cfg_dir)
        if not cfg_dir.exists():
            cfg_dir.mkdir()
        each_group_cfg_file = cfg_dir / f'{each_sbj}.dapars.cfg'
        with open(each_group_cfg_file, 'w') as cfg_inf:
            cfg_inf.write(each_group_cfg)


if __name__ == '__main__':
    fire.Fire(generate_dapars_cfg)
