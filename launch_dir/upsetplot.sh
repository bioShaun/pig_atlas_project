#Fri Oct 26 08:01:09 EDT 2018

alias vpy='pipenv run python'

# adipose expression upset plot data
adipose_dir='/project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/diff/adipose'
vpy ../expression/gene_detect_number_by_mean.py \
       ${adipose_dir}/rmpad_exp_table/exp.miRNA.tpm.txt \
       ${adipose_dir}/../muscle/cfg_dir/miRNA.type.txt \
       ${adipose_dir}/cfg_dir/rmpad.reorder.adipose.sample.ini \
       1 \
       ${adipose_dir}/exp_stats/exp_number/miRNA

vpy ../expression/gene_detect_number_by_mean.py \
       ${adipose_dir}/rmpad_exp_table/exp.circRNA.tpm.txt \
       ${adipose_dir}/../muscle/cfg_dir/circRNA.type.txt \
       ${adipose_dir}/cfg_dir/rmpad.reorder.adipose.sample.ini \
       0.05 \
       ${adipose_dir}/exp_stats/exp_number/circRNA

genes='lncRNA protein_coding TUCP'
for gene in ${genes}
do
    vpy ../expression/gene_detect_number_by_mean.py \
       ${adipose_dir}/rmpad_exp_table/exp.${gene}.tpm.txt \
       ${adipose_dir}/../muscle/cfg_dir/${gene}.type.txt \
       ${adipose_dir}/cfg_dir/rmpad.reorder.adipose.sample.ini \
       0.1 \
       ${adipose_dir}/exp_stats/exp_number/${gene}

done

