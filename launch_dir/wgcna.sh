# 骨骼肌的共表达网络就再用没有心肌的做一次
# 1. extract exprssion table without heart

muscle_exp_dir="/project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/diff/muscle/"
sk_muscle_sample=${muscle_exp_dir}/cfg_dir/rename.muscle.sample.ids
sk_muscle_group=${muscle_exp_dir}/cfg_dir/rename.muscle.group
sk_muscle_exp_dir=${muscle_exp_dir}/exp_table/exp_genes/
venv_dir="/public/pipenv_dir/pig_atlas_project/"

cd /public/pipenv_dir/pig_atlas_project


genes="lncRNA protein_coding TUCP"
for gene in ${genes}
do
    pipenv run python ~/scripts/reorder_table_col.py  \
        ${sk_muscle_exp_dir}/${gene}.exp.tpm.txt \
	${sk_muscle_sample} \
	${sk_muscle_exp_dir}/without_heart/${gene}.tpm.sk.txt 

    pipenv run python ${venv_dir}/expression/base.py \
        expression_filter \
	--exp-table ${sk_muscle_exp_dir}/without_heart/${gene}.tpm.sk.txt \
	--cutoff 0.1 \
	--outfile ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.txt

    pipenv run python ${venv_dir}/expression/base.py \
        exp_by_group \
	--exp-table ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.txt \
	--group-inf ${muscle_exp_dir}/cfg_dir/rename.muscle.sample.ini\
	--outfile ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.grp.txt

done

genes="miRNA"
for gene in ${genes}
do
    pipenv run python ~/scripts/reorder_table_col.py  \
        ${sk_muscle_exp_dir}/${gene}.exp.tpm.txt \
	${sk_muscle_sample} \
	${sk_muscle_exp_dir}/without_heart/${gene}.tpm.sk.txt 

    pipenv run python ${venv_dir}/expression/base.py \
        expression_filter \
	--exp-table ${sk_muscle_exp_dir}/without_heart/${gene}.tpm.sk.txt \
	--cutoff 1 \
	--outfile ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.txt

    pipenv run python ${venv_dir}/expression/base.py \
        exp_by_group \
	--exp-table ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.txt \
	--group-inf ${muscle_exp_dir}/cfg_dir/rename.muscle.sample.ini\
	--outfile ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.grp.txt

done
genes="circRNA"
for gene in ${genes}
do
    pipenv run python ~/scripts/reorder_table_col.py  \
        ${sk_muscle_exp_dir}/${gene}.exp.tpm.txt \
	${sk_muscle_sample} \
	${sk_muscle_exp_dir}/without_heart/${gene}.tpm.sk.txt 

    pipenv run python ${venv_dir}/expression/base.py \
        expression_filter \
	--exp-table ${sk_muscle_exp_dir}/without_heart/${gene}.tpm.sk.txt \
	--cutoff 0.05 \
	--outfile ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.txt

    pipenv run python ${venv_dir}/expression/base.py \
        exp_by_group \
	--exp-table ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.txt \
	--group-inf ${muscle_exp_dir}/cfg_dir/rename.muscle.sample.ini\
	--outfile ${sk_muscle_exp_dir}/without_heart/${gene}.exp.tpm.sk.grp.txt

done


#2. run wgcna

/public/scripts/sys/nohuprun.sh /project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/wgcna/wgcna/muscle/without_heart/wgcna.sh


