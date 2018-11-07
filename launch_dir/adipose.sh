alias vpy='pipenv run python'
echo start
date

proj_dir='/project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/'
script_dir='/public/pipenv_dir/pig_atlas_project/'

##cfg
adipose_dir=/project0/OM-mRNA-pig-limingzhou-P160901/sus_scrofa/diff/adipose/
adipose_estats_dir=${adipose_dir}/exp_stats/
adipose_exp_dir=${adipose_dir}/exp_table
adipose_r_exp_dir=${adipose_dir}/exp_table_all
adipose_cfg_dir=${adipose_dir}/cfg_dir/
rmpad_adipose_group=${adipose_cfg_dir}/rmpad.reorder.adipose.sample.ini
adipose_group=${adipose_cfg_dir}/reorder.adipose.sample.ini
l_adipose_group=${adipose_cfg_dir}/long.reorder.adipose.sample.ini
s_adipose_group=${adipose_cfg_dir}/small.reorder.adipose.sample.ini
adipose_mi_group=${adipose_cfg_dir}/rmpad.mi.reorder.adipose.sample.ini
adipose_sample=${adipose_cfg_dir}/sample.order
l_adipose_sample=${adipose_cfg_dir}/long.sample.order
rmpad_adipose_group_group=${adipose_cfg_dir}/rmpad.adipose.group.ini
adipose_group_group=${adipose_cfg_dir}/adipose.groups.ini



# extract sample exp table
gene_types="lncRNA TUCP protein_coding circRNA"
for gene_type in ${gene_types}
do
    outfile=${adipose_exp_dir}/${gene_type}.tpm.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/reorder_table_col.py \
	    ${adipose_r_exp_dir}/${gene_type}.tpm.txt \
	    ${l_adipose_sample} \
	    ${adipose_exp_dir}/${gene_type}.tpm.txt
    fi
    outfile=${adipose_exp_dir}/${gene_type}.count.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/reorder_table_col.py \
	    ${adipose_r_exp_dir}/${gene_type}.count.txt \
	    ${l_adipose_sample} \
	    ${adipose_exp_dir}/${gene_type}.count.txt
    fi
done

gene_types="miRNA"
for gene_type in ${gene_types}
do
    outfile=${adipose_exp_dir}/${gene_type}.tpm.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/reorder_table_col.py \
	    ${adipose_r_exp_dir}/${gene_type}.tpm.txt \
	    ${adipose_sample} \
	    ${adipose_exp_dir}/${gene_type}.tpm.txt
    fi
    outfile=${adipose_exp_dir}/${gene_type}.count.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/reorder_table_col.py \
	    ${adipose_r_exp_dir}/${gene_type}.count.txt \
	    ${adipose_sample} \
	    ${adipose_exp_dir}/${gene_type}.count.txt
    fi
done


# Sun Sep 30 05:38:05 EDT 2018

gene_types='miRNA'
cutoff=1

for gene_type in ${gene_types}
do
    outfile=${adipose_exp_dir}/exp.${gene_type}.tpm.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/base.py expression-filter \
	    --exp-table ${adipose_exp_dir}/${gene_type}.tpm.txt \
	    --cutoff ${cutoff} \
	    --outfile ${outfile}
    fi
    outfile=${adipose_exp_dir}/exp.${gene_type}.count.txt
    if [ ! -s ${outfile} ]
    then
	cut -f1 ${adipose_exp_dir}/exp.${gene_type}.tpm.txt > ${adipose_exp_dir}/exp.${gene_type}.gene.list
	vpy /public/scripts/omsTools/omstools/general/extract_info_by_id.py \
	    --id ${adipose_exp_dir}/exp.${gene_type}.gene.list \
	    --table ${adipose_exp_dir}/${gene_type}.count.txt \
	    --output ${outfile} \
	    --header yes
    fi
done

gene_types='circRNA'
cutoff=0.05

for gene_type in ${gene_types}
do
    outfile=${adipose_exp_dir}/exp.${gene_type}.tpm.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/base.py expression-filter \
	    --exp-table ${adipose_exp_dir}/${gene_type}.tpm.txt \
	    --cutoff ${cutoff} \
	    --outfile ${outfile}
    fi
    outfile=${adipose_exp_dir}/exp.${gene_type}.count.txt
    if [ ! -s ${outfile} ]
    then
	cut -f1 ${adipose_exp_dir}/exp.${gene_type}.tpm.txt > ${adipose_exp_dir}/exp.${gene_type}.gene.list
	vpy /public/scripts/omsTools/omstools/general/extract_info_by_id.py \
	    --id ${adipose_exp_dir}/exp.${gene_type}.gene.list \
	    --table ${adipose_exp_dir}/${gene_type}.count.txt \
	    --output ${outfile} \
	    --header yes
    fi
done

gene_types='lncRNA TUCP protein_coding'
cutoff=0.1

for gene_type in ${gene_types}
do
    outfile=${adipose_exp_dir}/exp.${gene_type}.tpm.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/base.py expression-filter \
	    --exp-table ${adipose_exp_dir}/${gene_type}.tpm.txt \
	    --cutoff ${cutoff} \
	    --outfile ${outfile}
    fi
    outfile=${adipose_exp_dir}/exp.${gene_type}.count.txt
    if [ ! -s ${outfile} ]
    then
	cut -f1 ${adipose_exp_dir}/exp.${gene_type}.tpm.txt > ${adipose_exp_dir}/exp.${gene_type}.gene.list
	vpy /public/scripts/omsTools/omstools/general/extract_info_by_id.py \
	    --id ${adipose_exp_dir}/exp.${gene_type}.gene.list \
	    --table ${adipose_exp_dir}/${gene_type}.count.txt \
	    --output ${outfile} \
	    --header yes
    fi
done


gene_types="lncRNA TUCP protein_coding miRNA circRNA"
cluster_method='pearson spearman euclidean'
cluster_dir=${proj_dir}/cluster/adipose

for gene_type in ${gene_types}
do
    #### cluster analysis
    for method in ${cluster_method}
    do
        outfile=${cluster_dir}/adipose.${gene_type}.${method}.pcluster.png
        if [ ! -s ${outfile} ]
        then
            Rscript ${script_dir}/expression/pcluster_cor.R \
        	    ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
        	    adipose.${gene_type} \
        	    ${method} \
        	    ${adipose_group} \
        	    ${cluster_dir}
        fi
    done

    outfile=${adipose_exp_dir}/exp.${gene_type}.tpm.grp.txt
    if [ ! -s ${outfile} ]
    then
        
        vpy ${script_dir}/expression/base.py exp_by_group \
            --exp-table ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
            --group-inf ${adipose_group} \
            --outfile ${outfile}
    fi
    
    #### cluster analysis for group mean expression
    
    for method in ${cluster_method}
    do
        outfile=${cluster_dir}/adipose.grp.${gene_type}.${method}.pcluster.png
        if [ ! -s ${outfile} ]
        then
            Rscript ${script_dir}/expression/pcluster_cor.R \
        	    ${adipose_exp_dir}/exp.${gene_type}.tpm.grp.txt \
        	    adipose.grp.${gene_type} \
        	    ${method} \
        	    ${adipose_group_group} \
        	    ${cluster_dir}
        fi
    done
done


# diff analysis
diff_dir=${adipose_dir}/rmpad_diff_dir/
muscle_data_dir=${proj_dir}/diff/muscle/
genes='protein_coding lncRNA TUCP circRNA'
for gene_type in ${genes}
do  
    outfile=${diff_dir}/${gene_type}/${gene_type}.de_number.matrix.txt
    if [ ! -s ${outfile} ]
    then
    vpy ${script_dir}/expression/launch_diff_coroutine.py \
	--sample-inf ${rmpad_adipose_group} \
	--counts ${adipose_exp_dir}/exp.${gene_type}.count.txt \
	--tpm-table ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
	--out-dir ${diff_dir}/${gene_type} \
	--gene_class ${muscle_data_dir}/cfg_dir/${gene_type}.type.txt \
	--lib-size ${muscle_data_dir}/cfg_dir/lib_size.all.txt \
	- launch-diff

    vpy ${script_dir}/expression/launch_diff.py \
        --sample_inf ${rmpad_adipose_group} \
        --counts ${adipose_exp_dir}/exp.${gene_type}.count.txt \
        --tpm_table ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
        --out_dir ${diff_dir}/${gene_type} \
        --gene_class ${muscle_data_dir}/cfg_dir/${gene_type}.type.txt
    fi
done

genes='miRNA'
for gene_type in ${genes}
do
    outfile=${diff_dir}/${gene_type}/${gene_type}.de_number.matrix.txt
    if [ ! -s ${outfile} ]
    then
    vpy ${script_dir}/expression/launch_diff_coroutine.py \
	--sample-inf ${rmpad_adipose_group} \
	--counts ${adipose_exp_dir}/exp.${gene_type}.count.txt \
	--tpm-table ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
	--out-dir ${diff_dir}/${gene_type} \
	--gene_class ${muscle_data_dir}/cfg_dir/${gene_type}.type.txt \
	- launch-diff

    vpy ${script_dir}/expression/launch_diff.py \
        --sample_inf ${rmpad_adipose_group} \
        --counts ${adipose_exp_dir}/exp.${gene_type}.count.txt \
        --tpm_table ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
        --out_dir ${diff_dir}/${gene_type} \
        --gene_class ${muscle_data_dir}/cfg_dir/${gene_type}.type.txt
    fi
done

# Mon Oct  1 21:34:59 EDT 2018

gene_types="lncRNA TUCP protein_coding circRNA miRNA"
cutoff=0

ts_dir=${proj_dir}/tissue_specific/

for gene_type in ${gene_types}
do
    outfile=${ts_dir}/tsi/adipose/${gene_type}/tsi.score.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/ts.py \
	    --matrix ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
	    --group ${adipose_group} \
	    --gene_classify ${muscle_data_dir}/cfg_dir/${gene_type}.type.txt \
	    --out_dir ${ts_dir}/tsi/adipose/${gene_type} \
	    --exp_cut ${cutoff} 
    fi
    outfile=${ts_dir}/shannon_entropy/adipose/${gene_type}/shannon_entropy.score.txt
    if [ ! -s ${outfile} ]
    then
	vpy ${script_dir}/expression/ts_shannon_entropy.py \
	    --matrix ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
	    --group ${adipose_group} \
	    --gene_classify ${muscle_data_dir}/cfg_dir/${gene_type}.type.txt \
	    --out_dir ${ts_dir}/shannon_entropy/adipose/${gene_type} \
	    --exp_cut ${cutoff} 
    fi
    
done

##adipose ts genes cluster
gene_types="lncRNA TUCP protein_coding circRNA miRNA"
for gene_type in ${gene_types}
do
    ts_genes=${ts_dir}/tsi/adipose/${gene_type}/tissue_specific.genes.list
    if [ ! -s ${ts_genes} ]
    then
	awk 'NR>1{print $1}' ${ts_dir}/tsi/adipose/${gene_type}/tissue_specific.genes.txt > ${ts_genes}
    fi

    for method in ${cluster_method}
    do
        outfile=${cluster_dir}/ts.adipose.${gene_type}.${method}.pcluster.png
	if [ ! -s ${outfile} ]
	then
	    Rscript ${script_dir}/expression/pcluster_cor_param.R \
		    --file_path ${adipose_exp_dir}/exp.${gene_type}.tpm.txt \
		    --name ts.adipose.${gene_type} \
		    --method ${method} \
		    --sample_inf ${adipose_group} \
		    --out_dir ${cluster_dir} \
		    --genes ${ts_genes}
	fi

        outfile=${cluster_dir}/ts.grp.adipose.${gene_type}.${method}.pcluster.png
	if [ ! -s ${outfile} ]
	then
	    Rscript ${script_dir}/expression/pcluster_cor_param.R \
		    --file_path ${adipose_exp_dir}/exp.${gene_type}.tpm.grp.txt \
		    --name ts.grp.adipose.${gene_type} \
		    --method ${method} \
        	    --sample_inf ${adipose_group_group} \
		    --out_dir ${cluster_dir} \
		    --genes ${ts_genes}
	fi
    done

done

# Fri Oct 26 07:25:49 EDT 2018


date
echo finished
