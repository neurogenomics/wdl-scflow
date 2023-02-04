version 1.0



workflow scflow_tasks {
}
# sdf  a struct is a user-defined type; it can contain members of any types

task check_inputs {
     input {
     File manifest_file
     File input_file
     }

     command <<<
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/check_inputs.r > check_inputs.r;
     chmod +x check_inputs.r
     ls -ld /testing/*
     ./check_inputs.r --input ~{input_file}  --manifest ~{manifest_file}
     cat ~{manifest_file}  | awk '(NR>1)' | awk {' print $1 '}  > keys.txt
     
>>>

     output {
     File checked_manifest = "checked_manifest.txt"
     Array[String] keys = read_lines("keys.txt")

     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.1"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}


task scflow_qc {
     input {
     String backend
     File input_file
     File ensembl_mappings
     File manifest_file
     String     qc_key_colname
     String     qc_key
     String     qc_factor_vars
     Int     qc_min_library_size
     String     qc_max_library_size
     Int     qc_min_features
     String     qc_max_features
     String     qc_max_mito
     Int     qc_min_ribo
     Int     qc_max_ribo
     Int     qc_min_counts
     Int     dge_min_counts
     Int     qc_min_cells
     Int     integ_min_cells
     Float     dge_min_cells_pc
     Boolean     qc_drop_unmapped
     Boolean     qc_drop_mito
     Boolean     qc_drop_ribo
     Float     qc_nmads
     Boolean     mult_find_singlets
     String     mult_singlets_method
     String   mult_vars_to_regress_out
     Int     mult_pca_dims
     Int     mult_var_features
     Int     mult_doublet_rate
     Float     mult_pK
     Float     mult_dpk
     Boolean     amb_find_cells
     Int     amb_lower
     String     amb_retain
     Float     amb_alpha_cutoff
     Int     amb_niters
     Int     amb_expect_cells
     String     species
     }

     command <<<
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_qc.r  > scflow_qc.r;
     chmod +x scflow_qc.r
     
     mat_path=`cat ~{manifest_file} | grep ~{qc_key} | awk {' print $2 '}`

     if [[ "$mat_path" == *"zip" ]]; then
          MATPATH=mat_folder_"${qc_key}"
          strato cp --backend ~{backend} -m "$mat_path" $MATPATH
          unzip $MATPATH/*zip -d $MATPATH
          # removing the assumed .zip extension
          OUTPUTPATH=${mat_path%.zip} 
          
     else
          mkdir mat_path
          strato sync --backend ~{backend} -m "$mat_path" "mat_path"
          MATPATH=mat_path

          OUTPUTPATH=$mat_path
     fi
     echo $OUTPUTPATH/~{qc_key}_sce > outputpath.txt

     #if [[ -d mat_path ]]; then
     #   echo "${mat_path} is a directory"
     #   MATPATH=mat_path
     #else
     #   echo "${mat_path} is not valid"
     #   MATPATH="${mat_path}"
     #   exit 1

     #wget mat_path
     #unzip individual_1.zip -d ./mat_folder

     ./scflow_qc.r --input ~{input_file} --mat_path ${MATPATH} --key ~{qc_key} --ensembl_mappings ~{ensembl_mappings} --key_colname ~{qc_key_colname} \
    --factor_vars ~{qc_factor_vars} \
    --min_library_size ~{qc_min_library_size} \
    --max_library_size ~{qc_max_library_size} \
    --min_features ~{qc_min_features} \
    --max_features ~{qc_max_features} \
    --max_mito ~{qc_max_mito} \
    --min_ribo ~{qc_min_ribo} \
    --max_ribo ~{qc_max_ribo} \
    --min_counts ~{qc_min_counts} \
    --min_cells ~{qc_min_cells} \
    --drop_unmapped ~{qc_drop_unmapped} \
    --drop_mito ~{qc_drop_mito} \
    --drop_ribo ~{qc_drop_ribo} \
    --nmads ~{qc_nmads} \
    --find_singlets ~{mult_find_singlets} \
    --singlets_method ~{mult_singlets_method} \
    --vars_to_regress_out ~{mult_vars_to_regress_out} \
    --pca_dims ~{mult_pca_dims} \
    --var_features ~{mult_var_features} \
    --doublet_rate ~{mult_doublet_rate} \
    --dpk ~{mult_dpk} \
    --pK ~{mult_pK} \
    --find_cells ~{amb_find_cells} \
    --lower ~{amb_lower} \
    --retain ~{amb_retain} \
    --alpha_cutoff ~{amb_alpha_cutoff} \
    --niters ~{amb_niters} \
    --expect_cells ~{amb_expect_cells} \
    --species ~{species} 
    
    for files in qc_plot_data qc_plots qc_report qc_summary ~{qc_key}_scflow_qc_report.html Rplots.pdf ~{qc_key}_sce; do
         strato sync --backend ~{backend} -m $files $OUTPUTPATH/$files; done 
    
    >>>

     output {
     File count_depth = "qc_plot_data/~{qc_key}_count_depth_distribution.tsv"
     File qc_summary = "qc_summary/~{qc_key}_qc_summary.tsv"
     String outputpath = read_string("outputpath.txt")
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "120G"
     bootDiskSizeGb: "16"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}


task merge_qc {
     input {
     Array[File] qc_summaries
     }

     command <<<
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/merge_tables.r > merge_tables.r;
     chmod +x merge_tables.r
     ./merge_tables.r --filepaths ~{sep="," qc_summaries}
>>>

     output {
     File merged_tsv = "merged.tsv"
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}

task merge_sce {
     input {
     File merged_tsv
     Array[String] sce_dirs 
     File ensembl_mappings 
     String backend
     String qc_key_colname 
     Array[String] plot_vars 
     Array[String] facet_vars 
     Array[String] outlier_vars 
     String species 
     }

     command <<<
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_merge.r > scflow_merge.r;
     chmod +x scflow_merge.r

     FILEPATHS=''
     for INPUTDIR in ~{sep=" " sce_dirs}; do 
          TARGETDIR=`basename $INPUTDIR`
          strato cp -r --backend ~{backend} -m ${INPUTDIR} .
          FILEDIRS=${FILEDIRS}${TARGETDIR},
          #if [[ -z $FILEPATHS ]] ;
          #     then FILEPATHS=${OUTPUTPATH}
          #     else FILEPATHS=${FILEPATHS},$OUTPUTPATH
          #fi
     done
     
     # remove final character
     FILEDIRS=${FILEDIRS:0:-1}
     # generate group outputdir (two dirs back should be variable)
     OUTPUTDIR=`dirname $INPUTDIR`
     OUTPUTDIR=`dirname $OUTPUTDIR`

     ./scflow_merge.r --sce_paths $FILEDIRS --ensembl_mappings ~{ensembl_mappings} \
     --unique_id_var ~{qc_key_colname} \
     --plot_vars ~{sep="," plot_vars} \
     --facet_vars ~{sep="," facet_vars} \
     --outlier_vars ~{sep="," outlier_vars} \
     --species ~{species}

     echo $OUTPUTDIR/merged_sce > outputpath.txt

     for files in merged_sce merge_plots merge_summary_plots merged_report; do
          strato sync --backend ~{backend} -m $files $OUTPUTDIR/$files; done 
    
>>>

     output {
     String merged_sce = read_string("outputpath.txt")    }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}

task scflow_integrate {
     input {
          String merged_sce
          String method
          String backend
		Int k
		String unique_id_var
		Boolean take_gene_union
		Boolean remove_missing
		Int num_genes
		String combine
		Boolean capitalize
		Boolean use_cols
		Float lambda
		Float thresh
		Int max_iters
		Int nrep
		Int rand_seed
		Int quantiles
		String ref_dataset
		Int min_cells
		Int knn_k
		Boolean center
		Float resolution 
     }

     command <<<
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_integrate.r > scflow_integrate.r;
     chmod +x scflow_integrate.r

     
     strato cp -r --backend ~{backend} -m ~{merged_sce} .
     INPUTDIR=`basename ~{merged_sce} `

     ./scflow_integrate.r --sce_path $INPUTDIR \
     --method ~{method} \
    --k ~{k} \
    --unique_id_var ~{unique_id_var} \
    --take_gene_union ~{take_gene_union} \
    --remove_missing ~{remove_missing} \
    --num_genes ~{num_genes} \
    --combine ~{combine} \
    --capitalize ~{capitalize} \
    --use_cols ~{use_cols} \
    --lambda ~{lambda} \
    --thresh ~{thresh} \
    --max_iters ~{max_iters} \
    --nrep ~{nrep} \
    --rand_seed ~{rand_seed} \
    --quantiles ~{quantiles} \
    --ref_dataset ~{ref_dataset} \
    --min_cells ~{min_cells} \
    --knn_k ~{knn_k} \
    --center ~{center} \
    --resolution ~{resolution}
     
     OUTPUTDIR=`dirname ~{merged_sce}`

     strato sync --backend ~{backend} -m integrated_sce $OUTPUTDIR/integrated_sce
     echo $OUTPUTDIR/integrated_sce > outputpath.txt
>>>

     output {
     String integrated_sce = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }

}



task scflow_reduce_dims {
     input {
    String integrated_sce
    Array[String] reddim_input_reduced_dim
    Array[String] reddim_reduction_methods
    Array[String] reddim_vars_to_regress_out
    Int reddim_umap_pca_dims
    Int reddim_umap_n_neighbors
    Int reddim_umap_n_components
    String reddim_umap_init
    String reddim_umap_metric
    Int reddim_umap_n_epochs
    Float reddim_umap_learning_rate
    Float reddim_umap_min_dist
    Float reddim_umap_spread
    Float reddim_umap_set_op_mix_ratio
    Int reddim_umap_local_connectivity
    Int reddim_umap_repulsion_strength
    Int reddim_umap_negative_sample_rate
    String reddim_umap_fast_sgd
    Int reddim_tsne_dims
    Int reddim_tsne_initial_dims
    Int reddim_tsne_perplexity
    Float reddim_tsne_theta
    Int reddim_tsne_stop_lying_iter
    Int reddim_tsne_mom_switch_iter
    Int reddim_tsne_max_iter
    Boolean reddim_tsne_pca_center
    Boolean reddim_tsne_pca_scale
    Boolean reddim_tsne_normalize
    Float reddim_tsne_momentum
    Float reddim_tsne_final_momentum
    Int reddim_tsne_eta
    Int reddim_tsne_exaggeration_factor
    String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_reduce_dims.r > scflow_reduce_dims.r;
     chmod +x scflow_reduce_dims.r  

     export MC_CORES=4
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{integrated_sce} .
     INPUTDIR=`basename ~{integrated_sce} `

     ./scflow_reduce_dims.r --sce_path $INPUTDIR \
    --input_reduced_dim ~{sep="," reddim_input_reduced_dim} \
    --reduction_methods ~{sep="," reddim_reduction_methods} \
    --vars_to_regress_out ~{sep="," reddim_vars_to_regress_out} \
    --pca_dims ~{reddim_umap_pca_dims} \
    --n_neighbors ~{reddim_umap_n_neighbors} \
    --n_components ~{reddim_umap_n_components} \
    --init ~{reddim_umap_init} \
    --metric ~{reddim_umap_metric} \
    --n_epochs ~{reddim_umap_n_epochs} \
    --learning_rate ~{reddim_umap_learning_rate} \
    --min_dist ~{reddim_umap_min_dist} \
    --spread ~{reddim_umap_spread} \
    --set_op_mix_ratio ~{reddim_umap_set_op_mix_ratio} \
    --local_connectivity ~{reddim_umap_local_connectivity} \
    --repulsion_strength ~{reddim_umap_repulsion_strength} \
    --negative_sample_rate ~{reddim_umap_negative_sample_rate} \
    --fast_sgd ~{reddim_umap_fast_sgd} \
    --dims ~{reddim_tsne_dims} \
    --initial_dims ~{reddim_tsne_initial_dims} \
    --perplexity ~{reddim_tsne_perplexity} \
    --theta ~{reddim_tsne_theta} \
    --stop_lying_iter ~{reddim_tsne_stop_lying_iter} \
    --mom_switch_iter ~{reddim_tsne_mom_switch_iter} \
    --max_iter ~{reddim_tsne_max_iter} \
    --pca_center ~{reddim_tsne_pca_center} \
    --pca_scale ~{reddim_tsne_pca_scale} \
    --normalize ~{reddim_tsne_normalize} \
    --momentum ~{reddim_tsne_momentum} \
    --final_momentum ~{reddim_tsne_final_momentum} \
    --eta ~{reddim_tsne_eta} \
    --exaggeration_factor ~{reddim_tsne_exaggeration_factor}
     
     OUTPUTDIR=`dirname ~{integrated_sce}`

     strato sync --backend ~{backend} -m reddim_sce $OUTPUTDIR/reddim_sce
     echo $OUTPUTDIR/reddim_sce > outputpath.txt
>>>

     output {
     String reddim_sce = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }

}

task scflow_cluster {
     input {
    String sce_path
    String clust_cluster_method
    String clust_reduction_method
    Float clust_res
    Int clust_k
    Int clust_louvain_iter
    String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_cluster.r > scflow_cluster.r;
     chmod +x scflow_cluster.r  

     export MC_CORES=4
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     ./scflow_cluster.r --sce_path $INPUTDIR \
     --cluster_method ~{clust_cluster_method} \
    --reduction_method ~{clust_reduction_method} \
    --res ~{clust_res} \
    --k ~{clust_k} \
    --louvain_iter ~{clust_louvain_iter}
     
     OUTPUTDIR=`dirname ~{sce_path}`

     strato sync --backend ~{backend} -m clustered_sce $OUTPUTDIR/clustered_sce
     echo $OUTPUTDIR/clustered_sce > outputpath.txt
>>>

     output {
     String clustered_sce = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }

}


task scflow_report_integrated {
     input {
    String sce_path
    Array[String] integ_categorical_covariates
    String integ_input_reduced_dim
    Float reddimplot_pointsize
    Float reddimplot_alpha
    String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_report_integrated.r > scflow_report_integrated.r ;
     chmod +x scflow_report_integrated.r   

     export MC_CORES=4
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     ./scflow_report_integrated.r --sce_path $INPUTDIR \
     --categorical_covariates ~{sep="," integ_categorical_covariates} \
    --input_reduced_dim ~{integ_input_reduced_dim} \
    --reddimplot_pointsize ~{reddimplot_pointsize} \
    --reddimplot_alpha ~{reddimplot_alpha}

     OUTPUTDIR=`dirname ~{sce_path}`

     strato sync --backend ~{backend} -m integration_report $OUTPUTDIR/integration_report
     echo $OUTPUTDIR/integration_report > outputpath.txt
>>>

     output {
     String integration_report = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}


task scflow_map_celltypes {
     input {
    String sce_path
    File ctd_path
    String cta_clusters_colname
    Int cta_cells_to_sample
    String species
    Float reddimplot_pointsize
    Float reddimplot_alpha
    String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_map_celltypes.r > scflow_map_celltypes.r ;
     chmod +x scflow_map_celltypes.r   

     export MC_CORES=4
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     mkdir ctd_folder && unzip ~{ctd_path} -d ./ctd_folder

     ./scflow_map_celltypes.r  --sce_path $INPUTDIR --ctd_folder ./ctd_folder \
     --clusters_colname ~{cta_clusters_colname} \
    --cells_to_sample ~{cta_cells_to_sample} \
    --species ~{species} \
    --reddimplot_pointsize ~{reddimplot_pointsize} \
    --reddimplot_alpha ~{reddimplot_alpha}

     OUTPUTDIR=`dirname ~{sce_path}`

     strato sync --backend ~{backend} -m celltype_mapped_sce $OUTPUTDIR/celltype_mapped_sce
     echo $OUTPUTDIR/celltype_mapped_sce > outputpath.txt
>>>

     output {
     String celltype_mapped_sce = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}


task scflow_finalize_sce {
     input {
    String sce_path
    String celltype_mappings
    String cta_clusters_colname
    String cta_celltype_var
    String cta_unique_id_var
    Array[String] cta_facet_vars
    String clust_reduction_method
    Array[String] cta_metric_vars
    Int cta_top_n
    Float reddimplot_pointsize
    Float reddimplot_alpha
    Int max_cores
    String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_finalize_sce.r > scflow_finalize_sce.r ;
     chmod +x scflow_finalize_sce.r   

     export MC_CORES=~{max_cores}
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     strato exists --backend ~{backend} ~{celltype_mappings} && strato cp --backend ~{backend} ~{celltype_mappings} celltype_mappings.tsv 

     ./scflow_finalize_sce.r  --sce_path $INPUTDIR --celltype_mappings celltype_mappings.tsv  \
     --clusters_colname ~{cta_clusters_colname} \
     --celltype_var ~{cta_celltype_var} \
     --unique_id_var ~{cta_unique_id_var} \
     --facet_vars ~{sep="," cta_facet_vars} \
     --input_reduced_dim ~{clust_reduction_method} \
     --metric_vars ~{sep="," cta_metric_vars} \
     --top_n ~{cta_top_n} \
     --reddimplot_pointsize ~{reddimplot_pointsize} \
     --reddimplot_alpha ~{reddimplot_alpha} \
     --max_cores ~{max_cores}

     OUTPUTDIR=`dirname ~{sce_path}`

     for DIR in celltype_metrics_report celltype_mapped_sce celltype_marker_plots  celltype_marker_tables final_sce; do
          strato sync --backend ~{backend} -m $DIR $OUTPUTDIR/$DIR
     done

     echo $OUTPUTDIR/final_sce > outputpath.txt

     cat celltypes.tsv | awk '(NR>1)' | awk {' print $1 '} > celltype.txt
  
>>>

     output {
     String final_sce = read_string("outputpath.txt") 
     Array[String] celltypes = read_lines("celltype.txt")
     File celltypes_n_cells = "celltypes.tsv"
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "24"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }

}


task scflow_dge {
     input {
    String sce_path
    String celltype
    File celltypes_n_cells
    String dge_de_method
    String dge_mast_method
    File ensembl_mappings
    Int dge_min_counts
    Float dge_min_cells_pc
    Boolean dge_rescale_numerics
    Boolean dge_force_run
    Boolean dge_pseudobulk
    String dge_celltype_var
    String dge_sample_var
    String dge_dependent_var
    String dge_ref_class
    String dge_confounding_vars
    String dge_random_effects_var
    Float dge_pval_cutoff
    Float dge_fc_threshold
    String species
    Int max_cores
    String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_dge.r > scflow_dge.r ;
     chmod +x scflow_dge.r   

     export MC_CORES=~{max_cores}
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     cat ~{celltypes_n_cells} | grep Micro | awk '{ print "celltype: " $1 " n_cells:" $2 }' 

     ./scflow_dge.r  --sce $INPUTDIR --celltype ~{celltype}  \
     --de_method ~{dge_de_method} \
     --ensembl_mappings ~{ensembl_mappings} \
     --mast_method ~{dge_mast_method} \
     --min_counts ~{dge_min_counts} \
     --min_cells_pc ~{dge_min_cells_pc} \
     --rescale_numerics ~{dge_rescale_numerics} \
     --force_run ~{dge_force_run} \
     --pseudobulk ~{dge_pseudobulk} \
     --celltype_var ~{dge_celltype_var} \
     --sample_var ~{dge_sample_var} \
     --dependent_var ~{dge_dependent_var} \
     --ref_class ~{dge_ref_class} \
     --confounding_vars ~{dge_confounding_vars} \
     --random_effects_var ~{dge_random_effects_var} \
     --pval_cutoff ~{dge_pval_cutoff} \
     --fc_threshold ~{dge_fc_threshold} \
     --species ~{species} \
     --max_cores ~{max_cores}

     OUTPUTDIR=`dirname ~{sce_path}`/DGE

     for DIR in *_DE.tsv *html *_volcano_plot.png ; do
          strato cp -r --backend ~{backend} -m $DIR $OUTPUTDIR/
     done

     echo $OUTPUTDIR/integration_report > outputpath.txt
     echo $OUTPUTDIR/`ls *_DE.tsv` > output_DE.txt
     echo $OUTPUTDIR/`ls *_de_report.html` > output_report.txt
     echo $OUTPUTDIR/`ls *_volcano_plot.png` > output_gg.txt
>>>

     output {
     String de_table = read_string("output_DE.txt") 
     String de_report = read_string("output_report.txt") 
     String de_plot = read_string("output_gg.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "36G"
     bootDiskSizeGb: "36"
     disks: "local-disk 100 HDD"
     cpu: 20
     preemptible: 1
     maxRetries: 5
     }

}


task scflow_ipa {
     input {
    String de_table
    String ipa_enrichment_tool
    String ipa_enrichment_method
    String ipa_enrichment_database
    Float dge_pval_cutoff
    Float dge_fc_threshold
    String species
    String backend
     }

     command <<<
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_ipa.r > scflow_ipa.r ;
     chmod +x scflow_ipa.r   

     strato cp -r --backend ~{backend} -m ~{de_table} .
     INPUTDIR=`basename ~{de_table} `

     export MC_CORES=4
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     ./scflow_ipa.r  --gene_file ${INPUTDIR} \
     --enrichment_tool ~{ipa_enrichment_tool} \
    --enrichment_method ~{ipa_enrichment_method} \
    --enrichment_database ~{ipa_enrichment_database} \
    --pval_cutoff ~{dge_pval_cutoff} \
    --fc_threshold ~{dge_fc_threshold} \
    --species ~{species}

     
     de_table_filename=`basename ${INPUTDIR}` 
     ipa_name="${de_table_filename%.*}"_ipa
     mv ipa ${ipa_name}
     OUTPUTDIR=`dirname ~{de_table}`

     strato cp -r --backend ~{backend} -m ${ipa_name} $OUTPUTDIR/${ipa_name}
     
     echo $OUTPUTDIR/${ipa_name} > outputpath.txt
>>>

     output {
     String celltype_mapped_sce = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }
}



task scflow_dirichlet {
     input {
     String sce_path
     String dirich_unique_id_var
     String dirich_celltype_var
     String dirich_dependent_var
     String dirich_ref_class
     String dirich_var_order
     String backend
     }

     command <<<
     
     
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_dirichlet.r > scflow_dirichlet.r ;
     chmod +x scflow_dirichlet.r   

     export MC_CORES=1
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     ./scflow_dirichlet.r  --sce_path $INPUTDIR  \
     --unique_id_var ~{dirich_unique_id_var} \
    --celltype_var ~{dirich_celltype_var} \
    --dependent_var ~{dirich_dependent_var} \
    --ref_class ~{dirich_ref_class} \
    --var_order ~{dirich_var_order}

     OUTPUTDIR=`dirname ~{sce_path}`

     strato sync --backend ~{backend} -m dirichlet_report $OUTPUTDIR/dirichlet_report

     echo $OUTPUTDIR/dirichlet_report > outputpath.txt
>>>

     output {
     String dirichlet_report = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }

}

task scflow_plot_reddim_genes {
     input {
     String sce_path
     String plotreddim_reduction_methods     
     Float reddimplot_pointsize
     Float reddimplot_alpha
     File reddim_genes_yml
     String backend
     }

     command <<<
      
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/scflow_plot_reddim_genes.r > scflow_plot_reddim_genes.r ;
     chmod +x scflow_plot_reddim_genes.r   

     export MC_CORES=1
     export MKL_NUM_THREADS=1
     export NUMEXPR_NUM_THREADS=1
     export OMP_NUM_THREADS=1
     export OPENBLAS_NUM_THREADS=1
     export VECLIB_MAXIMUM_THREADS=1

     strato cp -r --backend ~{backend} -m ~{sce_path} .
     INPUTDIR=`basename ~{sce_path} `

     ./scflow_plot_reddim_genes.r  --sce_path $INPUTDIR \
          --reddim_genes_yml ~{reddim_genes_yml} \
          --reduction_methods ~{plotreddim_reduction_methods} \
          --reddimplot_pointsize ~{reddimplot_pointsize} \
          --reddimplot_alpha ~{reddimplot_alpha}

     OUTPUTDIR=`dirname ~{sce_path}`/

     strato sync --backend ~{backend} -m reddim_gene_plots $OUTPUTDIR/reddim_gene_plots

     echo $OUTPUTDIR/reddim_gene_plots > outputpath.txt
>>>

     output {
     String reddim_gene_plots = read_string("outputpath.txt") 
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.11"
     memory: "12G"
     bootDiskSizeGb: "12"
     disks: "local-disk 100 HDD"
     cpu: 1
     preemptible: 1
     maxRetries: 0
     }

}