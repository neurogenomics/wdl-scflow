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
     echo $OUTPUTPATH > outputpath.txt

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
    
    for files in qc_plot_data qc_plots qc_report qc_summary ${qc_key}_scflow_qc_report.html Rplots.pdf; do
         strato cp -r --backend ~{backend} -m $files $OUTPUTPATH/$files; done 
    
    >>>

     output {
     File count_depth = "qc_plot_data/~{qc_key}_count_depth_distribution.tsv"
     File qc_summary = "qc_summary/~{qc_key}_qc_summary.tsv"
     String outputpath = read_string("outputpath.txt")
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:0.1"
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
     docker: "eugeneduff/scflow-wdl:0.1"
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
          OUTPUTDIR=`dirname $INPUTDIR`
          TARGETDIR=`basename $INPUTDIR`
          strato cp -r --backend ~{backend} -m ${INPUTDIR} .
          FILEDIRS=${FILEDIRS}${TARGETDIR},
          ls ${TARGETDIR}
          #if [[ -z $FILEPATHS ]] ;
          #     then FILEPATHS=${OUTPUTPATH}
          #     else FILEPATHS=${FILEPATHS},$OUTPUTPATH
          #fi
     done
     
     # remove final character
     FILEDIRS=${FILEDIRS:0:-1}

     echo ./scflow_merge.r $options --sce_paths $FILEDIRS \
     --ensembl_mappings ~{ensembl_mappings} \
     --unique_id_var ~{qc_key_colname} \
     --plot_vars ~{sep="," plot_vars} \
     --facet_vars ~{sep="," facet_vars} \
     --outlier_vars ~{sep="," outlier_vars} \
     --species ~{species}
     ls $TARGETDIR
     ./scflow_merge.r $options --sce_paths $FILEDIRS --ensembl_mappings ~{ensembl_mappings} \
     --unique_id_var ~{qc_key_colname} \
     --plot_vars ~{sep="," plot_vars} \
     --facet_vars ~{sep="," facet_vars} \
     --outlier_vars ~{sep="," outlier_vars} \
     --species ~{species}
     

     for files in merged_sce merge_plots merge_summary_plots merged_report; do
          strato cp -r --backend ~{backend} -m $files $OUTPUTDIR/$files; done 
    
>>>

     output {
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
