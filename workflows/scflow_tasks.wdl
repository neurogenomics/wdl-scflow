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
     echo hihi
     curl https://raw.githubusercontent.com/neurogenomics/wdl-scflow/master/workflows/r/check_inputs.r > check_inputs.r;
     chmod +x *.r
     echo ./check_inputs.r --input ~{input_file}  --manifest ~{manifest_file}
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
     chmod +x *.r

     
     mat_path=`cat ~{manifest_file} | grep ~{qc_key} | awk {' print $2 '}`
     echo mat_path $mat_path
     
     mkdir mat_path
     strato sync --backend gcp -m "$mat_path" "mat_path"

     if [[ -d mat_path ]]; then
        echo "${mat_path} is a directory"
        MATPATH=mat_path
    elif [[ -f mat_path ]]; then
        echo "${mat_path} is a file, unzipping"
        mkdir mat_folder && unzip mat_path -d ./mat_folder
        MATPATH=mat_folder
    else
        echo "${mat_path} is not valid"
        MATPATH="${mat_path}"
        exit 1
    fi


     echo MATPATH $MATPATH
     ls $MATPATH


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
    
    echo ls 
    ls

    strato cp -r --backend gcp -m qc_plot_data $mat_path/qc_plot_data

    >>>

     output {
     File count_depth = "qc_plot_data/"+ qc_key +"_count_depth_distribution.tsv"
     File tmp = "tmp.txt"
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


