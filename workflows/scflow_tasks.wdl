version 1.0



workflow scflow_tasks {
}
# a struct is a user-defined type; it can contain members of any types

task check_inputs {
     input {
     File manifest_file
     File input_file
     }

     command {
     curl https://raw.githubusercontent.com/combiz/nf-core-scflow/0.7.0dev/bin/check_inputs.r > check_inputs.r;
     chmod +x *.r
     ./check_inputs.r --samplesheet ~{input_file}  --manifest ~{manifest_file}
}

     output {
     File checked_manifest = "checked_manifest.txt"
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:1.0"
     }
}

task scflow_qc {
     input {
     File input_file
     File ensembl_mappings
     Directory mat_path
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

     command {
     curl https://raw.githubusercontent.com/combiz/nf-core-scflow/dev/bin/scflow_qc.r  > scflow_qc.r;
     chmod +x *.r

     # extract data 
     # pull in data 

     

     ./scflow_qc.r --input ~{input_file} --mat_path ~{mat_path} --key ~{qc_key} --ensembl_mappings ~{ensembl_mappings} --key_colname ~{qc_key_colname} \
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
    --species ~{species} }

     output {
     File checked_manifest = "checked_manifest.txt"
     }

     runtime {
     docker: "eugeneduff/scflow-wdl:1.0"
     }
}


