version 1.0



workflow scflow_tasks {
}

task check_inputs {
	input {
	File manifest
	File input
	}

	command {
	Rscript -e check_inputs.r --input input --manifest manifest
	}

	output {
	File checked_manifest = checked_manifest.txt
	}

	runtime {
	container: "combiz/scflow-docker:0.6.1"
	}
}

task scflow_qc {
	input {
	File  
	File input
	
	}

	command {
	Rscript -e scflow_qc.r --input input --mat_path matpath --key key --ensemble_mappings ensemble_mappings --key_colname params.qc_key_colname} \
    --factor_vars params.qc_factor_vars} \
    --min_library_size params.qc_min_library_size} \
    --max_library_size params.qc_max_library_size} \
    --min_features params.qc_min_features} \
    --max_features params.qc_max_features} \
    --max_mito params.qc_max_mito} \
    --min_ribo params.qc_min_ribo} \
    --max_ribo params.qc_max_ribo} \
    --min_counts params.qc_min_counts} \
    --min_cells params.qc_min_cells} \
    --drop_unmapped params.qc_drop_unmapped} \
    --drop_mito params.qc_drop_mito} \
    --drop_ribo params.qc_drop_ribo} \
    --nmads params.qc_nmads} \
    --find_singlets params.mult_find_singlets} \
    --singlets_method params.mult_singlets_method} \
    --vars_to_regress_out params.mult_vars_to_regress_out} \
    --pca_dims params.mult_pca_dims} \
    --var_features params.mult_var_features} \
    --doublet_rate params.mult_doublet_rate} \
    --dpk params.mult_dpk} \
    --pK params.mult_pK} \
    --find_cells params.amb_find_cells} \
    --lower params.amb_lower} \
    --retain params.amb_retain} \
    --alpha_cutoff params.amb_alpha_cutoff} \
    --niters params.amb_niters \
    --expect_cells params.amb_expect_cells} \
    --species params.species}
	}

	output {
	File checked_manifest = checked_manifest.txt
	}

	runtime {
	container: "combiz/scflow-docker:0.6.1"
	}
}


