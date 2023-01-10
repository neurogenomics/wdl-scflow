version 1.0

import "scflow_tasks.wdl" as tasks

struct BaseParams {
		File input_file
		File manifest
}

#struct PreProcParams {
#  String sample_name
#  File bam
#  File bam_index
#}

workflow scflow {
	input { 
        # Backend
        String backend
		File input_file 
		File manifest_file 
		File ensembl_mappings
		File ctd_path
		String     qc_key_colname
		String     qc_factor_vars
		Int qc_min_library_size
		String     qc_max_library_size
		Int     qc_min_features
		String     qc_max_features
		String     qc_max_mito
		Int     qc_min_ribo
		Int     qc_max_ribo
		Int     qc_min_counts
		Int     qc_min_cells
		Boolean     qc_drop_unmapped
		Boolean     qc_drop_mito
		Boolean     qc_drop_ribo
		Float     qc_nmads
		Boolean     amb_find_cells
		Int     amb_lower
		String     amb_retain
		Float     amb_alpha_cutoff
		Int     amb_niters
		Int     amb_expect_cells
		Boolean     mult_find_singlets
		String     mult_singlets_method
		String     mult_vars_to_regress_out
		Int     mult_pca_dims
		Int     mult_var_features
		Int     mult_doublet_rate
		Int     mult_dpk
		Float     mult_pK
		String     integ_method
		Int     integ_k
		String     integ_unique_id_var
		Boolean     integ_take_gene_union
		Boolean     integ_remove_missing
		Int     integ_num_genes
		String     integ_combine
		Boolean     integ_capitalize
		Boolean     integ_use_cols
		Float     integ_lambda
		Float     integ_thresh
		Int     integ_max_iters
		Int     integ_nrep
		Int     integ_rand_seed
		Int     integ_quantiles
		String     integ_ref_dataset
		Int     integ_min_cells
		Int     integ_knn_k
		Boolean     integ_center
		Int     integ_resolution
		Int     integ_max_cores
		Array[String]     integ_categorical_covariates
		String     integ_input_reduced_dim
		Array[String]     merge_plot_vars
		Array[String]     merge_facet_vars
		Array[String]     merge_outlier_vars
		Array[String]     reddim_input_reduced_dim
		Array[String]     reddim_reduction_methods
		Array[String]     reddim_vars_to_regress_out
		Int     reddim_umap_pca_dims
		Int     reddim_umap_n_neighbors
		Int     reddim_umap_n_components
		String     reddim_umap_init
		String     reddim_umap_metric
		Int     reddim_umap_n_epochs
		Int     reddim_umap_learning_rate
		Float     reddim_umap_min_dist
		Float     reddim_umap_spread
		Int     reddim_umap_set_op_mix_ratio
		Int     reddim_umap_local_connectivity
		Int     reddim_umap_repulsion_strength
		Int     reddim_umap_negative_sample_rate
		Boolean     reddim_umap_fast_sgd
		Int     reddim_tsne_dims
		Int     reddim_tsne_initial_dims
		Int     reddim_tsne_perplexity
		Float     reddim_tsne_theta
		Int     reddim_tsne_stop_lying_iter
		Int     reddim_tsne_mom_switch_iter
		Int     reddim_tsne_max_iter
		Boolean     reddim_tsne_pca_center
		Boolean     reddim_tsne_pca_scale
		Boolean     reddim_tsne_normalize
		Float     reddim_tsne_momentum
		Float     reddim_tsne_final_momentum
		Int     reddim_tsne_eta
		Int     reddim_tsne_exaggeration_factor
		String     clust_cluster_method
		String     clust_reduction_method
		Float     clust_res
		Int     clust_k
		Int     clust_louvain_iter
		String     cta_clusters_colname
		Int     cta_cells_to_sample
		String     cta_unique_id_var
		String     cta_clusters_colname
		String     cta_celltype_var
		Array[String]     cta_facet_vars
		Array[String]     cta_metric_vars
		Int     cta_top_n
		String     dge_de_method
		String     dge_mast_method
		Int     dge_min_counts
		Float     dge_min_cells_pc
		Boolean     dge_rescale_numerics
		Boolean     dge_pseudobulk
		String     dge_celltype_var
		String     dge_sample_var
		String     dge_dependent_var
		String     dge_ref_class
		String     dge_confounding_vars
		String     dge_random_effects_var
		Float     dge_fc_threshold
		Float     dge_pval_cutoff
		Int     dge_n_label
		Boolean     dge_force_run
		Int     dge_max_cores
		String     ipa_enrichment_tool
		String     ipa_enrichment_method
		String     ipa_enrichment_database
		String     dirich_unique_id_var
		String     dirich_celltype_var
		String     dirich_dependent_var
		String     dirich_ref_class
		String     dirich_var_order
		String     plotreddim_reduction_methods
		Float     reddimplot_pointsize
		Float     reddimplot_alpha
		String     species
		Int    max_cores
	}

	call tasks.check_inputs as check_inputs {
		input:
  			input_file = input_file,
  			manifest_file = manifest_file
			}


    scatter ( manifest_key in check_inputs.keys ) {

		call tasks.scflow_qc as scflow_qc {
			input:	
            backend = backend,
			input_file = input_file,
			manifest_file = manifest_file,
			ensembl_mappings = ensembl_mappings,
			qc_key_colname = qc_key_colname,
			qc_key = manifest_key,
			qc_factor_vars = qc_factor_vars,
			qc_min_library_size = qc_min_library_size,
			qc_max_library_size = qc_max_library_size,
			qc_min_features = qc_min_features,
			qc_max_features = qc_max_features,
			qc_max_mito = qc_max_mito,
			qc_min_ribo = qc_min_ribo,
			qc_max_ribo = qc_max_ribo,
			qc_min_counts = qc_min_counts,
			dge_min_counts = dge_min_counts,
			qc_min_cells = qc_min_cells,
			integ_min_cells = integ_min_cells,
			dge_min_cells_pc = dge_min_cells_pc,
			qc_drop_unmapped = qc_drop_unmapped,
			qc_drop_mito = qc_drop_mito,
			qc_drop_ribo = qc_drop_ribo,
			qc_nmads = qc_nmads,
			mult_find_singlets = mult_find_singlets,
			mult_singlets_method = mult_singlets_method,
			mult_vars_to_regress_out = mult_vars_to_regress_out,
			mult_pca_dims = mult_pca_dims,
			mult_var_features = mult_var_features,
			mult_doublet_rate = mult_doublet_rate,
			mult_pK = mult_pK,
			mult_dpk = mult_dpk,
			amb_find_cells = amb_find_cells,
			amb_lower = amb_lower,
			amb_retain = amb_retain,
			amb_alpha_cutoff = amb_alpha_cutoff,
			amb_niters = amb_niters,
			amb_expect_cells = amb_expect_cells,
			species = species
		}  

	}	

	call tasks.merge_qc as merge_qc {
		input: 
			qc_summaries = scflow_qc.qc_summary
			}
		
	call tasks.merge_sce as merge_sce {
		input:
		    merged_tsv = merge_qc.merged_tsv,
			sce_dirs = scflow_qc.outputpath,
			ensembl_mappings = ensembl_mappings,    
      		backend = backend,
     		qc_key_colname = qc_key_colname, 
     		plot_vars = merge_plot_vars ,
      		facet_vars = merge_facet_vars, 
     		outlier_vars = merge_outlier_vars, 
     		species =  species 
	}
	#call tasks.scflow_integrate {
    #    merged_sce.merged_sce
	#}

   
}

