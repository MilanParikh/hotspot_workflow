version 1.0
workflow hotspot {
    input {
    	String output_directory
        File anndata_file
        String hotspot_model = 'danb'
        Int n_neighbors = 30
        Int min_gene_threshold = 30
        Float hs_genes_fdr = 0.05
        Int hs_genes_num_genes = 5000
        Float hs_modules_fdr = 0.05
        #general parameters
        Int cpu = 24
        String memory = "128G"
        Int extra_disk_space = 32
        String docker = "mparikhbroad/hotspot:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call run_hotspot {
        input:
            output_dir = output_directory_stripped,
            anndata_file = anndata_file,
            hotspot_model = hotspot_model,
            n_neighbors = n_neighbors,
            min_gene_threshold = min_gene_threshold,
            hs_genes_fdr = hs_genes_fdr,
            hs_genes_num_genes = hs_genes_num_genes,
            hs_modules_fdr = hs_modules_fdr,
            cpu=cpu,
            memory=memory,
            extra_disk_space = extra_disk_space,
            docker=docker,
            preemptible=preemptible
    }

    output {
        File hotspot_object = run_hotspot.hotspot_object
    }
}

task run_hotspot {

    input {
        String output_dir
        File anndata_file
        String hotspot_model
        Int n_neighbors
        Int min_gene_threshold
        Float hs_genes_fdr
        Int hs_genes_num_genes
        Float hs_modules_fdr
        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        mkdir -p outputs

        python <<CODE
        import os
        import scanpy as sc
        import hotspot

        jobs = ~{cpu} * 2

        adata = sc.read_h5ad("~{anndata_file}")

        hs = hotspot.Hotspot(
            adata,
            layer_key="counts",
            model='~{hotspot_model}',
            latent_obsm_key="X_pca",
            umi_counts_obs_key="total_counts"
        )

        hs.create_knn_graph(weighted_graph=False, n_neighbors=~{n_neighbors})

        hs_results = hs.compute_autocorrelations(jobs=jobs)
        hs_results.to_csv('outputs/hs_results.csv')

        hs_genes = hs_results.loc[hs_results.FDR < ~{hs_genes_fdr}].head(~{hs_genes_num_genes}).index # Select genes
        local_correlations = hs.compute_local_correlations(hs_genes, jobs=jobs) # jobs for parallelization
        local_correlations.to_csv('outputs/local_correlations.csv')

        modules = hs.create_modules(
            min_gene_threshold=~{min_gene_threshold}, core_only=True, fdr_threshold=~{hs_modules_fdr}
        )

        modules.to_csv('outputs/modules.csv')

        hs.plot_local_correlations(vmin=-10,vmax=10)
        plt.savefig('outputs/local_correlations_plot.png', bbox_inches='tight')

        module_scores = hs.calculate_module_scores()
        module_scores.to_csv('outputs/module_scores.csv')

        with open('outputs/hotspot_object.pickle', 'wb') as f:
            pickle.dump(hs, f, pickle.HIGHEST_PROTOCOL)

        CODE

        gsutil -m rsync -r outputs ~{output_dir}
    >>>

    output {
        File hotspot_object = 'outputs/hotspot_object.pickle'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(anndata_file, "GB")*4) + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}