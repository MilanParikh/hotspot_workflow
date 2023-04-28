version 1.0

workflow cluster_robustness {
    input {
    	String output_directory
        File anndata
        Array[Float] resolutions
        Int num_iter
        Float subset_ratio = 0.8
        #general parameters
        Int cpu = 8
        String memory = "64G"
        String docker = "mparikhbroad/cluster_robustness:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    scatter(resolution in resolutions) {
        call run_cluster_robustness {
            input:
                anndata = anndata,
                resolution = resolution,
                subset_ratio = subset_ratio,
                num_iter = num_iter,
                cpu=cpu,
                memory=memory,
                docker=docker,
                preemptible=preemptible
        }
    }

    call compile_iterations {
        input:
            output_dir = output_directory_stripped,
            rand_index_score_files = run_cluster_robustness.rand_index_scores,
            resolutions = resolutions,
            num_iter = num_iter,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    output {
        File rand_index_score_plot_pdf = compile_iterations.rand_index_score_plot_pdf
        File rand_index_score_table_csv = compile_iterations.rand_index_score_table_csv
    }
}

task run_cluster_robustness {

    input {
        File anndata
        Float resolution
        Float subset_ratio
        Int num_iter
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command {
        set -e

        python << CODE
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import scanpy.external as sce
        import scipy as sp
        from sklearn.model_selection import train_test_split
        from sklearn.metrics.cluster import adjusted_rand_score
        from tqdm import tqdm
        import matplotlib.pyplot as plt

        adata = sc.read_h5ad('~{anndata}')
        sc.tl.leiden(adata, resolution=~{resolution})
        rand_score_arr = []
        for j in tqdm(range(0, ~{num_iter}), desc='iteration'):
            #get subset and save clustering solution from full dataset
            train, test = train_test_split(adata.obs, train_size=~{subset_ratio}, random_state=(j), stratify=adata.obs[['leiden']])
            temp = adata[train.index].copy()
            temp.obs['original_leiden'] = temp.obs['leiden']
            #rerun HVG, PCA, and leiden clustering on subset
            sc.pp.highly_variable_genes(temp, min_mean=0.0125, max_mean=3, min_disp=0.5)
            sc.tl.pca(temp, svd_solver='arpack')
            sc.pp.neighbors(temp, n_neighbors=10, n_pcs=40)
            sc.tl.leiden(temp, resolution=~{resolution})
            #calculate adjusted rand score between original cluster labels and new subset cluster labels
            rand_score_arr.append(adjusted_rand_score(temp.obs.original_leiden, temp.obs.leiden))
        np.save('~{resolution}.npy', rand_score_arr, allow_pickle=False)
        CODE
        
    }

    output {
        File rand_index_scores = "~{resolution}.npy"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(anndata, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}

task compile_iterations {

    input {
        String output_dir
        Array[File] rand_index_score_files
        Array[Float] resolutions
        Int num_iter
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command {
        set -e

        python << CODE
        import os
        import numpy as np
        import pandas as pd
        import scipy as sp
        import seaborn as sns
        import matplotlib.pyplot as plt

        list_of_files = ["${sep='","' rand_index_score_files}"]
        list_of_scorelists = []
        for score_file in list_of_files:
            temp = np.load(score_file)
            list_of_scorelists.append(temp)
        resolution_list = [~{sep=", " resolutions}]
        score_df = pd.DataFrame(list_of_scorelists, index=[str(i) for i in resolution_list], columns=[str(i) for i in range(0, ~{num_iter})])
        score_df.to_csv('rand_index_score_table.csv')

        score_df = score_df.T
        sns.set(rc={'figure.figsize':(10,7.5)})
        ax = sns.boxplot(data=score_df)
        ax.set(xlabel='Leiden Resolution', ylabel='Rand Index Score', title='Rand Index Scores per Leiden Resolution', ylim=(0,1))
        ax.figure.savefig('rand_index_score_plot.pdf')

        CODE
        
        gsutil -m cp rand_index_score_plot.pdf ~{output_dir}/
        gsutil -m cp rand_index_score_table.csv ~{output_dir}/
    }

    output {
        File rand_index_score_plot_pdf = 'rand_index_score_plot.pdf'
        File rand_index_score_table_csv = 'rand_index_score_table.csv'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(rand_index_score_files, "GB")*4) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}