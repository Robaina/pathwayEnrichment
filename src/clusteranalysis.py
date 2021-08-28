from sklearn.metrics import silhouette_samples
import pandas as pd


def computeGeneSilhouettes(clusters, data):
    """
    Compute gene silhouettes for each gene in clusters
    """
    genes_in_clusters, cluster_labels = [], []
    for cluster_id, cluster in clusters.items():
        genes_in_clusters.extend(cluster)
        cluster_labels.extend([cluster_id for _ in range(len(cluster))])

    X = data.loc[genes_in_clusters, :].values
    sil_values = silhouette_samples(X, cluster_labels)

    return {gene_id: sil_values[i] for i, gene_id in enumerate(genes_in_clusters)}


def rankGenesWithinClusters(clusters, data):
    """
    Rank genes within each cluster based on their silhouette
    """

    gene_sil = computeGeneSilhouettes(clusters, data)

    ranked_clusters = {}
    for cluster_id, cluster in clusters.items():
        sil_dict = {gene_id: gene_sil[gene_id] for gene_id in cluster}
        ranked_dict = dict(
            sorted(sil_dict.items(), key=lambda item: item[1], reverse=True)
        )
        ranked_clusters[cluster_id] = ranked_dict

    return ranked_clusters


def writeExcelOfClusterGenes(clusters, out_path, gene_info,
                             gene_pathways, gene_systems):
    """
    Write excel file with gene product and KEGG pathways for genes in each cluster.
    """
    silhouette = type(clusters[list(clusters.keys())[0]]) == dict
    writer = pd.ExcelWriter(out_path, engine='xlsxwriter')

    for cluster_id, cluster in clusters.items():
        cluster_info = {}
        for gene_id in cluster:
            if gene_id in gene_info.keys():
                if silhouette:
                    cluster_info[gene_id] = {
                        'Gene silhouette': cluster[gene_id],
                        'Gene product': gene_info[gene_id],
                        'KEGG system': ', '.join(gene_systems[gene_id]),
                        'KEGG subsystem': ', '.join(gene_pathways[gene_id]),
                    }
                else:
                    cluster_info[gene_id] = {
                        'Gene product': gene_info[gene_id],
                        'KEGG system': ', '.join(gene_systems[gene_id]),
                        'KEGG subsystem': ', '.join(gene_pathways[gene_id]),
                    }

        pd.DataFrame(cluster_info).transpose().to_excel(writer, sheet_name=f'Cluster {cluster_id}')
    writer.save()
