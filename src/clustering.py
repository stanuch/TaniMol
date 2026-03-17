import numpy as np
from scipy.spatial.distance import squareform
from rdkit.ML.Cluster import Butina

def cluster_similarity_matrix(similarity_matrix, threshold=0.6):
    """
    The Butina algorithm expects distances(Distance = 1 - Similarity). 
    Therefore, a similarity threshold of 0.6 means a distance threshold of 0.4.
    
    The algorithm also expects a flat, condensed distance matrix (only the upper 
    triangle without the diagonal) to avoid redundant calculations.
    
    Parameters
    ----------
    similarity_matrix : np.ndarray
        A square N x N Tanimoto similarity matrix.
    threshold : float
        The Tanimoto similarity threshold for two molecules to be considered
        "neighbors" in the same cluster. Typically 0.6 for ECFP4.
        
    Returns
    -------
    tuple[dict, list]
        - A dictionary where keys are the centroid molecule indices, and values
          are tuples of all molecule indices belonging to that cluster.
        - A list of singleton indices (molecules that do not belong to any cluster).
    """
    n_mols = similarity_matrix.shape[0]
    distance_matrix = np.clip(1 - similarity_matrix, 0, 1)  
    dists = [distance_matrix[i, j] for i in range(1, n_mols) for j in range(i)]
    clusters = Butina.ClusterData(dists, n_mols, 1-threshold, isDistData=True)
    
    cluster_dict = {}
    singletons = []
    for cluster in clusters:
        if len(cluster) == 1:
            singletons.append(cluster[0])
        else:
            cluster_dict[cluster[0]] = cluster
    
    return cluster_dict, singletons

def analyze_clusters(clusters, singletons):
    total_clustered_mols = sum(len(c) for c in clusters.values())
    biggest_cluster_size = max(len(c) for c in clusters.values()) if clusters else 0
    large_clusters = sum(1 for c in clusters.values() if len(c) > 50)
    
    print("\n--- Clustering Analysis ---")
    print(f"Total Clusters:           {len(clusters)}")
    print(f"Total Singletons:         {len(singletons)}")
    print(f"Molecules in Clusters:    {total_clustered_mols}")
    print(f"Biggest Cluster Size:     {biggest_cluster_size} molecules")
    print(f"Clusters > 50 molecules:  {large_clusters}")
    print("---------------------------\n")

    return {
        "num_clusters": len(clusters),
        "num_singletons": len(singletons),
        "total_clustered": total_clustered_mols,
        "max_size": biggest_cluster_size,
        "large_clusters": large_clusters
    }

