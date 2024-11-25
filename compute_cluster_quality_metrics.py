import sys
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

def read_clusters_from_file(file_path):
    """
    Reads cluster IDs from a file where each line corresponds to a node's cluster ID.

    Parameters:
        file_path (str): Path to the file containing cluster IDs.

    Returns:
        list: A list of cluster IDs.
    """
    try:
        with open(file_path, 'r') as file:
            return [int(line.strip()) for line in file]
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

def compute_ami_ari(predicted_clusters, ground_truth_clusters):
    """
    Compute Adjusted Mutual Information (AMI) and Adjusted Rand Index (ARI).

    Parameters:
        predicted_clusters (list): Cluster assignments for each node in the hypergraph as predicted by an algorithm.
        ground_truth_clusters (list): Ground truth cluster assignments for each node.

    Returns:
        dict: A dictionary containing AMI and ARI scores.
    """
    ami = adjusted_mutual_info_score(ground_truth_clusters, predicted_clusters)
    ari = adjusted_rand_score(ground_truth_clusters, predicted_clusters)
    return {"AMI": ami, "ARI": ari}

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <predicted_clusters_file> <ground_truth_clusters_file>")
        sys.exit(1)

    predicted_file = sys.argv[1]
    ground_truth_file = sys.argv[2]

    # Read clusters from files
    predicted_clusters = read_clusters_from_file(predicted_file)
    ground_truth_clusters = read_clusters_from_file(ground_truth_file)

    # Compute and print AMI and ARI
    scores = compute_ami_ari(predicted_clusters, ground_truth_clusters)
    print("Adjusted Mutual Information (AMI):", scores["AMI"])
    print("Adjusted Rand Index (ARI):", scores["ARI"])
