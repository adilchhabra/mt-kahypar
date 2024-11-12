import argparse
from collections import defaultdict

def convert_to_hmetis(input_file, output_file, cluster_file):
    # Initialize counters and storage for unique hypernodes
    num_hyperedges = 0
    unique_hypernodes = set()
    hyperedges = []

    # Read the input file and process each line for hyperedges
    with open(input_file, 'r') as infile:
        for line in infile:
            # Remove whitespace and split the line by commas
            pins = line.strip().split(',')
            # Store the processed hyperedge
            hyperedges.append(pins)
            # Update counters
            num_hyperedges += 1
            unique_hypernodes.update(pins)  # Add pins to the set of unique hypernodes

    # The number of unique hypernodes is the size of the set
    num_unique_hypernodes = len(unique_hypernodes)

    # Read the cluster file and gather statistics
    cluster_counts = defaultdict(int)
    with open(cluster_file, 'r') as cfile:
        for node_id, line in enumerate(cfile):
            cluster_id = line.strip()
            cluster_counts[cluster_id] += 1

    num_clusters = len(cluster_counts)

    # Write to the output file in hMETIS format
    with open(output_file, 'w') as outfile:
        # Write the first line: number of hyperedges and number of unique hypernodes
        outfile.write(f"{num_hyperedges} {num_unique_hypernodes}\n")
        # Write each hyperedge line with pins separated by a space
        for pins in hyperedges:
            outfile.write(" ".join(pins) + "\n")

    # Print cluster statistics
    print(f"Conversion complete! Output written to {output_file}")
    print(f"Cluster statistics:")
    print(f"Total clusters: {num_clusters}")
    print("Hypernodes per cluster:")
    for cluster_id, count in cluster_counts.items():
        print(f"Cluster {cluster_id}: {count} hypernodes")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Convert hypergraph edges format to hMETIS format with cluster statistics.")
    parser.add_argument("input_file", help="Path to the input file in hypergraph edges format")
    parser.add_argument("output_file", help="Path to the output file in hMETIS format")
    parser.add_argument("cluster_file", help="Path to the file with cluster IDs for each hypernode")

    # Parse the arguments
    args = parser.parse_args()

    # Run the conversion function with provided arguments
    convert_to_hmetis(args.input_file, args.output_file, args.cluster_file)
