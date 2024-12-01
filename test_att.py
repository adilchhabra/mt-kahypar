import re

def parse_line(line):
    """
    Parse a single line and extract relevant data.
    """
    match = re.match(r"Node\s+(\d+)\s+to cluster\s+(\d+)\s+has pi_mod gain\s+([\d\.\-e]+)\s+with from =\s+([\d\.\-e]+)\s+and to =\s+([\d\.\-e]+)", line)
    if match:
        node = int(match.group(1))
        cluster = int(match.group(2))
        pi_mod_gain = float(match.group(3))
        from_val = float(match.group(4))
        to_val = float(match.group(5))
        return node, cluster, pi_mod_gain, from_val, to_val
    return None

def check_discrepancies(non_att_data, att_sums, discrepancies):
    """
    Check discrepancies for a given set of ATT and non-ATT data.
    """
    for (node, cluster), (pi_mod_gain, from_val, to_val) in non_att_data.items():
        if (node, cluster) in att_sums:
            att_pi_mod_gain, att_from_val, att_to_val = att_sums[(node, cluster)]
            if not (
                    abs(att_pi_mod_gain - pi_mod_gain) < 1e-6 and
                    abs(att_from_val - from_val) < 1e-6 and
                    abs(att_to_val - to_val) < 1e-6
            ):
                discrepancies.append((node, cluster, att_sums[(node, cluster)], (pi_mod_gain, from_val, to_val)))

def check_file(file_path):
    """
    Read and validate the file based on the described logic.
    """
    non_att_data = {}
    att_sums = {}
    discrepancies = []

    current_node = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            is_att = line.startswith("ATT:")
            parsed_data = parse_line(line.replace("ATT:", "").strip())
            if parsed_data:
                node, cluster, pi_mod_gain, from_val, to_val = parsed_data

                if current_node is None:
                    current_node = node

                if node != current_node:
                    # New node starts, check discrepancies for the previous node
                    check_discrepancies(non_att_data, att_sums, discrepancies)
                    non_att_data = {}
                    att_sums = {}
                    current_node = node

                if is_att:
                    # Sum up values for ATT lines
                    if (node, cluster) not in att_sums:
                        att_sums[(node, cluster)] = [0.0, 0.0, 0.0]
                    att_sums[(node, cluster)][0] += pi_mod_gain
                    att_sums[(node, cluster)][1] += from_val
                    att_sums[(node, cluster)][2] += to_val
                else:
                    # Store non-ATT values
                    non_att_data[(node, cluster)] = (pi_mod_gain, from_val, to_val)

        # Check discrepancies for the last node
        check_discrepancies(non_att_data, att_sums, discrepancies)

    # Output discrepancies
    if discrepancies:
        print("Discrepancies found:")
        for node, cluster, att_totals, expected_values in discrepancies:
            print(f"Node {node}, Cluster {cluster}:")
            print(f"  ATT Totals = {att_totals}")
            print(f"  Expected = {expected_values}")
    else:
        print("All ATT lines match the non-ATT data.")

# Example usage
file_path = "build/see.txt"  # Replace with your actual file path
check_file(file_path)
