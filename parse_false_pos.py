import sys

def parse_file(filename, skip_columns):
    result = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            extra_values = parts[skip_columns:]
            result.append(set(extra_values))
    return result

def compare_lists(list1, list2):
    total_false_positives = 0
    # Iterate through both lists simultaneously using zip
    for set1, set2 in zip(list1, list2):
        # Calculate false positives: elements in set2 that are NOT in set1
        # Set difference (set2 - set1) gives exactly these elements
        false_positives_at_index = set2 - set1
        # Add the count of false positives for this index to the total
        total_false_positives += len(false_positives_at_index)

    # Calculate the average false positives per query (index)
    average_false_positives = total_false_positives / len(list1)

    return average_false_positives


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_files.py <ground_truth> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    list1 = parse_file(file1, skip_columns=4)
    list2 = parse_file(file2, skip_columns=2)

    mismatches = compare_lists(list1, list2)

    print(f"Average # of false positives per query: {mismatches}")
