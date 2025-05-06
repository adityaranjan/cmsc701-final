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
    mismatched_indices = []
    for idx, (set1, set2) in enumerate(zip(list1, list2)):
        if not set1.issubset(set2):
            mismatched_indices.append(idx)
    return mismatched_indices

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_files.py <file1> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    list1 = parse_file(file1, skip_columns=4)
    list2 = parse_file(file2, skip_columns=2)

    mismatches = compare_lists(list1, list2)

    print(len(mismatches))

    """
    if mismatches:
        print("Indices where file2 does NOT fully contain file1's extra values:")
        for idx in mismatches:
            print(f"Line {idx}")
    else:
        print("All lines in file2 fully contain the extra values from file1.")
    """
