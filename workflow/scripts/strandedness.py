import argparse

def determine_strandedness(file_path, cutoff1, cutoff2):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract the relevant values from the last two lines
    value1 = float(lines[-2].split(': ')[1])
    value2 = float(lines[-1].split(': ')[1])

    # Initialize the results dictionary
    results = {
        "default": "undetermined",
        "picard": "NONE",
        "featurecounts": 0,
        "biokit": 0
    }

    # Determine strandedness based on the conditions
    if value1 > cutoff1 and value2 > cutoff1:
        results["default"] = "unstranded"
    elif value1 > cutoff2 and value2 < (1 - cutoff2):
        results["default"] = "stranded"
        results["picard"] = "FIRST_READ_TRANSCRIPTION_STRAND"
        results["featurecounts"] = 1
        results["biokit"] = 1
    elif value2 > cutoff2 and value1 < (1 - cutoff2):
        results["default"] = "stranded"
        results["picard"] = "SECOND_READ_TRANSCRIPTION_STRAND"
        results["featurecounts"] = 2
        results["biokit"] = 2
 
    return results

def main():
    parser = argparse.ArgumentParser(description='Determine strandedness from a text file.')
    parser.add_argument('file_path', type=str, help='Path to the text file.')
    parser.add_argument('--cutoff1', type=float, help='Cutoff value for unstranded determination.', default=0.4)
    parser.add_argument('--cutoff2', type=float, help='Cutoff value for stranded determination.', default=0.9)

    args = parser.parse_args()

    results = determine_strandedness(args.file_path, args.cutoff1, args.cutoff2)

    # Print all results
    print(f"default={results['default']}")
    print(f"picard={results['picard']}")
    print(f"featurecounts={results['featurecounts']}")
    print(f"biokit={results['biokit']}")

if __name__ == "__main__":
    main()