import numpy as np
import pandas as pd

def compute_distribution(file_path):
    counts = []
    sequences = []

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()

            if len(parts) == 2:
                sequence = parts[0]
                count = int(parts[1])
                counts.append(count)
                sequences.append(sequence)

    mean_count = np.mean(counts)
    std_count = np.std(counts)
    upper_cut = round(mean_count + (3*std_count))
    lower_cut = round(mean_count - (3*std_count))

    df = pd.DataFrame({'Sequence': sequences, 'Count': counts})
    filtered_df = df[(df['Count'] > upper_cut) | (df['Count'] < lower_cut)]

    return mean_count, std_count, upper_cut, lower_cut, filtered_df


file_path = 'merged_sequence_counts.txt'

mean, std, upper, lower, filtered_df = compute_distribution(file_path)

output_file = 'filtered_sequences.txt'
filtered_df.to_csv(output_file, sep='\t', index=False, header=True)

print(f"Mean Count: {mean}")
print(f"Standard Deviation: {std}")
print(f"Upper Limit: {upper}")
print(f"Lower Limit: {lower}")
print(f"List of Over or Underrepresented Amino Acid Sequences: {filtered_df}")