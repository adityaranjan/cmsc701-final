import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


args = sys.argv[1:]
file_path = args[0]

# Read the CSV data
df = pd.read_csv(file_path)

print(df)

# Set the first column as index
df.set_index('w\\k', inplace=True)

mode = "lexmin" if "lexmin" in file_path else "hashmin"

# Create the heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(df, annot=True, cmap='viridis', fmt=".3f")
plt.title(file_path[file_path.rfind('/') + 1:][:-4])
plt.xlabel('k')
plt.ylabel('w')
plt.savefig("{}_{}.pdf".format(file_path[file_path.rfind('/') + 1:][:-4], mode))
