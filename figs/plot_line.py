import os
import pandas as pd
import matplotlib.pyplot as plt

base_path = "/Users/adi/School/CMSC_701/cmsc701-final/results"
orderings = ["lexmin", "hashmin"]
mode = "delta"
file_name_rate = "false_positive_rate"
file_name_count = "false_positives"

lists_fp_rate = [[0] * 5 for _ in range(len(orderings))]
lists_fp_count = [[0] * 5 for _ in range(len(orderings))]

i = 0
for ordering in orderings:
    ordering_path = os.path.join(base_path, ordering)

    if os.path.isdir(ordering_path):
        print(f"Navigating into ordering directory: {ordering_path}")
        for subdir in os.listdir(ordering_path):
            subdirectory_full_path = os.path.join(ordering_path, subdir)

            if subdir.startswith("default"):
                fp_rate_file_path = os.path.join(subdirectory_full_path, file_name_rate + ".csv")
                df_rate = pd.read_csv(fp_rate_file_path)
                df_rate.set_index('w\\k', inplace=True)
                lists_fp_rate[i][0] = df_rate.loc[11][0]

                fp_count_file_path = os.path.join(subdirectory_full_path, file_name_count + ".csv")
                df_count = pd.read_csv(fp_count_file_path)
                df_count.set_index('w\\k', inplace=True)
                lists_fp_count[i][0] = df_count.loc[11][0]

            if subdir.startswith(mode):
                fp_rate_file_path = os.path.join(subdirectory_full_path, file_name_rate + ".csv")
                df_rate = pd.read_csv(fp_rate_file_path)
                df_rate.set_index('w\\k', inplace=True)
                lists_fp_rate[i][int(subdir[-1])] = df_rate.loc[11][0]

                fp_count_file_path = os.path.join(subdirectory_full_path, file_name_count + ".csv")
                df_count = pd.read_csv(fp_count_file_path)
                df_count.set_index('w\\k', inplace=True)
                lists_fp_count[i][int(subdir[-1])] = df_count.loc[11][0]

    i += 1

fig, ax1 = plt.subplots()

ax1_color = 'tab:red'
ax1.set_xlabel("{} check count".format(mode))
ax1.set_ylabel("{}".format(file_name_count), color=ax1_color)

markers = ['o', '^']
linestyles = ['--', '-']

for i, ordering in enumerate(orderings):
    ax1.plot([j for j in range(len(lists_fp_count[i]))], lists_fp_count[i],
             label=f"{ordering} - {file_name_count}", color=ax1_color,
             marker=markers[i], linestyle=linestyles[i])
ax1.tick_params(axis='y', labelcolor=ax1_color)
ax1.grid(axis='y', linestyle='--', alpha=0.6) 


ax2 = ax1.twinx()

ax2_color = 'tab:blue'
ax2.set_ylabel("{}".format(file_name_rate), color=ax2_color)
for i, ordering in enumerate(orderings):
    ax2.plot([j for j in range(len(lists_fp_rate[i]))], lists_fp_rate[i],
             label=f"{ordering} - {file_name_rate}", color=ax2_color,
             marker=markers[i], linestyle=linestyles[i])
ax2.tick_params(axis='y', labelcolor=ax2_color)
ax2.grid(axis='y', linestyle='--', alpha=0.6) 


fig.tight_layout()
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.xticks([i for i in range(5)], [str(i) for i in range(5)])
plt.title(f"{mode} check count vs. {file_name_count} and {file_name_rate}")
plt.savefig(f"{mode}.pdf", bbox_inches='tight', dpi=300)
plt.tight_layout()
plt.show()