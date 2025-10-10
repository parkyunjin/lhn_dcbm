import os
import glob
import pandas as pd

base_dir = '/result/'

results = []
bad_files = []

for seed in range(123, 223):
    print(f"Processing seed: {seed}")
    for graph_num in range(1, 8):
        folder = os.path.join(
            base_dir,
            f's_{seed}',
            f'results_percentile_n{graph_num}'
        )
        if not os.path.isdir(folder):
            print(f"Warning: missing folder {folder}")
            continue

        # find all eval_*.txt files
        pattern = os.path.join(folder, 'eval_*_permloss.txt')
        for filepath in glob.glob(pattern):
            fname = os.path.basename(filepath)
            core = fname[len('eval_'):-len('_permloss.txt')]
            # split into metric and algorithm
            try:
                metric, algorithm = core.rsplit('_', 1)
            except ValueError:
                metric, algorithm = core, ''


            df = pd.read_csv(filepath, sep='\s+', header=None, names=['file', 'value'])
            if df.empty:
                continue

            first_val = df['value'].iloc[0]
            min_val   = df['value'].iloc[1:].min()
            diff      = first_val - min_val
            other_values = df['value'].iloc[1:]
            diff_all = other_values - first_val
            mean_difference = diff_all.mean()
            sd_difference = diff_all.std()
            is_neg = diff_all<0
            count_neg = is_neg.sum()
            ratio_neg = is_neg.mean()

            results.append({
                'seed':       seed,
                'graph':      f'n{graph_num}',
                'metric':     metric,
                'algorithm':  algorithm,
                'baseline': first_val,
                'minval': min_val,
                'difference': diff,
                'diff_all_mean': mean_difference,
                'diff_all_sd': sd_difference,
                'count_neg': count_neg,
                'ratio_neg': ratio_neg
            })

# assemble and save
out_df = pd.DataFrame(results, columns=['seed','graph','metric','algorithm','baseline', 'minval', 'difference','diff_all_mean','diff_all_sd','count_neg','ratio_neg'])
out_df.to_csv('permloss_100.csv', index=False)

