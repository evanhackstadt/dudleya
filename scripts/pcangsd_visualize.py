# Visualize PCAngsd covariance matrix as a 2D scatterplot of PCs (Principal Components)

# Evan Hackstadt
# Whittall Lab Dudleya Genomics
# 1 November 2025


import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2
import matplotlib.pyplot as plt
import seaborn as sns

# TODO: switch to seaborn plots


# --- CLI Args ---
parser = argparse.ArgumentParser()
parser.add_argument('dir', help='path to the directory containing pcangsd .cov and .info file, and optionally .selection', type=str)
parser.add_argument('-o', '--out', help='(optional) path to a dedicated output directory. if ommitted, dir is used.', type=str)
parser.add_argument('-x', '--x_axis_pc', help='(optional) integer value of the first of two PCs to plot', type=int)
parser.add_argument('-y', '--y_axis_pc', help='(optional) integer value of the second of two PCs to plot', type=int)
args = parser.parse_args()

# find files
dir_path = args.dir
out_path = args.out if args.out else args.dir
cov_path = info_path = selection_path = None
os.chdir(dir_path)
for f in os.listdir(dir_path):
    if f.endswith('.cov'):
        cov_path = os.path.join(dir_path, f)
    elif f.endswith('.info'):
        info_path = os.path.join(dir_path, f)
    elif f.endswith('.selection'):
        selection_path = os.path.join(dir_path, f)

# if we didn't find .cov and .info, raise an error
if not cov_path:
    raise(ValueError, "No .cov file found")
if not info_path:
    raise(ValueError, "No .info file found")

x_pc = args.x_axis_pc if args.x_axis_pc else 1
y_pc = args.y_axis_pc if args.y_axis_pc else 2


# --- Process Data ---
cov = np.loadtxt(cov_path)  # Reads estimated covariance matrix

info = np.loadtxt(info_path, dtype=str)    # Reads labels for our samples
pop_labels = []
point_labels = []

for sample in info:     # each sample = code_LP_xxx_Du-xxx
    substrings = sample.split('_')
    DU_num = substrings[-1].replace('Du-', '')
    if substrings[0] == 'CY':
        pop_labels.append('DUCY')
        point_labels.append(substrings[0] + '_' + substrings[1] + '_' + DU_num)
    else:
        pop_labels.append('DUSE')
        point_labels.append(substrings[0] + '_' + DU_num)

# process .cov into eigenvectors & eigenvalues
evals, evecs = np.linalg.eigh(cov)
# reverse order of evecs and evals so we have PC1, PC2, ... (descending variance)
evecs = evecs[:,::-1]
evals = evals[::-1]
evals_sum = np.sum(evals)

# put our data into a DataFrame
points = pd.DataFrame({   # each row is a point
    'pop': pop_labels,
    'point': point_labels,
})
# also store the variance explained by each PC
variance_explained = {}
# add info for each PC to the DataFrame and dict
for i in range(len(evals)):
    points[f'pc{i+1}'] = evecs[:,i]
    variance_explained[f'pc{i+1}'] = round(100*evals[i]/evals_sum, 2)


# --- PCA Plot ---
# one pop and one point at a time so we can label

unique_pops = points['pop'].unique()

fig, ax = plt.subplots(figsize=(8,6))

for pop in unique_pops:

    # TEMP manual color formatting
    if pop == 'DUCY':
        color = '#ff7f0e'
        marker = 'X'
    elif pop == 'DUSE':
        color = '#1f77b4'
        marker = 'o'
    else:
        color = None
        marker = None

    
    subset = points[points['pop'] == pop]
    ax.scatter(subset[f'pc{x_pc}'], subset[f'pc{y_pc}'], 
               label=pop, color=color, marker=marker)    #TODO: robust error handling
    for row in subset.itertuples():
        ax.annotate(row.point, xy=(getattr(row, f'pc{x_pc}'), getattr(row, f'pc{y_pc}')))

# decorate rest of plot
ax.set_xlabel(f"PC{x_pc} ({variance_explained[f'pc{x_pc}']}% of variance)",
              fontsize=18)
ax.set_ylabel(f"PC{y_pc} ({variance_explained[f'pc{y_pc}']}% of variance)",
              fontsize=18)
ax.set_title(f"{len(evals)}-Sample PCAngsd", fontsize=16)
ax.legend(fontsize=14)
png_path = os.path.join(out_path, f"pcangsd_pc{x_pc}_pc{y_pc}.png")
plt.savefig(png_path, dpi=400)
print("Plot saved to ", png_path)


# --- Extra Stuff ---

# Obtain p-values from PC-based selection scan, if provided
if selection_path:
    D = np.loadtxt(selection_path)    # Reads PC based selection statistics
    p = chi2.sf(D, 1)
    print("p-values from PC-based selection scan:\n", p)


# Save various info & metrics to a file
metrics_path = os.path.join(out_path, "metrics.txt")
if os.path.exists(metrics_path):
    os.remove(metrics_path)

with open(metrics_path, 'a') as f:
    f.write("——————\nScript Info\n——————\n")
    f.write("_arg_\t_received_\t_absolute_path_\n")
    f.write(f"dir:\t{args.dir}\t\t{os.path.abspath(args.dir)}\n")
    if args.out:
        f.write(f"out: {out_path}\t\t{os.path.abspath(args.out)}\n")
    f.write(f"x_pc:\t{x_pc}\n")
    f.write(f"y_pc:\t{y_pc}\n")
    f.write(f"Plot saved to: {png_path}\n")
    
    f.write("——————\nPopulation Info\n——————\n")
    f.write("_population_\t_sample_\n")
    for pop, point in zip(pop_labels, point_labels):
        f.write(f"{pop}\t\t{point}\n")

    f.write("——————\nPCA Summary\n——————\n")
    for i, val in enumerate(evals):
        f.write(f"PC{i+1} explains {100*evals[i]/evals_sum}% of variance\n")
    
    if selection_path:
        f.write("——————\nSelection Statistics\n——————\n")
        f.write(f"p-values from PC-based selection scan:\n{p}\n")

    # this will get really big with many samples. make optional?
    f.write("——————\nPCA Details\n——————\n")
    f.write("Eigenvalues:\t\tEigenvectors:\n")
    for val, vecs in zip(evals, evecs):
        f.write(f"{val}\t{vecs}\n")

    
print("Metrics file saved to ", metrics_path)
