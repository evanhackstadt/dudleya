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

'''
# --- OLD CLI Arg Handling ---
parser = argparse.ArgumentParser()
parser.add_argument('dir', help='path to the directory containing pcangsd .cov and .info file, and optionally .selection', type=str)
parser.add_argument('-o', '--out', help='(optional) path to a dedicated output directory. if ommitted, dir is used.', type=str)
parser.add_argument('-x', '--x_axis_pc', help='(optional) integer value of the first of two PCs to plot', type=int)
parser.add_argument('-y', '--y_axis_pc', help='(optional) integer value of the second of two PCs to plot', type=int)
args = parser.parse_args()

# find files
dir_path = os.path.abspath(args.dir)
out_path = os.path.abspath(args.out) if args.out else os.path.abspath(args.dir)
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
'''

# --- NEW Snakemake Arg Handling ---
# these are defined in the snakemake rule
cov_path = snakemake.input[0]
info_path = snakemake.input[1]
out_path = os.path.dirname(snakemake.output[0])     # output is .../pca.png; need the parent dir
x_pc = 1
y_pc = 2


# --- Process Data ---
cov = np.loadtxt(cov_path)  # Reads estimated covariance matrix
info = np.loadtxt(info_path, dtype=str)    # Reads sample names
species_labels = []
pop_labels = []
point_labels = []

# parse sample names assuming format: "POPCODE_LP_xxx_Du-xxx"
for sample in info:
    # split on _LP_ to isolate popcode
    substrings1 = sample.split('_LP_')
    popcode = substrings1[0]
    pop_labels.append(popcode)
    # extract Du# and LP# from split on '_'
    substrings2 = sample.split('_')
    for i, s in enumerate(substrings2):
        if 'Du-' in s:
            DU_label = s
        if s == 'LP':
            LP_num = substrings2[i+1]   # once we find 'LP', the next substring is the LP number
    point_labels.append(f"{DU_label}_LP_{LP_num}")
    # check if we have specified a non-setchelli (i.e. popcode = SPECIES_POPCODE)
    if '_' in popcode:
        popcode_substrings = popcode.split('_')
        species = popcode_substrings[0]
        if species in ['ABAB', 'ABBE', 'ABMU', 'CY']:
            if species == 'CY':
                species_labels.append('DUCY')
            else:
                species_labels.append(species)
    else:
        species_labels.append('DUSE')

# process .cov into eigenvectors & eigenvalues
evals, evecs = np.linalg.eigh(cov)
# reverse order of evecs and evals so we have PC1, PC2, ... (descending variance)
evecs = evecs[:,::-1]
evals = evals[::-1]
evals_sum = np.sum(evals)

# put our data into a DataFrame
points = pd.DataFrame({   # each row is a point
    'species': species_labels,
    'pop': pop_labels,
    'point': point_labels
})

# also store the variance explained by each PC
variance_explained = {}
# add info for each PC to the DataFrame and dict
for i in range(len(evals)):
    points[f'pc{i+1}'] = evecs[:,i]
    variance_explained[f'pc{i+1}'] = round(100*evals[i]/evals_sum, 2)


# --- PCA Plot ---
# one pop and one point at a time so we can label

'''
unique_pops = points['pop'].unique()

fig, ax = plt.subplots(figsize=(8,6))

# OLD MATPLOTLIB PLOTTING
for pop in unique_pops:
    subset = points[points['pop'] == pop]
    ax.scatter(subset[f'pc{x_pc}'], subset[f'pc{y_pc}'], label=pop)    #TODO: robust error handling
    # label points
    # for row in subset.itertuples():
    #     ax.annotate(row.point, xy=(getattr(row, f'pc{x_pc}'), getattr(row, f'pc{y_pc}')))

# decorate rest of plot
ax.set_xlabel(f"PC{x_pc} ({variance_explained[f'pc{x_pc}']}% of variance)",
              fontsize=18)
ax.set_ylabel(f"PC{y_pc} ({variance_explained[f'pc{y_pc}']}% of variance)",
              fontsize=18)
ax.set_title(f"{len(evals)}-Sample PCAngsd", fontsize=16)
ax.legend(fontsize=14)
'''

# NEW SEABORN PLOTTING
ax = sns.scatterplot(data=points, x=f"pc{x_pc}", y=f"pc{y_pc}",
                hue="species", style="species")                 # optional: legend=False
ax.legend(fontsize=10)
# sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
ax.set_title(f"{len(evals)}-Sample PCAngsd", fontsize=14)


png_path = os.path.join(out_path, f"pca.png")   # consider: if allowing modular PCs, include in filename
plt.savefig(png_path, dpi=600)
print("Plot saved to ", png_path)


# --- Extra Stuff ---

'''
# Obtain p-values from PC-based selection scan, if provided
if selection_path:
    D = np.loadtxt(selection_path)    # Reads PC based selection statistics
    p = chi2.sf(D, 1)
    print("p-values from PC-based selection scan:\n", p)
'''


# Save various info & metrics to a file
metrics_path = os.path.join(out_path, "metrics.txt")

with open(metrics_path, 'w') as f:
    # f.write("——————\nScript Info\n——————\n")
    
    f.write("——————\nPopulation Info\n——————\n")
    f.write("_species_\t_population_\t_sample_\n")
    for species, pop, point in zip(species_labels, pop_labels, point_labels):
        f.write(f"{species}\t{pop}\t\t{point}\n")

    f.write("——————\nPCA Summary\n——————\n")
    for i, val in enumerate(evals):
        f.write(f"PC{i+1} explains {100*evals[i]/evals_sum}% of variance\n")

    # this will get really big with many samples. make optional?
    f.write("——————\nPCA Details\n——————\n")
    f.write("Eigenvalues:\t\tEigenvectors:\n")
    for val, vecs in zip(evals, evecs):
        f.write(f"{val}\t{vecs}\n")

    
print("Metrics file saved to ", metrics_path)
