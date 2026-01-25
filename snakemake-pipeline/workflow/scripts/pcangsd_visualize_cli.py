# Command line version of pcangsd_visualize
# (Main script is used by Snakemake pipeline; to call manually with custom args, use this one)

# Evan Hackstadt
# Whittall Lab Dudleya Genomics
# 24 January 2026


import os
import argparse
from pcangsd_visualize import plot_pca


# --- CLI Arg Handling ---
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

plot_pca(cov_path, info_path, out_path, x_pc=x_pc, y_pc=y_pc)