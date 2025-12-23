# A streamlined way to update the sample config.yaml file needed for Snakemake pipeline runs
# Useful for folks who don't want to use VIM or a remote VSCode session to edit the file

# Evan Hackstadt
# Whittall Lab Dudleya Genomics
# December 2025



# TODO: REDO HOW WE WRITE TO CONFIG FILE!!! CONFIG FILE HAS NEW FORMAT!!!

# POTENTIALLY USEFUL OLD SAMPLE NAME PARSING FUNCTION:
'''
def shortname(sample):
    # code_LP_xxx_Du-xxx
    substrings = sample.split("_")
    LP = DU = ""
    code = substrings[0]
    code += substrings[1] if substrings[0] == 'CY' else ''
    for i, s in enumerate(substrings):
        if s == 'LP':
            LP = substrings[i+1]
        if 'Du' in s:
            DU = s.replace('Du-', '')
    return f"{code}{DU}"
'''


import sys
import os
import argparse

# --- CLI Args ---
parser = argparse.ArgumentParser()
parser.add_argument('data_dir', help='Path to the directory containing input data (raw reads) of interest. By default, adds all samples to a new config file.', type=str)
parser.add_argument('config_dir', help='Path to the directory where the config.yaml file should be saved.', type=str)
parser.add_argument('-n', '--n_samples', help='(optional) integer value --> add the first n samples to config file', type=int)
parser.add_argument('-c', '--custom_samples', help='(optional) allows you to enter custom samples for the config file', action='store_true')
parser.add_argument('-a', '--append_to_file', help='(optional) script will add samples to an existing "config.yaml" file in config_dir', action='store_true')
parser.add_argument('-q', '--quiet', help='(optional) script will not print modified config file contents after writing', action='store_true')
args = parser.parse_args()

# process args
data_dir = os.path.abspath(args.data_dir)   # ensure we have absolute paths
config_dir = os.path.abspath(args.config_dir)
n_samples = args.n_samples if args.n_samples else None
data_files = [f for f in os.listdir(data_dir) 
              if os.path.isfile(os.path.join(data_dir, f))]
data_files.sort()

if not os.path.isdir(data_dir):
    raise(ValueError, f"{data_dir} is not a directory.")
if not os.path.isdir(config_dir):
    raise(ValueError, f"{config_dir} is not a directory.")
if args.n_samples is not None and args.custom_samples is True:
    raise(ValueError, "Cannot use both --n_samples and --custom_samples flags. Choose one or neither.")
if len(data_files) == 0:
    raise(ValueError, f"{data_dir} does not contain any files.")

# parse filenames into samples (filenames up until _R1 or _R2)
data_samples = []
for f in data_files:
    if '_R1' in f:
        substrings = f.split('_R1')
        data_samples.append(substrings[0])  # left of the split (filename up until _R1)
    elif '_R2' in f:
        substrings = f.split('_R2')
        data_samples.append(substrings[0])  # left of the split (filename up until _R2)
data_samples = list(set(data_samples))    # remove duplicates
data_samples.sort()

print("\n----ARGS----")

print(f"Path to data directory: {data_dir}")

print(f"Found {len(data_files)} files in directory: ")
print(data_files) if not args.quiet else print("(message supressed)")

print(f"Parsed into {len(data_samples)} samples: ")
print(data_samples) if not args.quiet else print("(message supressed)")


# --- Select samples based on args provided ---
print("\n----SAMPLE SELECTION----")
selected_samples = []

if args.custom_samples:         # custom_samples
    print("Custom samples requested.")
    print("Please enter desired samples (not filenames) one-by-one, matching the options listed above.")
    print("Type 'save' to finish selection. Type 'quit' to cancel and exit.")
    print(">>>>>>>>")
    user_input = input()
    print("<<<<<<<<")
    while user_input != 'quit' and user_input != 'save':
        if user_input in data_samples and user_input not in selected_samples:
            selected_samples.append(user_input)
            print(f"Added {user_input}")
            print("All selected samples: ", selected_samples)
        else:
            print(f"Error: '{user_input}' is not a valid sample, or you've already selected it.")
            print(f"Valid samples: {[x for x in data_samples if x not in selected_samples]}") if not args.quiet else print()
        print("Enter another sample, or 'save' to finish, or 'quit' to exit.")
        print(">>>>>>>>")
        user_input = input()
        print("<<<<<<<<")
    
    if user_input == 'quit':
        print("Cancelling script. No config file written.")
        sys.exit()
    if user_input == 'save':
        print("--------")
        print("Custom sample entry complete.")
        print("Final custom samples:")
        print(selected_samples)
        
elif args.n_samples is not None:      # n_samples
    print("Specific number of samples requested.")
    print(f"Selecting the first {n_samples} from the directory.")
    if n_samples <= len(data_samples):
        selected_samples = data_samples[:n_samples]
    else:
        print(f"n={n_samples} provided, which is > number of samples present. ",
              f"Selecting all {len(data_samples)} samples instead.")
        selected_samples = data_samples
    
else:               # (default) all samples
    print("All samples will be added to config file.")
    selected_samples = data_samples

print(f"Selected {len(selected_samples)} samples to write: ", selected_samples)


# --- Write to the config file ---
print("\n----WRITING TO FILE----")
config_path = os.path.join(config_dir, "config.yaml")
print("Config file path: ", config_path)

if args.append_to_file:
    mode = "a"
    print("Append requested. Selected samples will be added to end of existing config file.")
else:
    mode = "w"
    print("Mode is write. Overwriting previous config file.")

with open(config_path, mode) as f:
    print(f"...Writing to {config_path} ...")
    
    # if we are writing a new file, add the path and samples header
    if mode == "w":
        f.write(f"# {config_path}")
        f.write(f"\nsamples_path: {data_dir}")
        f.write("\nsamples:")
    
    # now write the selected samples line-by-line
    for sample in selected_samples:
        # ensure sample is not already in file (skip duplicates)
        with open(config_path, 'r') as f2:
            contents = f2.read()
        if sample in contents:
            print(f">>>>WARNING: {sample} is already in config file. Skipping...")
        else:
            f.write(f"\n- {sample}")

print("Finished writing!")

if not args.quiet:
    print(f"Resulting contents of {config_path}:\n")
    with open(config_path, 'r') as f:
        print(f.read(), "\n")