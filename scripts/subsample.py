import argparse
import os
import random
import sys

# create parser
parser = argparse.ArgumentParser()

# add arguments to the parser
parser.add_argument("--replicates", help="How many subsampled datasets to generate from the original. Replicates are of equal size.")
parser.add_argument("-n", "--subset", help="Size of each replicate. Files are selected randomly WITHOUT REPLACEMENT.")
parser.add_argument("--with_replacement", action='store_true', help="Select files with replacement such that each replicate contains no duplicate files, but multiple replicates may contain the same file.")
parser.add_argument("-d", "--dir", help="Source directory. Contains all files in dataset.")
parser.add_argument('-o', '--outfile', help="What to name all subsampled datasets. Format: OUTFILE-REPLICATE#.txt")

# parse the arguments
args = parser.parse_args()
replicates = int(args.replicates)
subset = int(args.subset)

# randomly select a subset of files
source_files = [f for f in os.listdir(args.dir) if not f.startswith('.')]  # ignore hidden files, which Snakemake generates automatically
if not source_files:
    print("ERROR: no files in source directory.")
    sys.exit(1)
if subset > len(source_files):
    print("WARNING: Your subset size meets or exceeds the number of files in the source directory.")
    if not args.with_replacement:
        print("You do not have enough data to generate a replicate at the current size. Please reduce the size of your subset.")
        sys.exit(1)
if (replicates * subset) > (len(source_files)-1) and not args.with_replacement:
    print("You do not have enough data to generate all replicates at the current size. Please reduce the number of replicates or the subset size, or allow replacement.")
    

if args.with_replacement:
    for i in range(replicates):
        subset_files = random.sample(source_files, subset)
        
        # add absolute path to filenames
        subset_files = [args.dir + "/" + f for f in subset_files]

        # write current bootstrap file list to text file
        align_list = open(args.outfile + str(i+1) + ".txt", "w")
        print(*subset_files, sep = "\n", file = align_list)
        align_list.close()
else:
    random.shuffle(source_files)
    subset_files = [source_files[x:x+subset] for x in range(0, len(source_files), subset)]
    subset_files = subset_files[:replicates]

    i = 1
    for s in subset_files:
        s = [args.dir + "/" + f for f in s]

        # write current bootstrap file list to text file
        align_list = open(args.outfile + str(i) + ".txt", "w")
        print(*s, sep = "\n", file = align_list)
        align_list.close()

        i += 1
