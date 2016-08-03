import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import operator
import sys

# https://pypi.python.org/pypi/matplotlib-venn
from matplotlib_venn import venn2, venn3

# https://pythonhosted.org/pybedtools
from pybedtools import BedTool
from pybedtools.contrib.venn_maker import cleaned_intersect


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs="+")
    parser.add_argument("--names", nargs="+")
    parser.add_argument("--output")
    parser.add_argument("--title", default="")
    parser.add_argument("--textless", action="store_true", help="print text from Venn circles to the terminal and not in the circles themselves")
    parser.add_argument("--lists", action="store_true", help="input files are simple lists ready for comparison instead of BED files")
    parser.add_argument("--fontsize", type=int, help="font size for Venn labels", default=6)
    args = parser.parse_args()

    if args.files is None or len(args.files) not in (2, 3):
        sys.stderr.write("Error: Venns can only be built with 2 or 3 input files.\n")
        parser.print_usage(file=sys.stderr)
        sys.exit(1)

    # Set font size globally to keep Venn labels from overlapping.
    label_font_size = 12
    sublabel_font_size = args.fontsize
    font = {'size': args.fontsize}
    mpl.rc('font', **font)

    if not args.lists:
        # Create a BED instance for each input BED file.
        beds = [BedTool(file_name) for file_name in args.files]

        # Calculate the cleaned intersection of all BED files using set logic
        # (bedtools subtraction instead of intersection).
        intersection = cleaned_intersect(beds)

        # Create Python sets of all resulting intervals per file from the
        # intersection.
        sets = []
        for i in xrange(len(intersection)):
            sets.append(set([str(interval).rstrip().replace("\t", "-") for interval in intersection[i]]))
    else:
        # Inputs are simple lists to be compared directly in matplotlib-venn.
        sets = []
        for filename in args.files:
            with open(filename, "r") as fh:
                sets.append(set([line.strip() for line in fh]))

    # Build the Venn diagram given the calculated sets and given set names.
    f = plt.figure(figsize=(8,6))

    if len(sets) == 2:
        v = venn2(sets, args.names)

        named_sets = {
            "10": sets[0] - sets[1],
            "01": sets[1] - sets[0],
            "11": sets[0].intersection(sets[1])
        }
    elif len(sets) == 3:
        v = venn3(sets, args.names)

        named_sets = {
            "100": sets[0] - sets[1] - sets[2],
            "010": sets[1] - sets[0] - sets[2],
            "001": sets[2] - sets[0] - sets[1],
            "110": sets[0].intersection(sets[1]) - sets[2],
            "101": sets[0].intersection(sets[2]) - sets[1],
            "011": sets[1].intersection(sets[2]) - sets[0],
            "111": sets[0].intersection(sets[1]).intersection(sets[2])
        }

    # Update labels inside circles to include mean and median length of events
    # in the count.
    if not args.lists:
        for name, named_set in named_sets.iteritems():
            set_lengths = [operator.sub(*map(int, sorted(i.split("-")[1:3], reverse=True))) for i in named_set]
            new_label_text = "%s\n(%i bp,\n%i bp)" % (
                len(named_set),
                np.ceil(np.mean(set_lengths)),
                np.ceil(np.median(set_lengths)),
             )
            current_label_text = v.get_label_by_id(name).get_text()

            if args.textless:
                name_pieces = []
                for i in xrange(len(name)):
                    if name[i] == "1":
                        name_pieces.append(args.names[i])

                print "\t".join((" and ".join(name_pieces), new_label_text.replace("\n", " ")))
                v.get_label_by_id(name).set_text("")
            else:
                v.get_label_by_id(name).set_text(new_label_text)

    # Set font size.
    for l in v.subset_labels:
        l.set_fontsize(sublabel_font_size)

    for l in v.set_labels:
        l.set_fontsize(label_font_size)

    for p in v.patches:
        p.set_edgecolor("black")

    plt.title(args.title, fontsize=14)
    plt.savefig(args.output)
