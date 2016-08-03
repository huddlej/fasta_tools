import argparse
from pybedtools import BedTool


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bed12")
    parser.add_argument("introns")
    args = parser.parse_args()

    bed = BedTool(args.bed12)
    introns = bed.introns()
    introns.remove_invalid().saveas(args.introns)
