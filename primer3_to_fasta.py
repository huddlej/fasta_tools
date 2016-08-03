import argparse
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert primer3 output into FASTA format")
    parser.add_argument("primer3")
    parser.add_argument("fasta")
    parser.add_argument("--table", action="store_true", help="output sequence in tab-delimited format")
    args = parser.parse_args()

    with open(args.primer3, "r") as fh:
        with open(args.fasta, "w") as oh:
            if args.table:
                oh.write("%s\n" % "\t".join(("site_id", "primer_orientation", "primer_number", "primer_sequence")))

            for line in fh:
                pieces = [piece for piece in line.strip().split("=") if piece != ""]

                if len(pieces) > 0:
                    if pieces[0] == "SEQUENCE_ID":
                        sequence_id = pieces[1]
                    else:
                        match = re.search("PRIMER_(LEFT|RIGHT)_(\d+)_SEQUENCE", pieces[0])
                        if match is not None:
                            primer_orientation = match.groups()[0]
                            primer_number = match.groups()[1]

                            if args.table:
                                oh.write("%s\n" % "\t".join((sequence_id, primer_orientation, primer_number, pieces[1])))
                            else:
                                oh.write(">%s_%s_%s\n" % (sequence_id, primer_orientation, primer_number))
                                oh.write("%s\n" % pieces[1])
