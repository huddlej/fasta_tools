#!/usr/bin/env python
"""
Split a single FASTA file with multiple chromosome sequences into one file per
chromosome with the defline set to the name of the chromosome.
"""
import optparse
import os
import re
import sys


def split_fasta(handle, options):
    """
    Example regex patterns:

    >gi|253995373|gb|GL010027.1| Loxodonta africana unplaced genomic scaffold scaffold_0, whole genome shotgun sequence
    >gi|253588254|gb|AAGU03092498.1| Loxodonta africana cont3.92497, whole genome shotgun sequence
    name_regex = re.compile(r">.*?(chr\S+)")
    name_regex = re.compile(r">.*?(scaffold_\d+)")
    backup_name_regex = re.compile(r">.*?(cont\d+\.\d+)")
    """
    out_file = None
    previous_sequence_name = None
    name_regex = re.compile(options.pattern)
    for line in handle:
        # Look for a defline.
        if line.startswith(">"):
            match = name_regex.match(line)
            if not match:
                if not options.quiet:
                    print "Couldn't find a sequence name for defline: %s" % line

                continue

            sequence_name = "".join((options.prefix, match.groups()[0]))

            # If append mode is turned off or the sequence names differ and
            # there is an open file handle, close the open handle and prepare a
            # new one.
            if not options.append or previous_sequence_name != sequence_name:
                if not options.quiet:
                    print "Found sequence for %s" % sequence_name

                if not options.dryrun:
                    if out_file is not None:
                        out_file.close()

                    filename = "".join((re.sub("[^\w]+", options.delimiter, sequence_name), options.suffix))

                    if options.force_defline:
                        defline = ">%s\n" % filename
                    else:
                        defline = ">%s\n" % sequence_name

                    if options.output_directory:
                        filename = os.path.join(options.output_directory, filename)

                    out_file = open(filename, "w")
                    out_file.write(defline)

            previous_sequence_name = sequence_name
        else:
            if out_file is not None and not options.dryrun:
                out_file.write(line)

    if out_file:
        out_file.close()


if __name__ == "__main__":
    usage = "%prog [options] <input.fa or stdin>"
    parser = optparse.OptionParser(usage=usage)
    parser.set_defaults(pattern=r">(.+)", prefix="", suffix=".fasta", append=False, dryrun=False, delimiter="_")

    parser.add_option("-r", dest="pattern",
                      help="regular expression pattern to match each defline")
    parser.add_option("-d", dest="output_directory",
                      help="directory to place split files")
    parser.add_option("-p", dest="prefix",
                      help="string to prepend to each sequence name (and file name)")
    parser.add_option("-s", dest="suffix",
                      help="string to append to each output filename")
    parser.add_option("-a", dest="append", action="store_true",
                      help="append sequences from multiple defline matches to the same file")
    parser.add_option("-n", dest="dryrun", action="store_true",
                      help="print deflines of matching sequences but don't write out sequences themselves")
    parser.add_option("-q", dest="quiet", action="store_true",
                      help="don't print to screen")
    parser.add_option("-f", dest="force_defline", action="store_true",
                      help="force defline to be the same as filename")
    parser.add_option("--delimiter", dest="delimiter",
                      help="delimiter to replace special characters in output filenames")

    # The first positional argument is the name of the FASTA file to split.
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error("Specify the name of the FASTA file to split.")
        sys.exit(1)

    filename = args[0]
    if filename.lower() == "stdin":
        handle = sys.stdin
    else:
        handle = open(filename, "r")

    split_fasta(handle, options)
    handle.close()
