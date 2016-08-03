"""
Shuffle breakpoints for a given BED12 file into a given genomic space requiring
both ends to be shuffled into that space. This is distinct from bedtools shuffle
which places the entire event (not just the ends) into a given space.
"""
import argparse
from collections import OrderedDict
import copy
import csv
import pprint
import pybedtools
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("events", help="BED file of events whose breakpoints should be shuffled into the given genomic space")
    parser.add_argument("shuffle_space", help="BED file of coordinates representing the genomic space into which event breakpoints should be shuffled")
    parser.add_argument("genome", help="tab-delimited file with chromosome names and lengths (e.g., chromInfo.txt)")
    parser.add_argument("--max_attempts", help="maximum number of times to try to shuffle coordinates", type=int, default=1000)
    parser.add_argument("--iterations", help="number of times to repeat shuffling experiment for the given inputs", type=int, default=1)
    args = parser.parse_args()

    # Repeat the shuffling experiment as many times as the user has requested.
    for iteration in xrange(args.iterations):
        # Load events to shuffle.
        events = pybedtools.BedTool(args.events)

        # Load space to shuffle into.
        shuffle_space = pybedtools.BedTool(args.shuffle_space)

        # Load genomic coordinates.
        chromsizes = OrderedDict()
        with open(args.genome, "r") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                chromsizes[row[0]] = (0, int(row[1]))

        # Split BED12 events into BED6 coordinates and pair the results.
        split_events = events.bed6()
        split_events = split_events.set_chromsizes(chromsizes)

        # Calculate totals once.
        total_events = len(events)
        total_split_events = len(split_events)

        # While not all pairs have been shuffled or the max number of attempts has
        # been reached, shuffle all coordinates and find pairs that shuffle within
        # the given space.
        shuffled_pairs = {}
        attempts = 0

        while attempts < args.max_attempts and len(shuffled_pairs) < total_events:
            shuffled_events = split_events.shuffle(incl=args.shuffle_space)

            # For each shuffled pair, check whether the pairs map within the
            # required distance of each other on the same chromosome.
            for i in xrange(0, total_split_events, 2):
                # Skip this pair if it's already been shuffled.
                if i in shuffled_pairs:
                    continue

                first_event = shuffled_events[i]
                second_event = shuffled_events[i + 1]
                event_length = int(first_event.score)
                sys.stderr.write("Event length: %i\n" % event_length)

                # Handle two cases:
                # 1) distance from left breakpoint to putative right point maps within given genomic space
                # 2) distance from right breakpoint to putative left point maps within given genomic space

                # Case 1
                if first_event.start + event_length <= chromsizes[first_event.chrom][1]:
                    right_breakpoint_size = second_event.end - second_event.start
                    sys.stderr.write("Right breakpoint size: %i\n" % right_breakpoint_size)
                    first_event_right_breakpoint = copy.copy(first_event)
                    sys.stderr.write("Original right breakpoint end: %i\n" % first_event_right_breakpoint.end)
                    first_event_right_breakpoint.end = first_event.start + event_length
                    sys.stderr.write("Revised right breakpoint end: %i\n" % first_event_right_breakpoint.end)
                    first_event_right_breakpoint.start = first_event_right_breakpoint.end - right_breakpoint_size

                    # Intersect putative right breakpoint with genomic space.
                    intersection = shuffle_space.intersect([first_event_right_breakpoint], sorted=True)

                    if len(intersection) > 0:
                        shuffled_pairs[i] = (first_event.chrom, first_event.start, first_event_right_breakpoint.end, iteration + 1)
                        continue
                    else:
                        sys.stderr.write("First event's right breakpoint doesn't intersect segdups: %s and %s\n" % (first_event, first_event_right_breakpoint))
                else:
                    sys.stderr.write("First event's right breakpoint maps off the end: %s\n" % first_event)

                if second_event.end - event_length > 0:
                    left_breakpoint_size = first_event.end - first_event.start
                    sys.stderr.write("Left breakpoint size: %i\n" % left_breakpoint_size)
                    second_event_left_breakpoint = copy.copy(second_event)
                    sys.stderr.write("Original left breakpoint start: %i\n" % second_event_left_breakpoint.start)
                    second_event_left_breakpoint.start = second_event.end - event_length
                    sys.stderr.write("Revised left breakpoint start: %i\n" % second_event_left_breakpoint.start)
                    second_event_left_breakpoint.end = second_event_left_breakpoint.start + left_breakpoint_size

                    # Intersect putative right breakpoint with genomic space.
                    intersection = shuffle_space.intersect([second_event_left_breakpoint], sorted=True)

                    if len(intersection) > 0:
                        shuffled_pairs[i] = (second_event.chrom, second_event_left_breakpoint.start, second_event.end, iteration + 1)
                        continue
                    else:
                        sys.stderr.write("Second event's left breakpoint doesn't intersect segdups: %s and %s\n" % (second_event, second_event_left_breakpoint))
                else:
                    sys.stderr.write("Second event's left breakpoint maps off the end: %s\n" % second_event)

                if not i in shuffled_pairs:
                    sys.stderr.write("Couldn't shuffle events: %s\n" % split_events[i])

            attempts += 1

            if attempts % 20 == 0:
                sys.stderr.write("%i of %i pairs shuffled\n" % (len(shuffled_pairs), total_events))

        if attempts == args.max_attempts:
            sys.stderr.write("[WARNING] reached maximum attempts: %i\n" % args.max_attempts)

        for pair in shuffled_pairs.itervalues():
           print "\t".join(map(str, pair))

        # Clean up tmp files at the end of each iteration.
        pybedtools.cleanup()
