import argparse
import csv
import networkx as nx
import pybedtools


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("calls", help="BED file of calls from various callsets")
    parser.add_argument("--reciprocal_overlap", type=float, default=0.5, help="proportion of reciprocal overlap required to consider two calls the same")
    args = parser.parse_args()

    # Load all calls from the given BED file using only the first three
    # columns. Loading directly into a BedTool from the filename is not
    # guaranteed to work if the input BED file has the same number of columns as
    # a standard bedN format (e.g., bed12) and any of those columns contain
    # values of a type that differs from the standard.
    with open(args.calls, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        calls = pybedtools.BedTool([row[:3] for row in reader])

    # Intersect the given calls with themselves. Self-self overlaps will be
    # reported along with any other overlaps matching the given reciprocal
    # overlap proportion.
    intersected_calls = calls.intersect(b=calls, f=args.reciprocal_overlap, r=True, wao=True)

    # Create a graph connecting all calls that share a reciprocal overlap.
    graph = nx.Graph()
    for call in intersected_calls:
        left = tuple(call[:3])
        right = tuple(call[3:6])
        graph.add_edge(left, right)

    # Find all connected subgraphs such that calls with only self-self overlaps
    # will exist as singleton graphs while all nodes that overlap each other
    # directly or transitively will be clustered in the same graph.
    for subgraph in nx.connected_component_subgraphs(graph):
        # Identify the node in this subgraph that covers the most genomic space
        # (i.e., the largest call). This node represents all other nodes that
        # overlap in the same region.
        print "\t".join(max(subgraph.nodes(), key=lambda node: int(node[2]) - int(node[1])))
