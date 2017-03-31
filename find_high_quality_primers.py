import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sites", help="BED file of sites for which primers have been designed with site id in the name field (column 4)")
    parser.add_argument("primers", help="headered tab-delimited file of one or more primer pairs designed for each site in the format of site id, primer orientation, primer number, primer sequence")
    parser.add_argument("primer_alignments", help="header tab-delimited file with number of alignments per primer in the format of site id, primer number, primer orientation, number of alignments")
    parser.add_argument("high_quality_primers", help="tab-delimited file of high-quality primers")
    parser.add_argument("sites_without_primers", help="list of site ids without high-quality primers")
    args = parser.parse_args()

    # Load set of sites targeted for primer design.
    sites = pd.read_table(args.sites, header=None, names=("chromosome", "start", "end", "site_id"))
    site_ids = set(sites["site_id"].tolist())

    # Load primer names and sequences.
    primers = pd.read_table(args.primers)

    # Load alignment counts per primer.
    alignments = pd.read_table(args.primer_alignments)

    # Find all primers with distinct alignments.
    distinct_alignments = alignments.loc[alignments["alignments"] == 1].groupby(["site_id", "primer_number"])

    # Find all primer pairs (forward/reverse or left/right) for which both
    # primers in the pair have distinct alignments.
    distinct_alignment_counts = distinct_alignments.count()
    primer_pairs_with_distinct_alignments = distinct_alignment_counts.loc[distinct_alignment_counts["alignments"] == 2,].index

    # Take the first primer pair for each site. If primers were provided by
    # primer3, they should already be in order of highest to lowest score from
    # primer3's scoring system.
    primer_pairs_with_distinct_alignments_df = pd.DataFrame(primer_pairs_with_distinct_alignments.tolist(), columns=primer_pairs_with_distinct_alignments.names)
    high_quality_primers = primer_pairs_with_distinct_alignments_df.groupby("site_id").head(1)

    # Get sequences associated with high-quality primers, group sequences by
    # site id with both primers on the same line, and save.
    high_quality_primer_sequences = high_quality_primers.merge(primers, how="left", on=("site_id", "primer_number"))
    high_quality_primer_sequences_by_site = high_quality_primer_sequences.pivot("site_id", "primer_orientation", "primer_sequence")
    high_quality_primer_sequences_by_site.to_csv(args.high_quality_primers, sep="\t", index=True, header=True)

    # Find sites with high-quality primers.
    sites_with_high_quality_primers = set(high_quality_primer_sequences_by_site.index.tolist())

    # Find sites without high-quality primers.
    sites_without_high_quality_primers = list(site_ids - sites_with_high_quality_primers)

    with open(args.sites_without_primers, "w") as oh:
        for site in sites_without_high_quality_primers:
            oh.write("%s\n" % site)
