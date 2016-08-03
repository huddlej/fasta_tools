import argparse
import csv
import pprint


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", help="GFF3 file of clone end mappings from NCBI's Clone Database")
    parser.add_argument("accession_to_chromosome_mapping", help="Mapping of NCBI accessions for assembly contigs to chromosome names (e.g., ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/chr_accessions_GRCh38)")
    parser.add_argument("--prepend_chr", action="store_true")
    args = parser.parse_args()

    chromosomes_by_accession = {}
    with open(args.accession_to_chromosome_mapping, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if row[0].startswith("#"):
                continue

            chromosomes_by_accession[row[1]] = row[0]

    # KEYS_TO_EMIT = (
    #     "chromosome",
    #     "start",
    #     "end",
    #     "Name",
    #     "score",
    #     "strand",
    #     "thickStart",
    #     "thickEnd",
    #     "itemRgb",
    #     "blockCount",
    #     "blockSizes",
    #     "blockStarts"
    # )

    KEYS_TO_EMIT = (
        "chromosome",
        "start",
        "end",
        "Name",
        "score",
        "strand",
    )
    records_by_id = {}

    with open(args.gff, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")

        for row in reader:
            if len(row) == 0 or row[0].startswith("#") or row[0].strip() == "":
                continue

            attributes = dict([item.split("=") for item in row[8].split(";")])

            # Clean up the chromosome name from accession to standard "chr1"
            # style.
            attributes["chromosome"] = [item for item in row[0].split("|") if item != "ref"][0]
            attributes["chromosome"] = chromosomes_by_accession.get(attributes["chromosome"], attributes["chromosome"])
            if args.prepend_chr:
                attributes["chromosome"] = "chr%s" % attributes["chromosome"]

            attributes["start"] = row[3]
            attributes["end"] = row[4]
            attributes["score"] = row[5].replace(".", "0")
            attributes["strand"] = row[6]

            # If there isn't a Parent attribute, this is the parent record.
            if not "Parent" in attributes:
                #attributes["children"] = []
                #records_by_id[attributes["ID"]] = attributes

                print "\t".join([str(attributes.get(key)) for key in KEYS_TO_EMIT])
            # else:
            #     # If there is a Parent attribute, update the parent record to
            #     # include this clone end by end type (e.g, "clone_insert_start",
            #     # etc.).
            #     records_by_id[attributes["Parent"]]["children"].append(attributes)

            #     # Check the parent record for both start and end. If both exist,
            #     # emit the record and remove it from records to consider.
            #     if len(records_by_id[attributes["Parent"]]["children"]) == 2:
            #         parent = records_by_id[attributes["Parent"]]
            #         parent["thickStart"] = parent["start"]
            #         parent["thickEnd"] = parent["end"]
            #         parent["itemRgb"] = "0,0,0"

            #         # Calculate block sizes and starts then sort them by start position.
            #         block_sizes = [int(child["end"]) - int(child["start"]) for child in parent["children"]]
            #         block_starts = [int(child["start"]) - int(parent["start"]) for child in parent["children"]]

            #         # Resize starts when child start is placed before parent
            #         # start.
            #         block_starts = [start > 0 and start or 0 for start in block_starts]

            #         # Zip sizes and starts together into one tuple per clone end
            #         # and sort tuples by start position. Each entry in
            #         # ``blocks`` is a clone end.
            #         blocks = zip(block_sizes, block_starts)
            #         blocks = sorted(blocks, key=lambda x: x[1])

            #         #pprint.pprint(blocks)

            #         # Resize overlapping blocks by setting the first block's
            #         # size to the start of the second block.
            #         blocks_changed = False
            #         if len(blocks) > 1 and blocks[0][0] > blocks[1][1]:
            #             blocks[0] = (blocks[1][1], blocks[0][1])
            #             blocks_changed = True

            #             # If resized blocks are zero-sized, omit empty blocks or
            #             # this record altogether.
            #             blocks = [block for block in blocks if block[0] > 0]

            #         # Continue if we eliminated all blocks.
            #         if len(blocks) == 0:
            #             continue

            #         # Assign block count after filtering blocks.
            #         parent["blockCount"] = len(blocks)

            #         # Unzip block tuples into separate BED columns again. Each
            #         # entry in ``blocks`` is a BED column.
            #         blocks = zip(*blocks)

            #         # Assign block sizes and starts to parent record.
            #         parent["blockSizes"], parent["blockStarts"] = [",".join(map(str, block)) for block in blocks]

            #         print "\t".join([str(parent.get(key)) for key in KEYS_TO_EMIT])

            #         del records_by_id[attributes["Parent"]]

