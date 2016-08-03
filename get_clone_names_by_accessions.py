import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re

Entrez.email = ""
MAX_RECORDS_TO_RETURN = 50000

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", help="text file with one clone accession (e.g., AC254230) per line")
    parser.add_argument("annotated_accessions", help="tab-delimited file with clone name, nucleotide id, CloneDB id, and accessions")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    with open(args.accessions, "r") as fh:
        accessions = [line.strip() for line in fh]

    logger.info("Searching %i accessions", len(accessions))

    # Search NCBI nucleotide database for record ids associated with accessions.
    handle = Entrez.esearch(db="nuccore", term=" OR ".join(accessions), retmax=MAX_RECORDS_TO_RETURN)
    clone_record = Entrez.read(handle)
    clone_ids = clone_record["IdList"]

    logger.info("Found %i matching nucleotide records", len(clone_ids))

    # Download complete nucleotide records for each id returned in the search.
    handle = Entrez.esummary(db="nuccore", id=",".join(clone_ids), retmax=MAX_RECORDS_TO_RETURN)
    nuccore_summary = Entrez.read(handle)

    # Bind nucleotide record ids with accessions.
    logger.info("Downloaded %i nucleotide summaries", len(nuccore_summary))
    # ABC12_46660700_B3
    # CH17-124M21
    accession_nuccore_and_clone_id = pd.DataFrame([{"nuccore_id": record["Gi"], "accession": record["Caption"], "clone_id": record["Title"].split("clone ")[1].split(",")[0]}
                                                   for record in nuccore_summary])
    accession_nuccore_and_clone_id.to_csv(args.annotated_accessions, sep="\t", index=False)

    # # Look up the link between nucleotide ids and clone database ids.
    # logger.info("Search link between nuccore and clone with ids: %s", clone_ids)
    # handle = Entrez.elink(dbfrom="nuccore", db="clone", id=",".join(clone_ids))
    # link_record = Entrez.read(handle)

    # # Bind nucleotide record ids with clone db ids.
    # logger.debug(str(link_record))
    # clone_ids_in_clonedb = [record["Id"] for record in link_record[0]["LinkSetDb"][0]["Link"]]
    # logger.info("Found %i links from nucleotide to clone database", len(clone_ids_in_clonedb))
    # clonedb_and_nuccore_id = pd.DataFrame([{"clonedb_id": record[0], "nuccore_id": record[1]}
    #                                        for record in zip(clone_ids_in_clonedb, link_record[0]["IdList"])])

    # # Download complete clone records for each clone db id.
    # handle = Entrez.esummary(db="clone", id=",".join(clone_ids_in_clonedb), retmax=MAX_RECORDS_TO_RETURN)
    # clone_summary = Entrez.read(handle)
    # logger.info("Downloaded %i clone summaries", len(clone_summary["DocumentSummarySet"]["DocumentSummary"]))

    # # Bind clone db ids with clone names.
    # clonedb_id_and_clone_name = pd.DataFrame([{"clonedb_id": record.attributes["uid"], "clone_name": record["ClName"]} for record in clone_summary["DocumentSummarySet"]["DocumentSummary"]])

    # # Merge clone name with clone db and nucleotide ids.
    # clone_name_and_ids = clonedb_id_and_clone_name.merge(clonedb_and_nuccore_id, on="clonedb_id", how="left")

    # # Merge clone name and ids with accessions.
    # clone_name_ids_and_accessions = accession_and_nuccore_id.astype(str).merge(clone_name_and_ids, on="nuccore_id", how="left")
    # clone_name_ids_and_accessions.sort_values("clone_name", inplace=True)

    # # Save table of clone names, accessions, and ids.
    # clone_name_ids_and_accessions.to_csv(args.annotated_accessions, sep="\t", index=False)
