import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="FASTA sequences to reverse complement")
    parser.add_argument("output", help="reverse complemented FASTA sequences")
    args = parser.parse_args()

    with open(args.output, "w") as oh:
        for seq_record in SeqIO.parse(args.input, "fasta"):
            new_record = SeqRecord(seq_record.seq.reverse_complement(), id="%s_rc" % seq_record.id, description="")
            SeqIO.write([new_record], oh, "fasta")
