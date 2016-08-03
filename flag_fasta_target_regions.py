import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("slop", type=int)
    parser.add_argument("flagged")
    parser.add_argument("--primer3", action="store_true", help="output in primer3-compatible input format")
    parser.add_argument("--min_primer_size", type=int, default=18)
    parser.add_argument("--opt_primer_size", type=int, default=20)
    parser.add_argument("--max_primer_size", type=int, default=27)
    parser.add_argument("--max_primers_to_return", type=int, default=5)
    args = parser.parse_args()

    slop = args.slop
    alphabet = Alphabet()
    alphabet.letters = ["A", "T", "C", "G", "[", "]"]

    with open(args.flagged, "w") as fh:
        for seq_record in SeqIO.parse(args.fasta, "fasta"):
            if args.primer3:
                sequence = str(seq_record.seq)
                fh.write("SEQUENCE_ID=%s\n" % seq_record.id)
                fh.write("SEQUENCE_TEMPLATE=%s\n" % sequence)
                fh.write("SEQUENCE_TARGET=%i,%i\n" % (slop, len(sequence) - 2 * slop))
                fh.write("PRIMER_OPT_SIZE=%i\n" % args.opt_primer_size)
                fh.write("PRIMER_MIN_SIZE=%i\n" % args.min_primer_size)
                fh.write("PRIMER_MAX_SIZE=%i\n" % args.max_primer_size)
                fh.write("PRIMER_NUM_RETURN=%i\n" % args.max_primers_to_return)
                fh.write("=\n")
            else:
                sequence = Seq(str(seq_record.seq)[:slop] + "[" + str(seq_record.seq)[slop:-slop] + "]" + str(seq_record.seq)[-slop:])
                sequence.alphabet = alphabet
                flagged_seq_record = SeqRecord(sequence, id=seq_record.id, description="")
                SeqIO.write([flagged_seq_record], fh, "fasta")
