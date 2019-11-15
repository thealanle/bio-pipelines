import requests
import bs4
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# BLAST Functionality

INPUT_FILE = 'input_files/hbb.fasta'
BLAST_OUTPUT = 'output_files/my_blast.xml'

fasta_record = SeqIO.read(INPUT_FILE, format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_record.seq)

# NCBIXML.parse returns a generator. Here it is converted to a list for tests.
blast_records = list(NCBIXML.parse(result_handle))
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(hsp.expect)
print(blast_records[0].query)
# with open(BLAST_OUTPUT, "w") as f_out:
#     f_out.write(result_xml)


# Transcription and Translation
def get_transcript(sequence):
    return sequence.transcribe()


def get_translation(sequence):
    return sequence.translate()


print(get_transcript(fasta_record.seq))
print(get_translation(fasta_record.seq))


# EXPASY_TRANSLATE_URL = "https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
# query = {'dna_sequence': fasta_record, 'output_format': 'fasta'}
#
# expasy_response = requests.post(EXPASY_TRANSLATE_URL, params=query)
# print(expasy_response)
# print(expasy_response.text)
