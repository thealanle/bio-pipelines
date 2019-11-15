import requests
import bs4
from Bio.Blast import NCBIWWW, NCBIXML

INPUT_FILE = 'input_files/hbb.fasta'
BLAST_OUTPUT = 'output_files/my_blast.xml'
#
# with open('input_files/hbb.fasta', 'r') as f:
#     title = ''
#     sequence = ''
#     for line in f.readlines():
#         if line[0] == '>':
#             title = line.strip()
#         else:
#             sequence += line.strip()

with open(INPUT_FILE, 'r') as f:
    fasta_string = f.read()
print(fasta_string)

result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

result_xml = result_handle.read()

print(result_xml)
with open(BLAST_OUTPUT, "w") as f_out:
    f_out.write(result_xml)
