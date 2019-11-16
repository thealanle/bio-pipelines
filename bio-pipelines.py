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
# Search using NCBI Blast. hitlist_size determines the number of hits to return
result_handle = NCBIWWW.qblast(
    "blastn", "nt", fasta_record.seq, hitlist_size=10)

# NCBIXML.parse returns a generator. Here it is converted to a list for tests.
blast_records = list(NCBIXML.parse(result_handle))
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(hsp.expect)


# Transcription and Translation
def get_transcript(sequence):
    return sequence.transcribe()


def get_translation(sequence):
    return sequence.translate()


# Wikipedia Search
def wiki_search(query):
    WIKI_API_URL = "https://en.wikipedia.org/w/api.php"
    WIKI_PAGEID_URL = "https://en.wikipedia.org/?curid="

    params = {
        "action": "query",
        "format": "json",
        "list": "search",
        "srsearch": query
    }

    session = requests.Session()
    response = session.get(url=WIKI_API_URL, params=params)
    data = response.json()

    top_result = data['query']['search'][0]
    if top_result:
        print(f"Top result on Wikipedia: {top_result['title']} at {WIKI_PAGEID_URL + str(top_result['pageid'])}")
    else:
        print(f"{query} returned 0 results.")


def main():
    print(f"Original sequence:\n{fasta_record.seq}")
    print(f"Transcription:\n{get_transcript(fasta_record.seq)}")
    print(get_translation(fasta_record.seq))

    wiki_search("Python")

# EXPASY_TRANSLATE_URL = "https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
# query = {'dna_sequence': fasta_record, 'output_format': 'fasta'}
#
# expasy_response = requests.post(EXPASY_TRANSLATE_URL, params=query)
# print(expasy_response)
# print(expasy_response.text)
