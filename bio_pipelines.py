import requests
import bs4
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from app import app


# Transcription and Translation
def get_transcript(sequence):
    """
    Given a Seq object, return an mRNA transcript. Note that the given sequence is assumed to be DNA.
    """
    return sequence.transcribe()


def get_translation(sequence):
    """
    Given a Seq object, return a protein translation. Note that the given sequence is assumed to be mRNA.
    """
    return sequence.translate()


# BLAST Functionality
def blast_search(query):
    """
    Uses NCBI BLAST to perform a search on the input sequence(s).
    query: a file handle for a FASTA record
    """
    # INPUT_FILE = 'input_files/hbb.fasta'
    # BLAST_OUTPUT = 'output_files/my_blast.xml'

    fasta_record = SeqIO.read(query, format="fasta")
    # Search using NCBI Blast. hitlist_size determines the number of hits to return
    result_handle = NCBIWWW.qblast(
        "blastn", "nt", fasta_record.seq, hitlist_size=10)

    # NCBIXML.parse returns a generator. Here it is converted to a list for tests.
    blast_records = list(NCBIXML.parse(result_handle))
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print(hsp.expect)


# Wikipedia Search
class WikiSearch():
    """
    An object for performing a search and containing the results.
    """

    def __init__(self, search_term):
        """
        On instantiation, perform a search on Wikipedia via its API and store
        the returned articles in self.articles.
        """
        WIKI_API_URL = "https://en.wikipedia.org/w/api.php"
        WIKI_PAGEID_URL = "https://en.wikipedia.org/?curid="

        self.params = {
            "action": "query",
            "format": "json",
            "list": "search",
            "srsearch": search_term
        }

        session = requests.Session()
        response = session.get(url=WIKI_API_URL, params=self.params).json()
        self.results = response['query']['search']
        self.articles = []

        if self.results:
            for result in response['query']['search']:
                self.articles.append(
                    (result['title'], WIKI_PAGEID_URL + str(result['pageid'])))
        else:
            print(f"{search_term} returned 0 results.")

    def get_hrefs(self):
        """
        Return a list of HTML-formatted anchor tags.
        """
        return [f"""<a href="{article[1]}">{article[0]}</a>""" for article in self.articles]


def main():
    print(f"Original sequence:\n{fasta_record.seq}")
    print(f"Transcription:\n{get_transcript(fasta_record.seq)}")
    print(get_translation(fasta_record.seq))

    my_search = WikiSearch("Python")
    for each in my_search.get_hrefs():
        print(each)
