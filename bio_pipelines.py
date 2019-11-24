import requests
import bs4
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from app import app


class BLASTSearch():
    """
    Uses NCBI BLAST to perform a search on the input sequence.
    self.query: a nucleotide or protein sequence
    """

    def __init__(self, query):
        """
        Given a query in the form of a string, perform a BLAST search and store
        the result in self.hits. By default, only the first 10 hits are
        requested.
        """

        self.query = Seq(query)

        # Search using NCBI Blast. hitlist_size determines the number of hits to return
        result_handle = NCBIWWW.qblast(
            "blastn", "nt", self.query, hitlist_size=10)

        blast_record = NCBIXML.read(result_handle)
        self.hits = blast_record.alignments
        self.export_hits()

    def export_hits(self):
        for hit in self.hits:
            print(hit)


# Debugging
# my_blast = BLASTSearch('atggtgcacctgactcctgaggagaagtctgccgttactgccctgtggggcaaggtgaacgtggatgaagttggtggtgaggccctgggcaggttgctggtggtctacccttggacccagaggttctttgagtcctttggggatctgtccactcctgatgctgttatgggcaaccctaaggtgaaggctcatggcaagaaagtgctcggtgcctttagtgatggcctggctcacctggacaacctcaagggcacctttgccacactgagtgagctgcactgtgacaagctgcacgtggatcctgagaacttcaggctcctgggcaacgtgctggtctgtgtgctggcccatcactttggcaaagaattcaccccaccagtgcaggctgcctatcagaaagtggtggctggtgtggctaatgccctggcccacaagtatcactaa')
# my_blast.query
# my_blast.hits
# for hit in my_blast.hits:
#     print(hit.title)
#     for hsp in hit.hsps:
#         print("E-value:", hsp.expect)


class WikiSearch():
    """
    An object for performing a search of Wikipedia and containing the results.
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


def tests():
    """
    For debugging purposes only
    """
    my_search = WikiSearch("Python")
    for each in my_search.get_hrefs():
        print(each)


def main():
    pass
