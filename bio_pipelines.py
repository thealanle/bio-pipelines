import requests
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

    def __init__(self, query, data_type='nt'):
        """
        Given a query in the form of a string, perform a BLAST search and store
        the result in self.hits. By default, only the first 10 hits are
        requested.
        """

        # Dictionary of  services and their NSID abbreviations
        self.NSID_TABLE = {
            'bbm': 'GenInfo Backbone',
            'bbs': 'GenInfo Backbone',
            'dbj': 'DNA Database of Japan',
            'emb': 'EMBL',
            'gb': 'GenBank',
            'gi': 'GenInfo',
            'gim': 'GenInfo Import',
            'gnl': 'General',
            'gp': 'GenPept',
            'lcl': 'Local',
            'oth': 'Other',
            'pat': 'Patent',
            'pdb': 'Brookhaven Protein Database',
            'pir': 'Protein Information Resource International',
            'prf': 'Protein Research Foundation',
            'ref': 'RefSeq',
            'sp': 'SWISS-PROT',
            'tpd': 'Third-party annotation, DDBJ',
            'tpe': 'Third-party annotation, EMBL',
            'tpg': 'Third party annotation, GenBank'
        }

        # Specify which program to use for nucleotide vs. protein BLAST
        self.PROGRAM_TABLE = {
            'nt': 'blastn',
            'pro': 'blastp',
        }

        # Specify which database to use for nucleotide vs. protein BLAST
        self.DB_TABLE = {
            'nt': 'nt',
            'pro': 'swissprot',
        }

        self.titles = []

        self.query = Seq(query)

        # Search using NCBI Blast. hitlist_size determines the number of hits to return
        result_handle = NCBIWWW.qblast(
            program=self.PROGRAM_TABLE[data_type], database=self.DB_TABLE[data_type], sequence=self.query, hitlist_size=10)

        self.blast_record = NCBIXML.read(result_handle)
        # print(type(blast_record))  # Returns <class 'Bio.Blast.Record.Blast'>

        # Create a list of dicts containing the alignments and their properties
        self.hits = self.blast_record.alignments

        self.parse_headers(self.hits)
        self.build_table()

    def parse_headers(self, alignments):
        """
        Given a list of NCBI FASTA-formatted alignment, return a list of
        dicts for each alignment.
        """
        self.parsed_headers = []
        for alignment in alignments:
            self.parsed_headers.append(self.parse_header(alignment.title))

        # Print the resulting data to stdout
        print("Title strings parsed and stored in parsed_headers. Results: \n")
        for entry in self.parsed_headers:
            for key, value in entry.items():
                if value:
                    print(f"{key}: {value}")
            print()

    def parse_header(self, header):
        """
        Parse a single header string.
        """

        # Build a dict based off of all NSIDs, default value of None
        result = {id: None for id in self.NSID_TABLE.keys()}

        # Split header by '|'
        entries = [each.strip() for each in header.split('|')]

        # If an entry is an NSID, build an href to the appropriate search page.
        for i in range(len(entries)):
            entry = entries[i]
            # Check if entry is an NSID
            if entry in result.keys():
                try:
                    # TO-DO: This may vary depending on formatting of NSID.
                    # Consider making a table to indicate correct indices.
                    result[entry] = self.build_href(entry, entries[i + 1])
                    print(f"Parsed {header} into {result[entry]}")
                except Exception:
                    result[entry] = None

        try:
            # Most sequence titles are at position 4
            result['title'] = entries[4]
        except Exception:
            # Take the longest string as the title
            result['title'] = sorted(entries, key=len)[-1]

        self.titles.append(result['title'])
        return result

    def build_href(self, nsid, id):
        """
        Goal: Given an NSID type and ID number, return an href linking to an
        NCBI search for the id.
        """

        url = f"https://www.ncbi.nlm.nih.gov/search/all/?term={id}"
        return f"<a href=\"{url}\">{self.NSID_TABLE[nsid]}: {id}</a>"

    def build_table(self):
        self.data_table = [['Title', 'References']]

        # Iterate through each parsed header and add its contents to the table
        for d in self.parsed_headers:
            # Begin building a table row
            row = [d['title']]

            # Join references with a break and add to a single cell of the table
            refs = '<br>'.join(
                [value for key, value in d.items() if value and not key == 'title'])

            row.append(refs)
            self.data_table.append(row)


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
        Return a list of HTML - formatted anchor tags.
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
