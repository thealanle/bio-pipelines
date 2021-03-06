from flask import render_template, flash, redirect, url_for, request
from app import app
from app.forms import QueryForm
from Bio.Seq import Seq
from Bio import SeqIO
import bio_pipelines


@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='bio-pipelines')


@app.route('/query', methods=['GET', 'POST'])
def query():
    form = QueryForm()
    return render_template('query.html', title='Submit a Query', form=form)


@app.route('/result', methods=['POST'])
def result():
    """
    Based on form data, perform the appropriate search or conversion operation
    and return results.
    """
    if request.method == 'POST':
        print("\n>>>>>POST request received. Processing...")

        # Move form data into a dict
        form = {'query': str(request.form.get('query')),
                'have': request.form.get('have'),
                'want': request.form.get('want')}

        RESULT_TYPES = ['nt', 'pro', 'hits', 'wiki']

        # Create a dict of results, initialized to None
        results = {each: None for each in RESULT_TYPES}

        # If a radio button is unselected, return to a fresh Query page.
        if not (form['have'] and form['want']):
            return render_template('query.html', form=QueryForm(), error_message="You must choose a Have and a Want.")

        """
        # This functionality will be reworked into the Results table.
        # Based on the form entry, update the "results" accordingly
        if form['want'] == 'wiki':  # Wikipedia search results
            results['wiki'] = bio_pipelines.WikiSearch(
                form['query']).get_hrefs()
        """

        # Have nt, want nt (do nothing)
        if form['have'] == 'nt' and form['want'] == 'nt':
            results['nt'] = Seq(form['query'])

        # Have nt, want pro (translation)
        elif form['have'] == 'nt' and form['want'] == 'pro':  # Protein translation
            print(f"Translating nucleic acid string to protein...")
            form['query'] = form['query'].replace('\r', '')

            if form['query'][0] == '>':  # Check if FASTA-formatted
                print('>>>>>FASTA format detected')
                seq = ''.join([line.strip()
                               for line in form['query'].split('\n')][1:])
                form['query'] = Seq(seq)
            else:
                print('>>>>>Sequence not in FASTA format. Sanitizing...')
                seq = ''.join([line.strip()
                               for line in form['query'].split('\n')])
                form['query'] = Seq(seq)

            print(f">>>>>Now translating sequence:\n{form['query']}")
            results['pro'] = form['query'].translate()

        # Have nt, want BLAST (blastn)
        elif form['have'] == 'nt' and form['want'] == 'hits':  # BLAST hits
            print(f"Getting BLAST data using input: {form['query']}")
            results['hits'] = bio_pipelines.BLASTSearch(
                query=form['query'], data_type=form['have']).data_table

        # Have pro, want nt (tblastn)
        elif form['have'] == 'pro' and form['want'] == 'nt':  # BLAST hits
            print(f"Getting BLAST data using input: {form['query']}")
            results['hits'] = bio_pipelines.BLASTSearch(
                query=form['query'], data_type='pro-nt').data_table

        # Have pro, want BLAST (blastp)
        elif form['have'] == 'pro' and form['want'] == 'hits':  # BLAST hits
            print(f"Getting BLAST data using input: {form['query']}")
            results['hits'] = bio_pipelines.BLASTSearch(
                query=form['query'], data_type=form['have']).data_table

        print(f">>>>>POST request received. Rendering result.html...")
        # Return results to render_template for HTML display
        return render_template('result.html', title='Results', query=form['query'], have=form['have'], want=form['want'], nucleotide_seq=results['nt'], protein_seq=results['pro'], blast_results=results['hits'], wiki_results=results['wiki'])

    else:
        return redirect('index')


@app.route('/about')
def about():
    """
    Return the About page.
    """
    return render_template('about.html', title='About')


@app.route('/test')
def test():
    """
    Return the Test page. This should only be accessed directly via the URL.
    This page may be used for testing CSS or incomplete functions.
    """
    return render_template('test.html', title='Test')


# def sanitize(query):
#     if query[0] == '>':
#         return query.split('\n')[1:]
#     else:
#         return ''.join(query.split('\n'))
