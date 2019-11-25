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
    if form.validate_on_submit():
        data = {'query': form.query.data,
                'have': form.have.data, 'want': form.want.data}
        flash(f"Submitted the following string: {data['query']}")
        print(data)
        # return render_template('result.html', form=form, sequence=form.query.data, have=form.have.data, want=form.want.data)
    return render_template('query.html', title='Submit a Query', form=form)


@app.route('/result', methods=['POST'])
def result():
    if request.method == 'POST':
        print("\n>>>>>POST request received. Processing...")

        form = {'query': str(request.form.get('query')),
                'have': request.form.get('have'),
                'want': request.form.get('want')}

        RESULT_TYPES = ['gene', 'protein', 'hits', 'wiki']

        results = {each: None for each in RESULT_TYPES}

        if form['want'] == 'wiki':
            results['wiki'] = bio_pipelines.WikiSearch(
                form['query']).get_hrefs()

        elif form['want'] == 'hits':
            print(f"Getting BLAST data using input: {form['query']}")
            results['hits'] = bio_pipelines.BLASTSearch(form['query']).hits

        elif form['want'] == 'protein':
            print(f"Converting nucleic acid string to protein...")
            form['query'] = form['query'].replace('\r', '')

            if form['query'][0] == '>':
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
            results['protein'] = form['query'].translate()

        print(f">>>>>POST request received. Rendering result.html...")
        return render_template('result.html', title='Results', query=form['query'], have=form['have'], want=form['want'], protein_seq=results['protein'], blast_results=results['hits'], wiki_results=results['wiki'])

    else:
        return redirect('index')


@app.route('/about')
def about():
    return render_template('about.html', title='About')


@app.route('/test')
def test():
    return render_template('test.html', title='Test')


def sanitize(query):
    if query[0] == '>':
        return query.split('\n')[1:]
    else:
        return ''.join(query.split('\n'))
