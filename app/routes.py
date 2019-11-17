from app import app
from flask import render_template, flash, redirect, url_for, request
from app.forms import QueryForm


@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='bio-pipelines')


@app.route('/query', methods=['GET', 'POST'])
def query():
    form = QueryForm()
    if form.validate_on_submit():
        data = {'search_term': form.query.data,
                'have': form.have.data, 'want': form.want.data}
        flash(f"Submitted the following string: {data['search_term']}")
        print(data)
        # return redirect(url_for('result'))
        # return render_template('result.html', form=form, sequence=form.query.data, have=form.have.data, want=form.want.data)
    return render_template('query.html', title='Submit a Query', form=form)


@app.route('/result', methods=['POST'])
def result():
    form = request.form
    if request.method == 'POST':
        print(">>>>>POST request received. Rendering result.html...")
        return render_template('result.html', title='Results', form=form, sequence=form.get('query'), have=form.get('have'), want=form.get('want'))


@app.route('/about')
def about():
    return render_template('about.html', title='About')
