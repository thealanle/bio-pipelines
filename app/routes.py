from app import app
from flask import render_template, flash, redirect, url_for
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
        return redirect(url_for('result'))
    return render_template('query.html', title='Submit a Query', form=form)


@app.route('/result')
def result():
    return render_template('result.html', title='Result')
