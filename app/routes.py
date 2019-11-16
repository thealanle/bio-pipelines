from app import app
from flask import render_template, flash, redirect
from app.forms import QueryForm


@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='bio-pipelines')


@app.route('/query', methods=['GET', 'POST'])
def query():
    form = QueryForm()
    if form.validate_on_submit():
        flash(f"Submitted the following string: {form.query.data}")
        return redirect('/index')
    return render_template('query.html', title='Submit a Query', form=form)
