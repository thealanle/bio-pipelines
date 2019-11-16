from app import app
from flask import render_template
from app.forms import QueryForm


@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='bio-pipelines')


@app.route('/query')
def query():
    form = QueryForm()
    return render_template('query.html', title='Submit a Query', form=form)
