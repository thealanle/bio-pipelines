from flask_wtf import FlaskForm
from wtforms import StringField, RadioField, SubmitField
from wtforms.validators import DataRequired
from wtforms.widgets import TextArea


class QueryForm(FlaskForm):
    query = StringField(label='Query:', validators=[
                        DataRequired()], widget=TextArea())

    have = RadioField(label='Have:', choices=[
                      ('nt', 'Nucleotide Sequence'), ('pro', 'Protein Sequence')])

    want = RadioField(label='Want:', choices=[
                      ('nt', 'Nucleotide Sequence'), ('pro', 'Protein Sequence'), ('hits', 'BLAST Hits')])

    submit = SubmitField(label='Submit')
