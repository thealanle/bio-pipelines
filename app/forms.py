from flask_wtf import FlaskForm
from wtforms import StringField, RadioField, SubmitField
from wtforms.validators import DataRequired


class QueryForm(FlaskForm):
    query = StringField(label='Query:', validators=[DataRequired()])
    have = RadioField(label='Have:', choices=[
                      ('dna', 'DNA'), ('mrna', 'mRNA'), ('protein', 'Protein')])
    want = RadioField(label='Want:', choices=[
                      ('dna', 'DNA'), ('mrna', 'mRNA'), ('protein', 'Protein'),
                      ('wiki', 'Wikipedia Results'), ('diseases', 'Related Diseases'), ('hits', 'BLAST Hits')])
    submit = SubmitField(label='Submit')
