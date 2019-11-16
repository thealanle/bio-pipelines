from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired


class QueryForm(FlaskForm):
    query = StringField(label='Query', validators=[DataRequired()])
    submit = SubmitField(label='Submit')
