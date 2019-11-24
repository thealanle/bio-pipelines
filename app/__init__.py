from flask import Flask
from config import Config
from flask_bootstrap import Bootstrap
import os


app = Flask(__name__, static_url_path='/static', static_folder=os.getcwd())
app.config.from_object(Config)
app.config['BOOTSTRAP_SERVE_LOCAL'] = True  # This turns file serving static
bootstrap = Bootstrap(app)
from app import routes  # noqa:F401
