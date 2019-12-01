# bio-pipelines
Easily find information about nucleic acid and protein sequences.




# Features
Uses NCBI BLAST to gather and display relevant information about a sequence.

# Requirements
* Python 3.7.5
* [pipenv](https://github.com/pypa/pipenv)

# Setting Up the Web Application

1. Start by cloning into the repository:
```
git clone https://github.com/thealanle/bio-pipelines.git
```
Or by downloading [the latest release.](https://github.com/thealanle/bio-pipelines/releases)

2. Create a new environment in the application folder:
```
pipenv --python 3.7.5
```

3. Install the required packages using the repository's [Pipfile.lock](../Pipfile.lock) using `pipenv update`.
Then enter the environment using `pipenv shell`
4. **bio-pipelines** is built using Flask, so set FLASK_APP (required) and FLASK_DEBUG (optional).
* Windows
```
PS C:\...\bio-pipelines> Set-Variable -Name FLASK_APP -Value "bio_pipelines.py"
PS C:\...\bio-pipelines> Set-Variable -Name FLASK_DEBUG -Value 1
```

* Linux
```
(bio-pipelines): export FLASK_APP=bio-pipelines.py
(bio-pipelines): export FLASK_DEBUG=1
```
5. Start the server using `flask run`.

**bio-pipelines** should be up and running! Connect to `localhost:5000` using your browser of choice.
