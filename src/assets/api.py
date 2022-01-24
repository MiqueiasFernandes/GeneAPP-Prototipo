from flask import Flask
from flask import request
from flask import send_from_directory
from flask_cors import CORS
import os

app = Flask(__name__)
CORS(app)

@app.route("/efetch")
def ncbi_mirror():
    db = request.args.get('db')
    rettype = request.args.get('rettype')
    id = request.args.get('id')
    strand = request.args.get('strand')
    seq_start = request.args.get('seq_start')
    seq_stop = request.args.get('seq_stop')
    key = request.args.get('api_key')
    retmode = 'text' # request.args.get('retmode')
    suffix = "" if db != "nuccore" else  f".{seq_start}.{strand}.{seq_stop}"
    file = f"static/{db}/{rettype}/{id+suffix}"
    extra = "" if db != "nuccore" else f"&strand={strand}&seq_start={seq_start}&seq_stop={seq_stop}"
    
    cmd = f"wget -T 10 -O {file} 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db}&rettype={rettype}&retmode={retmode}{extra}&id={id}&api_key={key}'"
    if os.path.exists(file):
        if len(open(file).readlines()) < 3:
            os.remove(file)
        else:
            cmd = None

    if not cmd is None:
        os.system(cmd)

    return send_from_directory(f"static/{db}/{rettype}", id+suffix)

## https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_015449.3&strand=1&seq_start=2915348&seq_stop=2915569&rettype=fasta&retmode=text