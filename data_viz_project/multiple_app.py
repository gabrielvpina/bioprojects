from flask import Flask, render_template
import os

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

# Servir arquivos estáticos como JSON
@app.route("/static/<path:path>")
def send_static(path):
    return app.send_static_file(path)

if __name__ == "__main__":
    app.run(debug=True)

