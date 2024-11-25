from flask import Flask, render_template

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")  # Certifique-se de que o nome do template está correto


# Servir arquivos estáticos (JSON e CSS)
@app.route("/<path:path>")
def static_files(path):
    return app.send_static_file(path)

if __name__ == "__main__":
    app.run(debug=True)

