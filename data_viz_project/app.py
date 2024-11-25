from flask import Flask, render_template, jsonify
import os

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/get_orfs_data")
def get_orfs_data():
    # Optional: Add a route to debug if JSON is loading
    try:
        with open('static/orfs.json', 'r') as file:
            data = json.load(file)
        return jsonify(data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)