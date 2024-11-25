import json

# Dados simulados
data = [
    {"name": "ORF1", "start": 100, "end": 500},
    {"name": "ORF2", "start": 600, "end": 900},
    {"name": "ORF3", "start": 1000, "end": 1500},
]

# Salvar como arquivo JSON
with open("orfs.json", "w") as json_file:
    json.dump(data, json_file)

