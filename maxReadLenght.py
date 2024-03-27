from Bio import SeqIO

# Caminho para o arquivo FASTQ de entrada e saída
input_file = "entrada.fastq"
output_file = "saida.fastq"

# Lista para armazenar as reads filtradas
reads_filtradas = []

# Iterar sobre as reads do arquivo FASTQ de entrada
for record in SeqIO.parse(input_file, "fastq"):
    # Verificar se o comprimento da sequência é menor ou igual a 40
    if len(record.seq) <= 40:
        # Se sim, adicionar a read à lista de reads filtradas
        reads_filtradas.append(record)

# Escrever as reads filtradas em um novo arquivo FASTQ
with open(output_file, "w") as output_handle:
    SeqIO.write(reads_filtradas, output_handle, "fastq")

print("Filtragem concluída. As reads com 40 ou menos nucleotídeos foram salvas em", output_file)
