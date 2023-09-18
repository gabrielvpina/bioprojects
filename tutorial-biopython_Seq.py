# Importar todo o pacote biopython
import Bio

# Importar somente parte do pacote - ex. de manipulação de seq
from Bio.Seq import Seq

seq1 = Seq("ACTGTACGTTAC")

print(seq1)

# tipo da sequencia
print(type(seq1))
# classe 'Bio.Seq.Seq' - própria do biopython

# imprimir tamanho da sequencia
print(len(seq1))

# obter reverse complementar da sequencia
complementar = seq1.reverse_complement() 

# transcrever a sequencia 
rna = seq1.transcribe()
print("O rna é: ",rna)

# traduzir a sequencia
ptn = seq1.translate()
print("A proteína é: ",ptn)

#--------------------------------------------------

# Primeiro caractere de sequencia
print(seq1[1])

# Posiçao das bases da sequencia - i para posição e n para nucleotídeo
for i, n in enumerate(seq1):
    print(i, n)














