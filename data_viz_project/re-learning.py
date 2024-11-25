import os
import glob
import csv
import json
from Bio import SeqIO

def generate_orf_json(input_dir, output_file, suffixes, num_largest_orfs):

    def find_orf_files(suffixes):
        orf_fasta_files = []
        for suffix in suffixes:
            files = glob.glob(os.path.join(input_dir, f"*{suffix}"))
            if files:
                orf_fasta_files.append((files[0], f"Genetic Code {suffix.split('gc')[1].split('.')[0]}"))
        return orf_fasta_files
    

    def parse_original_fasta(input_dir):
        """
        Processa os arquivos FASTA originais e retorna os comprimentos das contigs.
        """
        original_lengths = {}
        fasta_files = glob.glob(os.path.join(input_dir, "*_merged.fasta"))

        for fastafile in fasta_files:
            for record in SeqIO.parse(fastafile, "fasta"):
                original_header = record.description
                contig_name = original_header.split()[0]  # Obtém o nome do contig
                original_length = len(record.seq)
                original_lengths[contig_name] = original_length  # Armazena o tamanho da sequência

        return original_lengths


    # Parse original lenghts
    original_lengths = parse_original_fasta(input_dir)

    def parse_orf_fastas(file_paths):
        orf_data_by_contig = {}
        for file_path, code in file_paths:
            for record in SeqIO.parse(file_path, "fasta"):
                header = record.description
                parts = header.split()
                
                contig_full = parts[0]
                if '_ORF.' in contig_full:
                    contig = contig_full[:contig_full.index('_ORF.')]
                else:
                    contig = contig_full
                
                coordinates_strand = parts[1].replace('[', '').replace(']', '')
                if '(' in coordinates_strand:
                    coordinates, strand = coordinates_strand.split('(')
                    coordinates = coordinates.strip()
                    strand = strand.replace(')', '').strip()
                    start, end = map(int, coordinates.split('-'))
                else:
                    continue

                additional_info = parts[2:]
                additional_info_dict = {info.split(':')[0]: info.split(':')[1] for info in additional_info}

                orf = {
                    'name': contig_full,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'sequence': str(record.seq),
                    'type': additional_info_dict.get('type'),
                    'length': int(additional_info_dict.get('length', 0)),
                    'frame': additional_info_dict.get('frame'),
                    'start_codon': additional_info_dict.get('start'),
                    'stop_codon': additional_info_dict.get('stop'),
                    'code': code,
                    'domains': []  # Será preenchido posteriormente
                }
                
                if contig not in orf_data_by_contig:
                    orf_data_by_contig[contig] = {
                        'contig': contig,
                        'length': original_lengths.get(contig, 0),
                        'orfs': []
                    }
                orf_data_by_contig[contig]['orfs'].append(orf)

        return orf_data_by_contig



    def parse_cdd_files(input_dir):
            cdd_data = {}
            cdd_suffixes = ['_ORFgc1_CDD.tsv', '_ORFgc5_CDD.tsv', '_ORFgc11_CDD.tsv']
            for suffix in cdd_suffixes:
                files = glob.glob(os.path.join(input_dir, f"*{suffix}"))
                for file in files:
                    with open(file, 'r') as f:
                        reader = csv.reader(f, delimiter='\t')
                        for row in reader:
                            orf_name = row[0]
                            database = row[3]
                            domain_code = row[4]
                            domain_name = row[5]
                            e_value = row[8]

                            # remover células em branco
                            if domain_name != '-' and e_value != '-' and database in ['CDD', 'Pfam', 'SUPERFAMILY']:
                                if orf_name not in cdd_data:
                                    cdd_data[orf_name] = []
                                cdd_data[orf_name].append({
                                    'database': database,
                                    'domain_code': domain_code,
                                    'domain_name': domain_name,
                                    'e_value': e_value
                                })
            return cdd_data

    # Encontrar e parsear arquivos ORF
    orf_files = find_orf_files(suffixes)
    orf_data_by_contig = parse_orf_fastas(orf_files)

    # Parsear dados CDD
    cdd_data = parse_cdd_files(input_dir)

    # Associar dados CDD aos ORFs
    for contig, contig_data in orf_data_by_contig.items():
        for orf in contig_data['orfs']:
            orf['domains'] = cdd_data.get(orf['name'], [])

    # Identificar e preparar JSON
    orf_json = []
    for contig, contig_data in orf_data_by_contig.items():
        # Ordenar ORFs por comprimento e selecionar as maiores
        orfs = sorted(contig_data['orfs'], key=lambda x: x['length'], reverse=True)[:num_largest_orfs]
        contig_length = contig_data['length']
        
        orf_json.append({
            'contig': contig,
            'length': contig_length,
            'orfs': orfs
        })

    # Salvar JSON
    with open(os.path.join(input_dir, output_file), 'w') as f:
        json.dump(orf_json, f, indent=4)

# Exemplo de uso
input_directory = "test_file/"
output_json = "correct.json"
suffixes = ["gc1.fasta", "gc5.fasta", "gc11.fasta"]
num_largest_orfs = 3  # Definir o número de maiores ORFs por contig
generate_orf_json(input_directory, output_json, suffixes, num_largest_orfs)