import itertools

# Genetic code dictionary
genetic_code = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def dna_to_mrna(dna):
    return dna.replace('T', 'U')

def complement_dna(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in dna)

def translate_dna_to_protein(dna):
    if len(dna) % 3 != 0:
        raise ValueError("DNA sequence length must be a multiple of 3")
    
    mrna = dna_to_mrna(dna)
    protein = ''.join(genetic_code[mrna[i:i+3]] for i in range(0, len(mrna), 3))
    return protein

def get_codon_frequencies(amino_acids):
    codon_freq = {}
    for aa in amino_acids:
        codons = [codon for codon, amino in genetic_code.items() if amino == aa]
        codon_freq[aa] = {codon: 0 for codon in codons}
    return codon_freq

def analyze_codon_frequency(dna, amino_acids):
    if len(amino_acids) > 3:
        raise ValueError("Input should be max. 3 Aminoacids")
    
    mrna = dna_to_mrna(dna)
    codon_freq = get_codon_frequencies(amino_acids)
    
    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]
        if codon in itertools.chain(*codon_freq.values()):
            for aa, codons in codon_freq.items():
                if codon in codons:
                    codon_freq[aa][codon] += 1
    
    return codon_freq

def get_all_codons(amino_acids):
    codon_dict = {}
    for aa in amino_acids:
        codons = [codon for codon, amino in genetic_code.items() if amino == aa]
        codon_dict[aa] = codons
    return codon_dict

def task1():
    print("1. Translate a DNA sequence -> mRNA (using 'U' instead of 'T') -> into an aminoacid sequence (protein)")
    dna = input("Input DNA = ").upper()
    
    if len(dna) % 3 != 0:
        print("Error: DNA sequence length must be a multiple of 3")
        return
    
    print(f"Input DNA = {dna}")
    print(f"mRNA = {dna_to_mrna(dna)}")
    print(f"Aminoacid = {translate_dna_to_protein(dna)}")

def task2():
    print("\n2. Provides the frequency of each RNA codon encoding given aminoacids, in a DNA sequence")
    amino_acids = input("Input Aminoacid = ").upper()
    
    if len(amino_acids) > 3:
        print("Error: Input should be max. 3 Aminoacids")
        return
    
    dna = input("Input DNA sequence = ").upper()
    mrna = dna_to_mrna(dna)
    
    print(f"Input Aminoacid = {amino_acids}")
    print(f"mRNA = {mrna}")
    
    all_codons = get_all_codons(amino_acids)
    codon_freq = analyze_codon_frequency(dna, amino_acids)
    
    for aa, codons in all_codons.items():
        print(f"\nCodons for {aa}:")
        for codon in codons:
            count = codon_freq[aa][codon] if aa in codon_freq and codon in codon_freq[aa] else 0
            print(f"{codon} = {count}")

if __name__ == "__main__":
    task1()
    task2()