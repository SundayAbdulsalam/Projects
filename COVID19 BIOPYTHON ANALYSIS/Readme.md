## 🧬 COVID-19 Genome and Protein Structure Analysis Using Python

## Overview

This project performs a bioinformatics analysis of a viral genome using Python and Biopython. The workflow retrieves a nucleotide sequence from the NCBI database, performs genome-level analysis, transcribes and translates the sequence, analyzes amino acid composition, identifies long protein sequences, performs a BLAST search against protein structures, and visualizes a selected protein structure.

The analysis is implemented in the Jupyter Notebook `covid19.ipynb`.

## Project Workflow

The project follows these major steps:

1. Retrieve the genome sequence from NCBI.
2. Calculate basic genome statistics.
3. Calculate molecular weight.
4. Calculate GC content.
5. Analyze nucleotide distribution.
6. Visualize nucleotide composition.
7. Transcribe the DNA sequence into mRNA.
8. Translate the mRNA sequence into an amino acid sequence.
9. Analyze and visualize amino acid distribution.
10. Remove stop codons and split the translated sequence into protein sequences.
11. Filter protein sequences based on length.
12. Identify the longest protein sequence.
13. Save the selected protein sequence in FASTA format.
14. Perform a BLAST search against the Protein Data Bank.
15. Download the selected protein structure.
16. Parse the PDB structure using Biopython.
17. Identify protein chains and residues.
18. Visualize the protein structure interactively using NGLView.

## Dataset

The nucleotide sequence is retrieved from the NCBI Nucleotide database using the accession ID:

```text
MT121215.1
```

The sequence is retrieved in GenBank format using Biopython's `Entrez` module.

## Installation

Install the required Python libraries before running the notebook:

```bash
pip install biopython matplotlib nglview
```

Alternatively, inside a Jupyter Notebook:

```python
%pip install biopython matplotlib nglview
```

## Required Libraries

The project uses the following libraries:

* **Biopython** — sequence retrieval, sequence analysis, BLAST, FASTA, and PDB parsing.
* **Matplotlib** — visualization of nucleotide and amino acid distributions.
* **NGLView** — interactive three-dimensional visualization of protein structures.
* **Collections Counter** — counting amino acid frequencies.
* **urllib.request** — downloading the PDB structure file.

## Genome Retrieval

The genome sequence is downloaded from the NCBI Nucleotide database using `Entrez.efetch()`.

```python
from Bio import Entrez, SeqIO

Entrez.email = "your_email@example.com"

handle = Entrez.efetch(
    db="nucleotide",
    id="MT121215.1",
    rettype="gb",
    retmode="text"
)

record = SeqIO.read(handle, "genbank")
handle.close()
```

The notebook then displays:

* Sequence ID
* Sequence description
* Genome length
* Number of genomic features

## Genome Analysis

### Molecular Weight

The molecular weight of the genome is calculated using:

```python
from Bio.SeqUtils import molecular_weight

mw = molecular_weight(cold_dna, seq_type="DNA")
```

### GC Content

GC content is calculated using:

```python
from Bio.SeqUtils import gc_fraction

gc_content = gc_fraction(cold_dna)
```

### Nucleotide Distribution

The frequency of the four DNA nucleotides is calculated:

* Adenine (A)
* Thymine (T)
* Guanine (G)
* Cytosine (C)

The results are visualized using a bar chart.

## Transcription and Translation

The DNA genome is transcribed into mRNA:

```python
cold_mrna = cold_dna.transcribe()
```

The mRNA sequence is then translated into an amino acid sequence:

```python
protein_sequence = cold_mrna.translate()
```

## Amino Acid Analysis

The project counts the occurrence of each amino acid using Python's `Counter`:

```python
from collections import Counter

common_amino_acid = Counter(protein_sequence)
```

Stop codons are removed from the amino acid count before visualization:

```python
del common_amino_acid["*"]
```

The amino acid distribution is then displayed as a bar chart.

## Protein Sequence Extraction

The translated sequence is split at stop codons:

```python
protein_sequence_no_stop = protein_sequence.split("*")
```

Protein sequences shorter than 20 amino acids are removed:

```python
protein_sequence_no_stop = [
    protein
    for protein in protein_sequence_no_stop
    if len(protein) >= 20
]
```

The remaining sequences represent protein candidates that are at least 20 amino acids long.

## Identification of the Longest Protein

The protein sequences are sorted by length:

```python
top_5_proteins = sorted(
    protein_sequence_no_stop,
    key=len
)
```

The longest protein sequence is selected for further analysis:

```python
longest_protein = top_5_proteins[-1]
```

## FASTA File Generation

The selected protein sequence is saved in FASTA format:

```python
with open("protein_seq.fasta", "w") as file:
    file.write(f">cold protein\n{longest_protein}\n")
```

## BLAST Search

The selected protein sequence is submitted to NCBI BLAST using the `blastp` program against the Protein Data Bank database:

```python
from Bio.Blast import NCBIWWW

result_handle = NCBIWWW.qblast(
    "blastp",
    "pdb",
    protein_seq.seq
)
```

The BLAST results are read and analyzed using Biopython's `SearchIO` module.

The notebook displays information such as:

* Hit ID
* Hit description
* E-value
* Bit score
* Hit length
* Sequence alignment

## Protein Structure Analysis

A protein structure with the PDB ID `9SAP` is selected for structural visualization.

### Downloading the PDB File

The structure is downloaded directly from the RCSB Protein Data Bank:

```python
from urllib.request import urlretrieve

url = "https://files.rcsb.org/download/9SAP.pdb"
filename = "9SAP.pdb"

urlretrieve(url, filename)
```

## Reading the PDB File

The PDB structure is parsed using Biopython:

```python
from Bio.PDB import PDBParser

parser = PDBParser()

structure = parser.get_structure(
    "9SAP",
    "9SAP.pdb"
)
```

## Protein Chain Analysis

The notebook identifies the chains and counts the residues in each chain:

```python
for chain in structure[0]:
    print(
        f"Chain ID: {chain.id}, "
        f"Number of residues: {len(chain)}"
    )
```

## Protein Structure Visualization

The protein structure is visualized interactively using NGLView:

```python
import nglview as nv

view = nv.show_biopython(
    structure,
    gui=True
)

view
```

This allows the user to rotate, zoom, and inspect the three-dimensional protein structure.

## Project Structure

```text
COVID19-Genome-Analysis/
│
├── covid19.ipynb
├── protein_seq.fasta
├── 9SAP.pdb
└── README.md
```



## Requirements

* Python 3.x
* Jupyter Notebook or Visual Studio Code
* Biopython
* Matplotlib
* NGLView
* Internet connection for:

  * NCBI sequence retrieval
  * BLAST searches
  * Downloading PDB files

## Outputs

The project produces the following outputs:

* Genome information and sequence statistics
* Molecular weight of the genome
* GC content
* Nucleotide distribution plot
* Amino acid distribution plot
* Protein sequences after stop-codon separation
* Filtered protein sequences of at least 20 amino acids
* FASTA file containing the selected protein sequence
* BLAST search results
* PDB structure file for `9SAP`
* Protein chain and residue information
* Interactive 3D protein structure visualization

## Technologies Used

* Python
* Biopython
* Matplotlib
* NCBI Entrez
* NCBI BLAST
* Protein Data Bank
* NGLView
* Jupyter Notebook

## Author

**Sunday Abdulsalam**
* Linkedin: https://www.linkedin.com/in/sundayabdulsalam/

* Email: abdulsalamsunday@yahoo.com

* Bioinformatics | Computational Biology | Data Analytics

Interested in:
•	Transcriptomics 
•	Cancer genomics 
•	Computational biology 
•	AI for health 
•	Data science 
•	Molecular biology 

# Disclaimer

The genome and protein data used in this project are obtained from publicly available biological databases and resources. The accuracy and completeness of results depend on the quality of the source data, selected parameters, software versions, and computational methods used.

This project is intended to demonstrate bioinformatics concepts such as genome analysis, sequence transcription and translation, amino acid analysis, protein sequence extraction, BLAST searches, PDB structure analysis, and molecular visualization.


