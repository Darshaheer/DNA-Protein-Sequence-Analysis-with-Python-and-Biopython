# DNA-Protein-Sequence-Analysis-with-Python-and-Biopython
This project provides a simple pipeline to analyze DNA sequences and translate them into protein sequences using **Python** and the **Biopython** library.  
It also supports visualization of sequence features and prepares the output for further downstream analysis (e.g., BLAST).  
## 📂 Repository Structure
DNA-Protein-Sequence-Analysis-with-Python-and-Biopython/
│── sequence_analysis.py # Main analysis script
│── data/ # Input folder
│ └── Fasta_Seq.fasta # Example FASTA file
│── results/ # Output folder
│ └── (generated plots and analysis results)
│── README.md # Project documentation
│── LICENSE # License file (if included)
│── .gitignore # Ignore unnecessary files

## Features

- Reads **FASTA** files  
- Performs **basic DNA analysis** (length, nucleotide composition)  
- Transcribes DNA → mRNA  
- Translates DNA → Protein sequence  
- Provides **reverse complement**  
- Saves results in the `results/` folder for easy access  

---

## Requirements

- Python 3.13.5  
- [Biopython](https://biopython.org/)  

Install dependencies with:

```bash
pip install biopython matplotlib
``` 
  
## Usage
Place your DNA sequences in data/Fasta_Seq.fasta
Run the script:
python sequence_analysis.py
Output files (translated protein, summary, plots) will be saved in results/

## Example Output
Nucleotide composition statistics
Transcribed RNA sequence
Protein sequence from translation
Reverse complement of the DNA strand
Visual plots (saved in results/)

## Notes
Replace the sample FASTA file with your own sequence(s).
You can analyze multiple sequences at once if they are stored in FASTA format.
This project can be extended to integrate with NCBI BLAST for homology searches.

## License
MIT License
This project is open-source. You are free to use and modify it for research or educational purposes.
