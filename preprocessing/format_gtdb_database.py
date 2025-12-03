
"""
Formats the headers of gzipped GTDB database 16S gene FASTA files ready for DADA2/Phyloseq. Without this minor formatting, the taxonomic assignment will fail.

The input files originate from here: https://figshare.scilifelab.se/articles/dataset/SBDI_Sativa_curated_16S_GTDB_database/14869077/9
 - The taxonomy file is used for assignTaxonomy
    - We use the '1 genome' file, as recommended by the publishers
 - The species file is used for addSpecies
    - We use the '20 genomes' file, as recommended by the publishers

Example of input TAXONOMY_TRAIN_SET_UNFORMATTED header:
>Archaea;Archaea;Methanobacteriota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobacterium_A;Methanobacterium_A sp002494495

Example of input SPECIES_ASSIGNMENT_SET_UNFORMATTED header:
>RS_GCF_000189915.1~NZ_AELQ01000009.1 Methanocatella smithii_A

Formatting rules:
- Taxonomy headers (assignTaxonomy):
  - Input is expected to be 8 semicolon-delimited ranks.
  - Drop the first rank and the last rank.
  - Join with semicolons and append a trailing semicolon. Output must therefore contain exactly 6 semicolons and end with `;`.
- Species headers (addSpecies):
  - Nothing to do here.
"""


# Imports
import gzip
import os


# Input files
TAXONOMY_TRAIN_SET_UNFORMATTED = "/mnt/secondary/16S_databases/sbdi-gtdb-sativa.r10rs226.1genome.assignTaxonomy.fna.gz"    #  <--- GTDB r226 (https://figshare.scilifelab.se/articles/dataset/SBDI_Sativa_curated_16S_GTDB_database/14869077/9)
SPECIES_ASSIGNMENT_SET_UNFORMATTED = "/mnt/secondary/16S_databases/sbdi-gtdb-sativa.r10rs226.20genomes.addSpecies.fna.gz"    #  <--- GTDB r226
# Output files
TAXONOMY_TRAIN_SET_FORMATTED = "/mnt/secondary/16S_databases/sbdi-gtdb-sativa.r10rs226.1genome.assignTaxonomy_formatted.fna.gz"
SPECIES_ASSIGNMENT_SET_FORMATTED = "/mnt/secondary/16S_databases/sbdi-gtdb-sativa.r10rs226.20genomes.addSpecies_formatted.fna.gz"


def process_tax_file():
    """
    Processes the taxonomy FASTA file.
    - Formats headers by removing the first and last labels.
    - Gathers statistics on header label counts.
    """
    print(f"Processing taxonomy file: {TAXONOMY_TRAIN_SET_UNFORMATTED}")
    
    headers_processed = 0
    incorrect_label_count = 0
    mismatched_first_labels_count = 0

    try:
        with gzip.open(TAXONOMY_TRAIN_SET_UNFORMATTED, 'rt') as infile, gzip.open(TAXONOMY_TRAIN_SET_FORMATTED, 'wt') as outfile:
            for line in infile:
                if line.startswith('>'):
                    headers_processed += 1
                    header = line.strip()[1:]
                    labels = header.split(';')
                    
                    if len(labels) != 8:
                        incorrect_label_count += 1

                    if len(labels) >= 2 and labels[0] != labels[1]:
                        mismatched_first_labels_count += 1
                    
                    if len(labels) > 1:
                        new_labels = labels[1:-1]
                        new_header = '>' + ';'.join(new_labels) + ';\n'
                        outfile.write(new_header)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)

        print("Taxonomy file processing complete.")
        print("\n--- Report for Taxonomy File ---")
        print(f"Total headers processed: {headers_processed}")
        print(f"Headers with number of labels not equal to 8: {incorrect_label_count}")
        print(f"Headers where the first two labels are not the same: {mismatched_first_labels_count}")
        print("--------------------------------\n")

    except FileNotFoundError:
        print(f"Error: Input file not found at {TAXONOMY_TRAIN_SET_UNFORMATTED}")
    except Exception as e:
        print(f"An error occurred during taxonomy file processing: {e}")

def process_species_file():
    """
    Processes the species FASTA file.
    - Reads and writes the file without modification to maintain consistent processing style.
    """
    print(f"Processing species file: {SPECIES_ASSIGNMENT_SET_UNFORMATTED}")
    try:
        with gzip.open(SPECIES_ASSIGNMENT_SET_UNFORMATTED, 'rt') as infile, gzip.open(SPECIES_ASSIGNMENT_SET_FORMATTED, 'wt') as outfile:
            for line in infile:
                outfile.write(line)
        print("Species file processing complete.")
    except FileNotFoundError:
        print(f"Error: Input file not found at {SPECIES_ASSIGNMENT_SET_UNFORMATTED}")
    except Exception as e:
        print(f"An error occurred during species file processing: {e}")

if __name__ == "__main__":

    # Check if the output files already exist
    if os.path.exists(TAXONOMY_TRAIN_SET_FORMATTED) and os.path.exists(SPECIES_ASSIGNMENT_SET_FORMATTED):
        input("Output files already exist. Press enter to overwrite them.")
        os.remove(TAXONOMY_TRAIN_SET_FORMATTED)
        os.remove(SPECIES_ASSIGNMENT_SET_FORMATTED)

    process_tax_file()
    process_species_file()
    print("\nAll files have been processed.")


