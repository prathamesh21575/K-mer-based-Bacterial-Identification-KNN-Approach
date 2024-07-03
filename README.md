# K-mer-based-Bacterial-Identification-KNN-Approach

This project introduces a k-mer based prediction model for bacterial identification, offering computational efficiency compared to traditional alignment techniques. It leverages Biopython and collections for k-mer matching and features an interactive Tkinter GUI for user-friendly DNA sequence comparison. The model uses a FASTA-formatted reference database and applies the K-Nearest Neighbors (KNN) algorithm to achieve high accuracy. Users can input query sequences, view the top matching sequences with their percentages, and access entire sequences through an intuitive interface.


Explanation of the Code:
This Python code creates a graphical user interface (GUI) for comparing a DNA sequence against a reference database using k-mer matching. The application uses tkinter for the GUI and BioPython for handling DNA sequences. Here are the main features and components:

Imports and Global Variables:

Imports necessary libraries: tkinter, SeqIO from BioPython, and some tools from scikit-learn.
Defines a global dictionary reference_sequences to store reference DNA sequences.
K-mer Matching Function:

kmer_matching: Compares a query DNA sequence with reference sequences by breaking them into smaller fragments called k-mers.
Returns a dictionary with the percentage of shared k-mers and the matching fragments for each reference sequence.
Comparison Function:

compare_sequences: Uses kmer_matching to find the top 10 matches from the reference sequences.
Returns a list of matches including sequence ID, description, match percentage, and top matching k-mer fragments.
GUI Components:

Main Window: Contains entry boxes for the DNA sequence and k-mer length.
Labels and Buttons: Display instructions and initiate the sequence comparison.
Event Handling Functions:

compare_sequences_button_click: Triggered when the "Compare sequences" button is clicked. Loads reference sequences, validates the query sequence, performs k-mer matching, and displays results.
show_comparison_results: Displays comparison results in a table format.
on_hit_double_click: Shows the full DNA sequence of a selected match when a result entry is double-clicked.
show_full_sequence: Creates a new window to display the full sequence of a selected DNA record.
Utility Functions:

load_fasta_file: Loads reference DNA sequences from a FASTA file.
is_dna_sequence: Checks if a given sequence is a valid DNA sequence (only A, C, G, T).
remove_placeholder: Removes placeholder text from the query entry box when clicked.
Main GUI Execution:

Creates the main window, entry boxes for the query sequence and k-mer length, and buttons for comparing sequences.
When the "Compare sequences" button is clicked, it loads reference sequences, validates the query, performs the comparison, and displays the results in a new window.
