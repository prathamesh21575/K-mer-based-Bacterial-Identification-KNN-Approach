# K-mer-based-Bacterial-Identification-KNN-Approach

This project introduces a k-mer based prediction model for bacterial identification, offering computational efficiency compared to traditional alignment techniques. It leverages Biopython and collections for k-mer matching and features an interactive Tkinter GUI for user-friendly DNA sequence comparison. The model uses a FASTA-formatted reference database and applies the K-Nearest Neighbors (KNN) algorithm to achieve high accuracy. Users can input query sequences, view the top matching sequences with their percentages, and access entire sequences through an intuitive interface.


Explanation of the Code:

This Python code creates a graphical user interface (GUI) for comparing a DNA sequence against a reference database using a technique called k-mer matching. The application uses the tkinter library for the GUI and BioPython for handling DNA sequences. Here's a step-by-step explanation of how the code works:

Imports and Global Variables:
The necessary libraries are imported: tkinter for the GUI, SeqIO from BioPython for reading DNA sequences, and some machine learning tools from scikit-learn.
A global dictionary, reference_sequences, is defined to store the reference DNA sequences.

K-mer Matching Function:
The kmer_matching function takes a query DNA sequence and compares it to reference sequences by breaking them into smaller fragments called k-mers (of length k).
It returns a dictionary where each reference sequence ID is mapped to the percentage of k-mers it shares with the query sequence and a list of these matching k-mers.

Comparison Function:
The compare_sequences function uses kmer_matching to compare the query sequence against all reference sequences.
It returns a list of the top 10 matches, including the sequence ID, description, match percentage, and top 1 or 2 matching k-mer fragments.

GUI Components:
Main Window: The main window contains an entry box for the user to input a DNA sequence and another entry box to specify the k-mer length.
Labels and Buttons: There are labels to display instructions and buttons to initiate the sequence comparison.

Event Handling Functions:
compare_sequences_button_click: This function is triggered when the "Compare sequences" button is clicked. It loads the reference sequences from a FASTA file, validates the query sequence, performs the k-mer matching, and displays the results.
show_comparison_results: This function displays the comparison results in a new window using a table format.
on_hit_double_click: This function is triggered when a user double-clicks on a result entry, showing the full DNA sequence of the selected match.
show_full_sequence: This function creates a new window to display the full sequence of a selected DNA record.

Utility Functions:
load_fasta_file: This function loads the reference DNA sequences from a FASTA file.
is_dna_sequence: This function checks if a given sequence is a valid DNA sequence (containing only A, C, G, or T).
remove_placeholder: This function removes placeholder text from the query entry box when the user clicks on it.

Main GUI Execution:
The main part of the code creates the GUI, including the main window, entry boxes for the query sequence and k-mer length, and buttons for comparing sequences.
When the user clicks the "Compare sequences" button, the application loads the reference sequences, validates the query sequence, performs the comparison, and displays the results in a new window.
