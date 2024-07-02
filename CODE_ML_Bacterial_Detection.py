import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from Bio import SeqIO
from collections import defaultdict
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import KNeighborsClassifier
import pandas as pd

reference_sequences = {}

def kmer_matching(sequence, reference_sequences, k=25):
    """
    Calculates the number of shared k-mers between a query sequence and reference sequences.

    Args:
        sequence: The query sequence (DNA string).
        reference_sequences: A dictionary mapping sequence IDs to their DNA sequences.
        k: The k-mer size (number of nucleotides).

    Returns:
        A dictionary mapping sequence IDs to a tuple containing the percentage of shared k-mers with the query sequence
        and a list of matching k-mer fragments.
    """
    matches = defaultdict(lambda: [0, []])
    input_kmers = set([sequence[i:i + k] for i in range(len(sequence) - k + 1)])
    query_kmer_count = len(input_kmers)

    for seq_id, ref_seqrecord in reference_sequences.items():
        # Extract the DNA sequence for k-mer calculation
        ref_sequence = ref_seqrecord.seq
        ref_kmers = set([ref_sequence[i:i + k] for i in range(len(ref_sequence) - k + 1)])
        common_kmers = input_kmers.intersection(ref_kmers)
        matches[seq_id][0] = (len(common_kmers) / query_kmer_count) * 100
        matches[seq_id][1] = list(common_kmers)

    return matches

def compare_sequences(query_sequence, reference_sequences, k=25):
    """
    Compares a query sequence to the reference database using k-mer matching.

    Args:
        query_sequence: The query sequence (DNA string).
        reference_sequences: A dictionary mapping sequence IDs to their DNA sequences.
        k: The k-mer size (number of nucleotides).

    Returns:
        A list of tuples containing ID, name, % match, and top 1 or 2 matching k-mer fragments of up to 10 top matching sequences from the reference database.
    """
    matches = kmer_matching(query_sequence, reference_sequences, k)
    sorted_matches = sorted(matches.items(), key=lambda x: x[1][0], reverse=True)
    top_matches = []
    for seq_id, match_data in sorted_matches[:10]:
        seq_description = reference_sequences[seq_id].description
        match_percentage = match_data[0]
        matching_fragments = match_data[1][:2]  # Extract top 1 or 2 matching fragments
        top_matches.append((seq_id, seq_description, match_percentage, matching_fragments))
    return top_matches

def compare_sequences_button_click():
    # Load the FASTA file
    load_fasta_file()

    query_sequence = entry_query.get().upper()
    if not query_sequence:
        messagebox.showerror("Error", "Please enter a query sequence.")
        return

    if not is_dna_sequence(query_sequence):
        messagebox.showerror("Error", "Invalid DNA sequence. Please enter a sequence with only A, C, G, or T.")
        return

    kmer_size = entry_kmer.get()
    if not kmer_size:
        kmer_size = 25
    else:
        kmer_size = int(kmer_size)

    top_matches = compare_sequences(query_sequence, reference_sequences, kmer_size)
    show_comparison_results(top_matches)

def show_comparison_results(top_matches):
    comparison_window = tk.Toplevel(root)
    comparison_window.title("Comparison Results")
    # Set geometry to full-screen
    width = comparison_window.winfo_screenwidth()
    height = comparison_window.winfo_screenheight()
    comparison_window.geometry("%dx%d+0+0" % (width, height))

    tree = ttk.Treeview(comparison_window, columns=('Sr. No.', 'ID', 'Name', 'Match%', 'Kmer', 'Position'), show='headings')
    tree.heading('Sr. No.', text='Sr. No.', anchor=tk.CENTER)
    tree.heading('ID', text='ID', anchor=tk.CENTER)
    tree.heading('Name', text='Name')
    tree.heading('Match%', text='Match %', anchor=tk.CENTER)
    tree.heading('Kmer', text='Kmer')
    tree.heading('Position', text='Position', anchor=tk.CENTER)
    tree.pack(fill='both', expand=True)

    # Adjust column widths
    tree.column('Sr. No.', width=10, anchor=tk.CENTER)
    tree.column('ID', width=50, anchor=tk.CENTER)
    tree.column('Name', width=400)
    tree.column('Match%', width=10, anchor=tk.CENTER)
    tree.column('Kmer', width=500)
    tree.column('Position', width=10, anchor=tk.CENTER)

    for i, (seq_id, seq_description, match_percentage, matching_fragments) in enumerate(top_matches, 1):
        if match_percentage == 0.00:
            continue  # Skip displaying this entry

        # Insert blank row to increase distance between hits (ID)
        tree.insert('', 'end', values=('', '', '', '', '', ''))

        # Display each hit once with its top two matching k-mers
        tree.insert('', 'end', values=(i, seq_id, seq_description, f"{match_percentage:.2f}%", matching_fragments[0], reference_sequences[seq_id].seq.find(matching_fragments[0])), tags=(seq_id,))
        tree.insert('', 'end', values=('', '', '', '', matching_fragments[1], reference_sequences[seq_id].seq.find(matching_fragments[1])), tags=(seq_id,))

    # Bind double click event to show full sequence
    tree.bind("<Double-1>", on_hit_double_click)

def on_hit_double_click(event):
    item = event.widget.focus()
    values = event.widget.item(item, "values")
    if values:
        seq_id = values[1]
        seq_description = values[2]
        sequence_record = reference_sequences.get(seq_id)
        if sequence_record:
            show_full_sequence(seq_description, sequence_record)
        else:
            messagebox.showerror("Error", "Sequence record not found.")
    else:
        messagebox.showerror("Error", "No values associated with the selected item.")

def show_full_sequence(seq_description, sequence_record):
    sequence_window = tk.Toplevel(root)
    sequence_window.title("Full Sequence")
    # Set geometry to full-screen
    width = sequence_window.winfo_screenwidth()
    height = sequence_window.winfo_screenheight()
    sequence_window.geometry("%dx%d+0+0" % (width, height))

    sequence_text_widget = scrolledtext.ScrolledText(sequence_window, wrap=tk.WORD, width=60, height=20, font=("Arial", 12))
    sequence_text_widget.pack(expand=True, fill=tk.BOTH)
    sequence_text_widget.insert(tk.END, f"Name: {seq_description}\n", "bold")  # Display only the name
    sequence_text_widget.insert(tk.END, f"ID: {sequence_record.id}\n")
    sequence_text_widget.insert(tk.END, f"Sequence:\n{sequence_record.seq}")
    sequence_text_widget.tag_configure("bold", font=("Arial", 14, "bold"))
    sequence_text_widget.configure(state=tk.DISABLED)

def load_fasta_file():
    global reference_sequences
    file_path = r"C:\Users\prath\Desktop\freshwater chip\seq_EX.fasta" # Provide the file path here
    try:
        reference_sequences = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
        # messagebox.showinfo("Success", "FASTA file loaded successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"Error loading FASTA file: {e}")

def is_dna_sequence(sequence):
    """
    Checks if a sequence contains only valid DNA characters (A, C, G, T).

    Args:
        sequence: The DNA sequence to check.

    Returns:
        True if the sequence contains only valid DNA characters, False otherwise.
    """
    return all(base in ('A', 'C', 'G', 'T') for base in sequence)

def remove_placeholder(event):
    if entry_query.get() == "Enter DNA sequence here":
        entry_query.delete(0, tk.END)

if __name__ == '__main__':
    root = tk.Tk()
    root.title("DNA Sequence Comparison")

    # Add a Label widget to display the homepage content
    homepage_text = """
  <<Welcome to DNA Sequence Comparison Tool !!!>>

This tool is designed to compare a query DNA sequence against a reference database of DNA sequences.
It utilizes k-mer matching, a powerful technique for sequence analysis, to identify similarities between sequences.

<< Key features >>
	Fast and efficient comparison using k-mer matching algorithm.
	Multi-threaded processing for parallel comparison of large datasets.
	User-friendly interface for easy input and visualization of results.

<< How k-mers are selected and calculated >>
	K-mers are short subsequences of length 'k' extracted from DNA sequences.
	The tool calculates the number of shared k-mers between the query sequence and reference sequences.
	The percentage of shared k-mers is used to determine the similarity between sequences.

Please upload your reference sequences in FASTA format and provide a query sequence to start comparing.

<< Thank you for using our tool >>

    """
    homepage_label = tk.Label(root, text=homepage_text, justify='left', padx=20, pady=20, font=('Arial', 12))
    homepage_label.pack()

    # Entry for entering query sequence
    entry_query = ttk.Entry(root, width=100)
    entry_query.pack(pady=30)
    entry_query.insert(0, "Enter DNA sequence here")
    entry_query.bind("<FocusIn>", remove_placeholder)

    # Entry for entering k-mer size
    kmer_frame = ttk.Frame(root)
    kmer_frame.pack()
    kmer_label = ttk.Label(kmer_frame, text="Enter k-mer length:")
    kmer_label.pack(side=tk.LEFT)
    entry_kmer = ttk.Entry(kmer_frame, width=15)
    entry_kmer.pack(side=tk.LEFT)

    # Compare sequences button
    search_button = ttk.Button(root, text="Compare sequences", command=compare_sequences_button_click, width=50)
    search_button.pack(pady=10)

    # Display the GUI
    root.mainloop()