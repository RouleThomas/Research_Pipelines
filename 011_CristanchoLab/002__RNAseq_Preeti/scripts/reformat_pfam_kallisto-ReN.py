import csv

# Input and output file paths
input_file = "output/pfam/pfam_results_kallisto-ReN.txt"  # Replace with your input file name
output_file = "output/pfam/pfam_results_kallisto_reformat-ReN.txt"

# Open the input file and write the reformatted output
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.DictReader(infile)  # Reads input as dictionary
    outfile.write("# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>\n")

    for row in reader:
        # Format each row
        reformatted_row = f"{row['seq_id']:<20} {row['aln_start']:>4} {row['aln_end']:>8} {row['env_start']:>8} {row['env_end']:>8} {row['hmm_acc']:<12} {row['hmm_name']:<16} {row['type']:<10} {row['hmm_start']:>4} {row['hmm_end']:>4} {row['hmm_length']:>4} {row['score']:>10} {row['evalue']:>10} {row['significance']:>4} {row['clan']}\n"
        outfile.write(reformatted_row)

print("Reformatting complete. Output saved to:", output_file)
