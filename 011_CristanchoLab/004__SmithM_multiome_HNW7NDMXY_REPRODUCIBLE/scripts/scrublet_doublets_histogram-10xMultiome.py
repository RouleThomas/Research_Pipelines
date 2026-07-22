import sys
import pandas as pd
from scipy.io import mmread
import scrublet as scr
import matplotlib.pyplot as plt

def run_scrublet(cellranger_path, output_path, manual_threshold=None):
    
    matrix = mmread(f"{cellranger_path}/matrix.mtx.gz").T.tocsr()
    barcodes = pd.read_csv(f"{cellranger_path}/barcodes.tsv.gz", header=None, sep="\t")[0]

    features = pd.read_csv(
        f"{cellranger_path}/features.tsv.gz",
        header=None,
        sep="\t"
    )

    # Keep only RNA / Gene Expression features
    gene_mask = features[2] == "Gene Expression"
    matrix_rna = matrix[:, gene_mask.values]

    print(f"Total features: {matrix.shape[1]}")
    print(f"Gene Expression features used for Scrublet: {matrix_rna.shape[1]}")
    print(f"Cells: {matrix_rna.shape[0]}")

    scrub = scr.Scrublet(matrix_rna)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    if manual_threshold is not None:
        predicted_doublets = doublet_scores > manual_threshold

    results = pd.DataFrame({
        "Doublet_score": doublet_scores,
        "Is_doublet": predicted_doublets
    }, index=barcodes)

    results.to_csv(output_path, sep="\t")

    plt.figure()
    plt.hist(doublet_scores, bins=50)
    plt.xlabel("Doublet Score")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Doublet Scores")

    if manual_threshold is not None:
        plt.axvline(manual_threshold, linestyle="--")

    plt.savefig(output_path.replace(".tsv", "_doublet_scores_histogram.png"))
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) not in [3, 4]:
        print("Usage: python3 script.py <input_path> <output_path> [manual_threshold]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    manual_threshold = float(sys.argv[3]) if len(sys.argv) == 4 else None

    run_scrublet(input_path, output_path, manual_threshold)