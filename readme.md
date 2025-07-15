# JCAP CRISPR Mixscape Pipeline

Welcome to the **JCAP CRISPR Mixscape Pipeline**, a user-friendly Shiny application for interactive single-cell CRISPR analysis using the [Mixscape](https://satijalab.org/seurat/articles/mixscape_vignette.html) workflow (Seurat). This app enables wet-lab biologists and computational researchers to analyze CRISPR perturbation screens with visual and statistical clarity.

---

## 🌱 Features

- **Upload your own count matrix and metadata**
- **Automated QC, normalization, HVG selection**
- **UMAP visualization**
- **Mixscape perturbation score calculation**
- **KO/NP/NT class assignment**
- **Interactive plots (barplots, violin plots, heatmaps, etc)**
- **Downloadable differential gene expression (KO vs NT) tables**
- **Summary statistics for your dataset**
- **Responsive design, with soft pastel CSS inspired by _Les Carnets de l’Apothicaire_**

---

## 🚀 Quick Start

### 1. Install R dependencies

```r
# In R/RStudio, run:
install.packages(c("shiny", "Seurat", "dplyr", "ggplot2", "patchwork", "plotly", "DT"))
# For most, CRAN is fine. If you use advanced Seurat features, consider Bioconductor.
2. Clone this repository
bash
Copy
Edit
git clone https://github.com/yourusername/jcap-crispr-mixscape.git
cd jcap-crispr-mixscape
3. Run the app
r
Copy
Edit
shiny::runApp()   # or specify the path if outside the directory
📂 Input Data Format
Counts matrix: CSV, cells as columns, genes as rows.

Metadata: CSV, one row per cell, with required columns: gene, replicate, guide_ID, etc.

Example files are provided in the example_data/ folder.

🖱️ Using the App
Upload your count matrix and metadata CSVs.

Click Run Normalization & UMAP to process the data.

Adjust number of neighbors if desired.

Click Run Mixscape Analysis.

Explore the results in the different tabs:

UMAP Plots

Summary Table

KO Genes Table (downloadable)

KO % Barplots

Perturbation Score

Posterior Violin

Heatmap

How To instructions

Download gene tables and summary stats as needed.

For detailed step-by-step instructions, see the How To tab inside the app.

📊 Example Data
You’ll find example data in example_data/ to help you get started.
Note: For real Mixscape results (KO/NP/NT separation, DGE, etc.), use a dataset with at least ~100 cells per group and >500 variable genes.

🎨 Customization
Theme/CSS: The app features a pastel garden theme (see www/custom.css).

Modular Code: Feel free to fork and adapt for your own CRISPR or single-cell workflows!

📖 Documentation
See How_To.txt or the in-app “How To” tab for a walkthrough of all features.

🙏 Credits
Original method: Seurat Mixscape vignette

UI inspiration: Les Carnets de l’Apothicaire (Lumen editions art)

App author: Your Name Here

Contributors: Open to pull requests!

🛟 Support
Open a GitHub Issue for bugs or feature requests
