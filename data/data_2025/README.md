
This folder contains cleaned and raw datasets used in the study:  
**"Evolution-Based Deep Generative Design of SH3 Domains"**  
📄 Published in *Cell Systems* (2024)  
👉 [Link to paper](https://www.cell.com/cell-systems/abstract/S2405-4712(24)00204-7)

---

## 📁 File Overview

| File                         | Description                                   |
| ---------------------------- | --------------------------------------------- |
| `LibraryDesign.csv`          | Raw design library data (with NaNs)           |
| `LibraryDesign_cleaned.csv`  | Cleaned version (NaNs in `[RE_norm]` removed) |
| `LibraryNatural.csv`         | Raw natural sequence data (with NaNs)         |
| `LibraryNatural_cleaned.csv` | Cleaned version (NaNs in `[RE_norm]` removed) |
| `dropnan.py`                 | Python script used to clean the `.csv` files  |

Use the `_cleaned.csv` files for analysis. These have missing values removed from the key column `[RE_norm]`, a functional fitness proxy.

---

## 🧪 What is `[RE_norm]`?

`[RE_norm]` stands for normalized relative enrichment — a fitness-like score indicating how well a specific SH3 sequence performs its biological function in a cellular context.

In the paper, it reflects how each sequence complements yeast growth or signaling, giving insight into protein functionality.

## 🔹 What is `Sequence_aligned` and `Sequences_unaligned` in each `.csv`?
`Sequence_aligned` contains sequences that were both aligned and trimmed using methods described in the paper, making homologous positions directly comparable across sequences. This processing is essential for position-specific analyses, such as mapping `[RE_norm]` values to individual amino acid sites. In contrast, `Sequences_unaligned` are the original raw sequences, provided for reference and comparison.
## 📊 Project Ideas

#### 1. **Draw a heatmap of `[RE_norm]` vs amino acid positions**
- Use `LibraryDesign_cleaned.csv`
- Visualize fitness variation across positions
- Learn: plotting, data wrangling, amino acid indexing

#### 2. **t-SNE or PCA visualization of protein sequences**
- Use `LibraryNatural_cleaned.csv`
- Compare sequence embeddings (e.g., k-mer counts or pretrained models like ESM) across orthologous groups
- Color by taxonomy or functional group

#### 3. **Compare the fitness distributions**
- Plot `[RE_norm]` distributions for `LibraryDesign_cleaned.csv` vs `LibraryNatural_cleaned.csv`

#### 4. **Track conservation vs function**
- Identify conserved positions among functional proteins
- Highlight positions where mutations lead to large drops in `[RE_norm]`

#### 5. **Protein logo plots**
- Create sequence logos of high-performing proteins
- Compare with low-performing variants

---

## 🔧 Tools You Can Use

- Python (Jupyter Notebooks)
- pandas
- seaborn / matplotlib
- scikit-learn (for PCA/t-SNE)
- biopython (for sequence handling)

---

## 💡 Want more ideas?

- **Build a sequence classifier**: Train a simple machine learning model to predict `[RE_norm]` from sequence.
- **Phylogenetic clustering**: Group sequences based on similarity and see if they match functional clusters.
- **Mutation scanning**: Simulate or analyze effects of single amino acid mutations.

