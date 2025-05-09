# Linear Discriminant Analysis (LDA) for Species Delimitation

This repository contains an R script to conduct **Linear Discriminant Analysis (LDA)** for species delimitation using morphometric measurements of macromycetes. The workflow includes data standardization, dimensionality reduction, variable contribution analysis, Mahalanobis distance calculation, and high-quality visualizations.

---

## 📊 Overview

The script performs the following steps:

1. **Install and load required packages**
2. **Transform raw morphometric measurements into Z-scores**
3. **Run LDA and extract projected coordinates and variable contributions**
4. **Calculate Mahalanobis distances (individual and inter-group)**
5. **Optionally filter top contributing variables for plot clarity**
6. **Visualize group separation and variable influence (biplot and barplot)**
7. **Export summary tables and graphics**

---

## ⚙️ Why Standardize with Z-scores?

Micromorphological variables often have different measurement scales  
(e.g., *basidia length* in 80–100 µm vs. *spine length* in 0.8–1.2 µm).  
Z-score transformation allows:
- **Scaling all variables to a common range**
- **Avoiding numerical bias in multivariate analyses**
- **Combining features with contrasting units (e.g., macroscopic and microscopic)**

The script automatically handles this transformation and outputs a standardized `.csv`.

---

## 🔧 Usage

### 📂 Input File

Your input `.csv` must:
- Include one column named `Species`
- Contain only numeric morphometric variables in the remaining columns

### ⚙️ Customization Tips

- **Arrow scaling**: You can adjust `scale_factor` to modify arrow length in the biplot.
- **Variable filtering**: Use `top_vars` to show only variables in the top 25% of contribution to LD1 or LD2.

---

## 📂 Output

The script generates the following files:

- `Micromorphology_Measurements_Example_Z-Scores.csv`: Standardized measurements
- `VariableContribution_LDA_Example.csv`: Contribution of each variable to LD1 and LD2
- `Mahalanobis_Distances_LDA.csv`: Distance of each specimen to the centroid of its predicted group
- `Mahalanobis_CentroidDistances.csv`: Pairwise Mahalanobis distances between species centroids

In addition, it produces a **combined LDA figure** composed of:

- A **biplot** showing species separation in LDA space (LD1 vs LD2), including:
  - Colored ellipses per species
  - Morphometric variables represented as arrows scaled by their contribution
- A **barplot** showing the percentage contribution of each variable to LD1 and LD2 axes

These plots are rendered side by side in a single figure for visual comparison.

![Example_LDA+Contributions](https://github.com/user-attachments/assets/ac2091b2-7dfd-4320-848f-0002b7c17162)


---

## 📈 Interpretation of Mahalanobis Distances

- **Individual distances** (`Mahalanobis_Distances_LDA.csv`):  
  Indicate how far each specimen is from the centroid of its group in LDA space.  
  Higher values may suggest morphological outliers or weak group membership.

- **Centroid distances** (`Mahalanobis_CentroidDistances.csv`):  
  Represent multivariate morphological divergence between species.  
  Larger values reflect stronger group separation and may support taxonomic distinctiveness.

---

## 📘 Dependencies

The script installs and loads the following R packages:

- `ggplot2`
- `MASS`
- `dplyr`
- `tidyr`
- `gridExtra`
- `showtext`
- `sysfonts`
- `grid`

Fonts are loaded from Google Fonts using `Roboto` for a consistent, professional aesthetic.

---

## 🧪 Citation

If you use this script in your research, please cite:

**Cruz-Campuzano, E.A. et al. (2025)**  
*Jigsaw Falling into Place: species delimitation in Thelephora, with insights into a cryptic Neotropical species complex*  
_DOI: [to be added]_

---

## 📫 Contact

For questions, feedback, or contributions, please contact:  
**Ezequiel Alberto Cruz-Campuzano**  
ezequielcruz.1997@gmail.com

**Julieta Alvarez-Manjarrez, PhD** 
julieta.alvarez@ib.unam.mx 

---
