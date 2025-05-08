# DNA Sequence Similarity Analysis across Species

This repository contains an R script to visualize genetic similarity among species based on pairwise comparison matrices (e.g., ITS, LSU, or other loci). The output is a density plot focused on a **focal species**, allowing clear visualization of intra- and interspecific variation (e.g., a candidate new species vs. similar taxa).

---

## 📊 Overview

The script performs the following steps:

1. **Load similarity matrix and required packages**
2. **Reshape the matrix into long format**
3. **Extract species names from specimen labels**
4. **Filter pairwise comparisons involving a focal species**
5. **Generate a density plot showing intra- and interspecific similarities**
6. **Customize axis scaling and styling for publication-ready output**

---

## 🔧 Usage

### Input File
Provide a `.csv` file containing a **square similarity matrix**, where:
- The **first column** includes specimen names (e.g., `"Species A.001"`)
- The **column headers** match the row labels
- **Cell values** represent similarity in percentage (%)

> 💡 The script automatically extracts species names by removing everything after the first period (`.`) in each specimen label (e.g., `"Species A.001"` → `"Species A"`).

### Optional Renaming
If your species names are long (e.g., `"Thelephora renispora"`), you can optionally abbreviate them by uncommenting and customizing the lines in **STEP 3**.

---

## ⚙️ Customization

- Set your **focal species** in `STEP 4` using the variable:  
  focal_species <- "Species A"
- Modify the legend colors in STEP 5: 
  color_palette <- c("Species A" = "#9933FF", ...)
- Adjust X-axis breaks in STEP 7 to match your expected similarity range:
  breaks = seq(90, 100, by = 2.5)  # Adjust based on focal species' lowest interspecific score

## 📂 Output

The script produces:

- A density plot with custom colors per species with the following:
  - Dashed lines at thresholds or breaks set manually (e.g. 92.5%, 95.0%, 97.5%)
  - Automatic axis scaling based on your data
  - Clean layout using the Roboto font

The output can be further customized to:
- Export plots as .png or .pdf
- Combine multiple genes or loci in a composite figure

---

## 📘 Dependencies

Required R packages:
- `ggplot2`
- `dplyr
- `tidyr`
- `stringr`
- `readr`
- `showtext`
- `sysfonts`
- `scales`

The script automatically installs any missing packages on first run.

---

## 🧪 Citation

If you use this script in your research, please cite:

**Cruz-Campuzano, E.A. et al. (2025)**  
*Jigsaw Falling into Place: species delimitation in Thelephora, with insights into a cryptic Neotropical species complex*  
_DOI: [to be added]_

---

## 📫 Contact

For questions, feedback, or contributions, please contact:  
**Ezequiel Alberto Cruz Campuzano**  
ezequielcruz.1997@gmail.com

**Julieta Alvarez Manjarrez** 
julieta.alvarez@ib.unam.mx 

---
