# Micromorphological Character Comparison in Macromycetes

This repository contains an R script for performing statistical analysis and Boxplot+Violin for data visualization of micromorphological characters in macromycetes. The workflow supports both **parametric** and **non-parametric** approaches depending on the assumptions met by the input data.

---

## 📊 Overview

The script performs the following steps:

1. **Load data and required packages**
2. **Test for normality (Shapiro-Wilk) and homogeneity of variances (Levene’s test)**
3. **Automatically suggest the appropriate statistical test**
4. **Run either ANOVA + Tukey or Kruskal-Wallis + Dunn tests**
5. **Visualize the results with violin and box plots**
6. **Annotate significance levels on the plots**
7. **Export plots and statistical results**

---

## 🔧 Usage

### Input File
Provide your data as a `.csv` file where:
- One column must be named `Species`
- Other columns should be numeric variables representing micromorphological characters (e.g., `Spore.length`, `Basidia.width`, etc.)

### Default Configuration
The script is set to run **non-parametric analyses** (Kruskal-Wallis + Dunn) by default.  
To switch to **parametric analysis** (ANOVA + Tukey), follow the instructions and comments in:

- `STEP 4`: Modify `parametric = TRUE` for relevant variables
- `create_violin_plot()` function

---

## 📂 Output

The script produces:

- A combined figure with all plots (`cowplot::plot_grid`)
- CSV files:
  - `Diagnostic_Summary.csv`: Normality and homogeneity test results
  - `Global_Tests.csv`: ANOVA or Kruskal-Wallis p-values
  - `Multiple_Comparisons.csv`: Pairwise test results (Tukey or Dunn)

Plots are annotated with:
- Means
- Significant group comparisons
- Test statistics

---

## 📘 Dependencies

Required R packages:
- `ggplot2`
- `dplyr`
- `rstatix`
- `car`
- `cowplot`
- `showtext`
- `sysfonts`

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
**Ezequiel Alberto Cruz-Campuzano**  
ezequielcruz.1997@gmail.com

**Julieta Alvarez-Manjarrez, PhD** 
julieta.alvarez@ib.unam.mx 

---
