# WBS Use-Resistance Modeling

R code for modeling the relationship between outpatient antibiotic use and antibiotic resistance gene (ARG) abundance in wastewater, as described in Brown et al. 2024: https://www.medrxiv.org/content/10.1101/2024.12.11.24318846v1

## Demo Script

`mag_seasonality_modeling_demo.R` is a fully self-contained script that:

1. **Generates synthetic data** — 52 weeks of simulated MAG abundance (RPKM), antibiotic use (SXT, folate antagonists), influent flow, crAssphage abundance, and day-of-week, all with realistic seasonal structure
2. **Fits Models 2, 3, and 4** across three simulated *Pseudomonas* MAGs (`gs43_SIM`, `gs77_SIM`, `gs12_SIM`)
3. **Compares models** via AIC
4. **Plots seasonal deviate panels** — observed vs. model-predicted deviates with cosine fit and amplitude confidence ribbon, arranged as a multi-panel figure via `ggarrange`

should take less than 1 minute if all packages are installed. 

### Dependencies

Install all required packages from CRAN:

```r
install.packages(c(
  "tidyr", "dplyr", "data.table", "plyr",
  "ggplot2", "ggpubr", "ggthemes", "purrr", "magrittr"
))
```

### Running the demo

```r
source("mag_seasonality_modeling_demo.R")
```

No input files are required — all data are generated within the script. However, the scripts used to run the real data are also provided (plotting_pseudomonas_modeling_results.R and pseudomonas-MAG-resistance-modeling.R). 

---

## License

MIT License — see [`LICENSE`](LICENSE) for details.
