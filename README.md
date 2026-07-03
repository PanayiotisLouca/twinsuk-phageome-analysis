# TwinsUK Phageome Analysis

### Population-level characterisation of the human gut phageome and its associations with bacterial communities and microbial metabolites

*Panayiotis Louca, Mohammadali Khan Mirzaei, Afroditi Kouraki, Erfan Khamespanah, Yu Lin, Xue Peng, Robert Pope, Alessia Visconti, Francesco Asnicar, Daniel Kirk, Ricardo Costeira, Nicola Segata, Mario Falchi, Jordana T. Bell, Tim D. Spector, Lindsey A. Edwards, Li Deng, Ana M. Valdes, Cristina Menni*

---

This repository contains the R scripts used for the downstream statistical analyses presented in *Population-level characterisation of the human gut phageome and its associations with bacterial communities and microbial metabolites*.

## Citation
> *If you use this code, please cite:*
> > Louca, P. et al. (2025). *Population-level characterisation of the human gut phageome and its associations with bacterial communities and microbial metabolites*. [DOI: To be updated]

---

  ## 📂 Repository Structure
```
.
├── 1.Phage_Heritability
│   └── 1.1.Phage_alpha_div_heritability_SCRIPT.R
│   └── 1.2.Phage_beta_div_heritability_SCRIPT.R
│   └── 1.3.Individual_phage_heritability_SCRIPT_CREATE.R
├── 2.Phage_Taxa_Network_Analysis
│   ├── 2.1.Phage_species_ggraph_network_POSITIVE_SCRIPT.R
│   ├── 2.2.Phage_species_ggraph_network_NEGATIVE_SCRIPT.R
├── 3.Phage_Bacterial_Metabolite_Associations_and_Network
│   ├── 3.1.Phage_bile_acids_SCRIPT_CREATE_ARRAY.R
│   ├── 3.2.Phage_SCFA_SCRIPT_CREATE_ARRAY.R
│   └── 3.3.Phage_SCFA_bile_acid_network_SCRIPT.R
├── 4.Phage_Diet_Associations
│   └── 4.1.Phages_diet_SCRIPT_CREATE_ARRAY.R
└── README.md

```

---

## ℹ️ Repository Information

### 1. Phage Heritability

- `1.1.Phage_alpha_div_heritability_SCRIPT.R`: Estimates the heritability of phage alpha diversity (Shannon diversity and observed richness) using the `mets` package.
- `1.2.Phage_beta_div_heritability_SCRIPT.R`: Estimates the heritability of phage beta diversity using principal coordinates derived from Bray–Curtis dissimilarities.
- `1.3.Individual_phage_heritability_SCRIPT_CREATE.R`: Estimates the heritability of individual viral operational taxonomic units (vOTUs).

### 2. Phage Taxa Network Analysis

- `2.1.Phage_species_ggraph_network_POSITIVE_SCRIPT.R`: Constructs networks of positively correlated vOTUs & bacterial species.
- `2.2.Phage_species_ggraph_network_NEGATIVE_SCRIPT.R`: Constructs networks of negatively correlated vOTUs & bacterial species.

### 3. Phage - Bacterial Metabolite Associations and Network

- `3.1.Phage_bile_acids_CREATE_ARRAY_SCRIPT.R`: Tests associations between vOTUs and circulating/faecal bile acids.
- `3.2.Phage_SCFA_CREATE_ARRAY_SCRIPT.R`: Tests associations between vOTUs and short-chain fatty acids.
- `3.3.Phage_SCFA_bile_acid_network_SCRIPT.R`: Builds correlation networks between vOTUs, bile acids, and short-chain fatty acids.

### 4. Phage Diet Associations

- `4.1.Phages_diet_CREATE_ARRAY_SCRIPT.R`: Tests associations between vOTUs and dietary variables using linear mixed-effects models.

---

## Data Availability

The analyses in this repository were performed using data from the TwinsUK cohort. Access to the underlying datasets is managed by the Department of Twin Research at King's College London.

Data can be made available to bona fide researchers following the TwinsUK data access procedures:

https://twinsuk.ac.uk/resources-for-researchers/access-our-data/

The viral abundance tables analysed in this repository were generated using the [TwinsUK Viromics Profiling Pipeline](https://github.com/PanayiotisLouca/twinsuk-viromics-pipeline).

---

## License

This repository is distributed under the MIT License. See the `LICENSE` file for details.