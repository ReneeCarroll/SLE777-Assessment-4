# SLE777-Assessment-4
Bioinformatics Assessment 4 Project File
Deakin University

## Author: Renee Carroll
## Deakin ID: s225723235

# Summary

This is the README document for SLE777 Applied Bioinformatics' Assessment 4.

This Assessment contains 2 Parts.

### Part 1

Part 1 is the exploration of two different data sets.
Questions 1-5 are focused on analyzing gene expression data from a .tsv file
Questions 6-10 are focused on analyzing tree growth data from a .csv file.

### Part 2

Part 2 is a deeper exploration of biological sequence diversity between two organisms, E.coli, and separate organism allocated specifically to me, namely, Acetobacter aceti (A.aceti, GCA_002005445).
This part comprises 6 sections, each with a number of questions to guide the data exploration.

### Packages required - Summary

The packages required for this Assessment are:

"car"
"seqinr"
"R.utils"
"ggplot2"
"ggpubr"

The first code chunk in the RMarkdown document provides the code to install and run these packages, albeit with the install code inactivated by #.
Should the user require these packages installed, simply remove the # in front of the install code to activate this code again.

### Data downloads required - Summary

Within the R.Markdown document, there is code to download files from the 4 following URLs.

URL1 = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"
URL2 = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"
URL3 = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
URL4 = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_80_collection/acetobacter_aceti_gca_002005445/cds/Acetobacter_aceti_gca_002005445.ASM200544v1.cds.all.fa.gz"

## Please note: In section 2.1.1, the code to download an unzip the bacteria sequence files has been supressed with #, as it was inhibiting R.Markdown's Knit function during project development.

Should the user not already have this data downloaded and available, the code will need to be activated again by removing the # prior to running the R.Markdown file.

# All Scripts, with Explanations

























## License

MIT License

## Citations

"car" package: Fox, J. & Weisberg, S. (2019) An R Companion to Applied Regression, 3rd edn. Sage, Thousand Oaks, CA. Available at: https://www.john-fox.ca/Companion/
 (Accessed: 27 September 2025) 

"seqinr" package: Charif, D. & Lobry, J. (2007) “SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.” In Bastolla, U., Porto, M., Roman, H. & Vendruscolo, M. (eds.) Structural approaches to sequence evolution: Molecules, networks, populations, Biological and Medical Physics, Biomedical Engineering series, Springer Verlag, New York, pp. 207–232. Available at: https://cran.r-project.org/web/packages/seqinr/citation.html
 (Accessed: 25 September 2025) 

"R.utils" package: R.utils: Various programming utilities for R [Computer software]. Available at: https://CRAN.R-project.org/package=R.utils
 (Accessed: 25 September 2025)

"ggplot2" package: Wickham, H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York. Available at: https://ggplot2.tidyverse.org
 (Accessed: 01 October 2025) 

"ggpubr" package: Kassambara, A. (2025) ggpubr: ’ggplot2’ Based Publication Ready Plots [R package version 0.6.1]. Available at: https://rpkgs.datanovia.com/ggpubr/
 (Accessed: 01 October 2025)