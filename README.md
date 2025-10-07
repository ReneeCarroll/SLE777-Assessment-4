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

# Part 1
## Analysis of RNA-seq Expression Data

#### 1.1. Read in the file, making the gene identifiers the row names. Show a table of values for the first six genes.

To begin the analysis, the RNA-seq expression data file must first be read into R. 

The following code will download the required data file, and will read it into R. 
The code instructs that the first row of data contains the column headings, and the first column is to be used as row names. 
It also specifies that the column separation is by tabs, and that all text columns should stay as character data, and not be converted to categorical variables.

```{r download files from GitHub, and show first 6 genes}
URL1 = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"
download.file(URL1, destfile="gene_expression.tsv")
gene_expression <- read.delim("gene_expression.tsv", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
```

The gene expression data for the first six genes for each of the three samples of interest is shown in the table below. Using the 'knitr::kable()' code ensures the document output is in a nice table format, with caption included.

```{r first 6 genes}
knitr::kable(head(gene_expression),
             caption = "Table 1. Gene expression levles for the first six genes")
```

#### 1.2. Make a new column which is the mean of the other columns. Show a table of values for the first six genes.

A new column can easily be created using the $ code. The new column is added to the end of the data frame.
To this new column, we can save the row mean calculation of the expression values from the three samples, for each gene.
The head() command displays a table of the first 6 genes from the now updated data frame.

```{r mean saved as new column}
gene_expression$mean_expression <- rowMeans(gene_expression)
knitr::kable(head(gene_expression),
             caption = "Table 2. Gene expression levles for the first six genes, including mean values")
```

#### 1.3. List the 10 genes with the highest mean expression

In order to determine the highest mean expression, the data first must be put in order according to the column that contains the mean expression values.
By including the - sign, it is ensured that the data is ordered from largest to smallest, and that is is ordering based on the column 'mean_expression'.
The comma followed by the space prior to the closing square bracket instructs that all the columns are to be kept.

With the data now ordered, the top 10 rows containing the genes with the highest mean expression can be subsetted into a new data frame.
These genes can then be listed by retrieving the row names of this subsetted data frame, in line with the text, using 'r (code)' inside the specialised R inverted commas.

```{r 10 genes with highest mean}
ordered_expression <- gene_expression[order(-gene_expression$mean_expression),]
genes_highest_expression <- head(ordered_expression, 10)
```

(Note: The following sentences contain in line code)

The 10 genes with the highest mean expression are:

1. `r rownames(genes_highest_expression)[1]`

2. `r rownames(genes_highest_expression)[2]`

3. `r rownames(genes_highest_expression)[3]`

4. `r rownames(genes_highest_expression)[4]`

5. `r rownames(genes_highest_expression)[5]`

6. `r rownames(genes_highest_expression)[6]`

7. `r rownames(genes_highest_expression)[7]`

8. `r rownames(genes_highest_expression)[8]`

9. `r rownames(genes_highest_expression)[9]`

10. `r rownames(genes_highest_expression)[10]`


#### 1.4. Determine the number of genes with a mean < 10

In order to determine the number of genes with mean expression <10, a subset is created from the main data frame with 'mean_expression<10' as the criteria for subsetting.
The total number of genes in this subset can therefore be determined by retrieving the number of rows of this new subset data frame, and the value displayed in line with text using 'r (code)' with the specialised inverted commas.

```{r number of genes with low expression < 10}
low_exp_genes <- subset(gene_expression, mean_expression<10)
```

(Note: The following sentences contain in line code)

The number of genes with a mean expression <10 is `r nrow(low_exp_genes)`.

#### 1.5. Make a histogram plot of the mean values and include it into your report.

Histograms area useful way to visualize expression data. 

Before we construct the graph, it is beneficial to first remove genes from the analysis where zero expression was detected in these samples. 
In this way, we can better focus on the expression data from the genes that are expressed in these samples.

Next, the histogram is created, with the column containing the mean expression values is selected from the appropriate data frame using the $ code.
Additional bins are included to aid in visualising the data spread. Titles for the histogram and X axis, plus the colour, are also specified.

```{r histogram of expression values, fig.cap="*Figure 1. Mean global gene expression from three samples of interest, showing severe right skewing in the data, indicating a large proportion of genes with low expression levels.*"}
gene_exp_greaterthan0 <- subset(gene_expression, mean_expression>0)
hist(gene_exp_greaterthan0$mean_expression, breaks=100, main = "Mean Gene Expression", xlab = "Mean Gene Expression", col = "dodgerblue2")
```

From this histogram, we can see that the vast majority of genes have very low expression. This severe right skewing of the data makes it difficult to glean useful information from this graph, as the large number of lowly expression genes dominates.

Data with a distribution such as this is better visualized on a logarithmic scale.

As there would be genes in this data with mean expression levels below 1, we can add 1 to all mean expression values, to ensure that the log of the values stays positive.

Here is the resulting histogram.

```{r histogram of log mean values, fig.cap="*Figure 2. Mean global gene expression from three samples of interest, on a logarithmic scale.*"}
log_gene_exp_greaterthan0 <- log(gene_exp_greaterthan0+1)
hist(log_gene_exp_greaterthan0$mean_expression, breaks=50, main = "Mean Gene Expression (log)", xlab = "Log Mean Gene Expression", col = "dodgerblue2")
```
























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