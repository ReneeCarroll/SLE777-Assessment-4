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

Note that in the R.Markdown report, all scripts have "echo=FALSE" included, so the code chuncks do not appear in the final report.

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

To generate the required table, the 'knitr::kable()' code was used, with caption included.

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
The following uses in-line code: 

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

# Part 1
## Tree Circumference Data Analysis

#### 1.6. Import this csv file into an R object. 

The file “growth_data.csv” contains measurements for the circumferences of trees growing at two sites, a control site and a treatment site, which were planted 20 years ago.

To begin the analysis, the tree circumference data file must first be read into R. 

The following code will download the required data file, and will read it into R. 
As the file is a .csv file, and as the R code read.csv() automatically treats the first row as headings, we don't need to specify additional information. 

```{r Import csv file}
URL2 = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"
download.file(URL2, destfile="growth_data.csv")
growth_data <- read.csv("growth_data.csv") # read.csv automatically treats 
```

#### What are the column names?

(Note: The following sentences contain in line code)

The column names of this data set are:

1. `r colnames(growth_data)[1]`

2. `r colnames(growth_data)[2]`

3. `r colnames(growth_data)[3]`

4. `r colnames(growth_data)[4]`

5. `r colnames(growth_data)[5]`

6. `r colnames(growth_data)[6]`


#### 1.7. Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites.

To show these important data statistics, we can create a new data frame and specify the data that should populate it.
Here, the new data frame is called 'girth_stats_df'. The first column is given the name 'Location', and contains a character vector providing the names of the two sites.
Columns 2 to 5 are similarly provided first with a name, then the vector of data that is to populate those columns. 
The subset code used here requires us to specify the data frame, then the variable to subset by, which in this case is each of the locations in the Site column, followed by the column of data that's required to be copied into the new data frame, in this case, the various circumference measurements. The subset values are aggregated using the mean() or sd() function, and saved to the appropriate column. 

In this way, the data frame is built with the required information. 

For ease of data viewing, all values contained in columns 2 to 5 have been rounded to 2 decimal places.

```{r Create new data frame for the mean and sd}
girth_stats_df <- data.frame(
  Location = c("Northeast", "Southwest"),
  Mean_2005_cm = c(mean(subset(growth_data, Site == "northeast")$Circumf_2005_cm), mean(subset(growth_data, Site == "southwest")$Circumf_2005_cm)),
  SD_2005_cm = c(sd(subset(growth_data, Site == "northeast")$Circumf_2005_cm), sd(subset(growth_data, Site == "southwest")$Circumf_2005_cm)),
  Mean_2020_cm = c(mean(subset(growth_data, Site == "northeast")$Circumf_2020_cm), mean(subset(growth_data, Site == "southwest")$Circumf_2020_cm)),
  SD_2020_cm = c(sd(subset(growth_data, Site == "northeast")$Circumf_2020_cm), sd(subset(growth_data, Site == "southwest")$Circumf_2020_cm))
)
girth_stats_df[,2:5] <- round(girth_stats_df[,2:5], digits = 2)
```

The following table shows the mean and standard deviation of tree circumference at the start and end of the study at both sites

``` {r Print the new data frame of mean and sd}
knitr::kable(girth_stats_df,
             caption = "Table 3. Girth statistics at the start and end of study")
```

#### 1.8. Make a box plot of tree circumference at the start and end of the study at both sites.

To create the boxplot, we'll again use the subset code to select the data we want plotted in the boxplot, from the original growth data frame. 

```{r boxplot of tree circumference at the start and end of the study at both sites, fig.cap="*Figure 3. Boxplot of tree girth measurements (in cm) taken in 2005 and 2020 at the Norheast and Southwest study sites. The median value is shown by the black horizontal central line. The interquartile range, comprised of the upper quartile above the median, and the lower quartile below the median, is represented by the coloured box. The maximum and minimum values are shown by the ends of the whiskers.*"}
boxplot(
  subset(growth_data, Site == "northeast")$Circumf_2005_cm, 
  subset(growth_data, Site == "northeast")$Circumf_2020_cm, 
  subset(growth_data, Site == "southwest")$Circumf_2005_cm, 
  subset(growth_data, Site == "southwest")$Circumf_2020_cm, 
  names = c("Northeast 2005", "Northeast 2020", "Southwest 2005", "Southwest 2020"), 
  col = c("forestgreen"),
  ylab = "Circumference (cm)", 
  main = "Increase in tree girth between 2005 and 2020, at each study site")
```

#### 1.9. Calculate the mean growth over the last 10 years at each site.

To calculate and display the mean growth from 2010 to 2020, use the $ code to add a column to the original data frame to contain this information. 
To this column, assign the difference in circumference measurements from 2010 to 2020.
Next, create a new data frame, and save into this the site names, and the calculated mean of the difference in measurements from the main data frame, using the subset() code to select the appropriate data as per each location. Including the code 'check.names=FALSE' ensures that the column names appear as intended.

```{r mean growth over the last 10 years}
# Add a column to the existing data frame, to contain the calculation of 2020 measurements minus the 2010 measurements.
growth_data$growth_2010_2020 <- growth_data$Circumf_2020_cm-growth_data$Circumf_2010_cm 

# Create a new data frame to present the mean of the 10 year growth data, for each location.
mean_10yr_growth_df <- data.frame(
  Location = c("Northeast", "Southwest"),
  "Mean 10 year growth (cm)" = c(mean(subset(growth_data, Site == "northeast")$growth_2010_2020), mean(subset(growth_data, Site == "southwest")$growth_2010_2020)),
  check.names = FALSE 
)
```

Here is the mean growth, measured in cm, over the last 10 years at each of the sites.

```{r table showing mean growth over 10 years}
knitr::kable(mean_10yr_growth_df,
             caption = "Table 4. Mean growth over last 10 years")
```

#### 1.10. Use the t.test to estimate the p-value that the 10 year growth is different at the two sites.

Prior to conducting a t-test, it is important to fist check for equal variance in the data.
Equal variance can be calculated using leveneTest() code as part of the "car" package.

``` {r}
leveneTest(growth_2010_2020 ~ Site, data = growth_data)
```
As the Levene's Test returned a P-value of 0.006116, it indicates the variances in this data are not equal, and therefore the appropriate t-test to use is a Welch's t-test, as it adjusts for unequal variances.

The output for Welch's t-test is below:

```{r t.test)}
t.test(growth_2010_2020 ~ Site, data = growth_data, var.equal = FALSE) # var.equal = FALSE isn't required here, as it is the default for t.test() in R.
```
This t-test provided a p-value of 0.06229, indicating the difference in 10 year growth between the two sites is not statistically significant, however as this value is only slightly above the cutoff of 0.05, is is close to being significant.

The recommendation is to keep monitoring tree growth at the two sites, and repeat the analysis in the future.


# Part 2
## Examining biological sequence diversity between E.coli and A.aceti


Here in this report, we will perform a comparative analysis on the coding sequences of two bacteria, Escherichia coli (E.coli) and Acetobacter aceti (A.aceti).

The analysis utilizes two packages installed and loaded earlier; the Seqinr package, to process and analyse sequence data, and the R.utils package, which is a general utilities package enabling files to be unzipped. 


#### 2.1.1. Download the whole set of coding DNA sequences for E. coli and your organism of interest.

The coding DNA sequence files (cds files) can be downloaded from Ensembl's database, using the respective URLs in the code chunk below. 
Once downloaded and saved, the FASTA file can then be unzipped using the gunzip() code from the "R.utils" package.

Both FASTA files are read into R using the read.fasta() code as part of the "seqinr" package.

```{r Download E.coli coding DNA sequences}
# Download the E.coli cds, and unzip the file.
URL3 = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
#download.file(URL3, destfile="ecoli_cds.fa.gz")
#gunzip("ecoli_cds.fa.gz")

# Download the cds for my organism of interest, Acetobacter aceti GCA_002005445, and unzip the file.
URL4 = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_80_collection/acetobacter_aceti_gca_002005445/cds/Acetobacter_aceti_gca_002005445.ASM200544v1.cds.all.fa.gz"
#download.file(URL4, destfile="aaceti_cds.fa.gz")
#gunzip("aaceti_cds.fa.gz")

# Read in the two sequence files.
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
aaceti_cds <- seqinr::read.fasta("aaceti_cds.fa")
```

#### 2.1.2. How many coding sequences are present in these organisms? Present this in the form of a table. 

Within the FASTA files, each line contains the coding DNA sequence for a different gene. Therefore, to number of coding seuqences present in each organism is given by the number of rows in the FASTA file. In this context, the length() code provides us with this answer. By creating a new data frame with this data, we can then easily present it in the form of a table.

```{r Table of the number of coding sequences present in these organisms}
cds_df <- data.frame(
  Organism = c("E.coli", "A.aceti"),
  "Number of Coding Sequences" = c(length(ecoli_cds), length(aaceti_cds)),
  check.names = FALSE 
  )
knitr::kable(cds_df,
             caption = "Table 5. Number of coding sequences present in E.coli and A.aceti")
```

#### 2.1.3. Describe any differences between the two organisms. 

(Note: The following sentences contain in line code)

The organism E.coli has `r (length(ecoli_cds)-length(aaceti_cds))` additional genes in its coding DNA, which is an increase of `r round((length(ecoli_cds)-length(aaceti_cds))/length(aaceti_cds)*100, 1)`%. 

This data suggests that E.coli is a more complex organism than A.aceti, as a greater number of genes would result in a greater number of proteins, and hence additional or more complex biological functionality.

#### 2.2.1 How much coding DNA is there in total for these two organisms? Present this in the form of a table.

In order to determine the total amount of coding DNA for each organism, we can perform the 'summary()' code to show the length of each coding sequence (presented in the first column). Then, we can take the values from this first column and save this information as a new number vector. Taking the sum of these new vectors gives us the total amount of coding DNA for each organism, which can then be saved as a data frame and presented in the form of a table, as shown below. 

```{r Table of the total amount of coding DNA for the two organisms}
length_ecoli_genes <- as.numeric(summary(ecoli_cds)[,1])
length_aaceti_genes <- as.numeric(summary(aaceti_cds)[,1])
cds_total_df <- data.frame(
  Organism = c("E.coli", "A.aceti"), 
  "Total Amount of Coding Sequence" = c(sum(length_ecoli_genes), sum(length_aaceti_genes)),
  check.names = FALSE 
) # Including the check.names = FALSE code so that R doesn't add '.' in between the words of the Column names.
knitr::kable(cds_total_df,
             caption = "Table 6. Total quantity of coding DNA in E.coli and A.aceti")
```

#### 2.2.2. Describe any differences between the two organisms. 

(Note: The following sentences contain in line code)

In line with the previous observation that E.coli contains `r round((length(ecoli_cds)-length(aaceti_cds))/length(aaceti_cds)*100, 1)`% more genes than A.aceti, we can further determine that E.coli has a proportinal increase in the total amount of coding sequence as well, namely, that E.coli contains `r  round((sum(length_ecoli_genes)-sum(length_aaceti_genes))/sum(length_aaceti_genes)*100, 1)`% more coding DNA in total. This further supports the hypothesis that E.coli is a more complex organism than A.aceti.

#### 2.3.1. Calculate the length of all coding sequences in these two organisms. Make a boxplot of coding sequence length in these organisms. 

Using the vectors of cDNA lengths just created, the boxplot can be easily generated. 

```{r Boxplot of coding sequence length in these organisms, fig.cap="*Figure 4. Distributions of coding sequence lengths between two organisms, E.coli (red) and A.aceti (blue). The median values are shown by the black horizontal central lines. The interquartile range, comprised of the upper quartile above the median, and the lower quartile below the median, are represented by the coloured boxes. The maximum and minimum values are shown by the ends of the whiskers, with outlying values shown by points above the whisker border line.*"}
boxplot(length_ecoli_genes, length_aaceti_genes, names = c("E.coli", "A.aceti"), main = "Coding DNA Sequence Lengths in E.coli and A.aceti", ylab="Length of cDNA (bp)", col=c("red", "blue"))
```

#### 2.3.2 What is the mean and median coding sequence length of these two organisms? Describe any differences between the two organisms. 

(Note: The following sentences contain in line code)

The average length of a coding sequence in E.coli is `r round(mean(length_ecoli_genes), 0)` bp long, and in A.aceti is `r round(mean(length_aaceti_genes), 0)` bp long. The median coding sequence length in E.coli is `r median(length_ecoli_genes)`, and in A.aceti is `r median(length_aaceti_genes)`. 

Together with the boxplot shows in Figure 4 above, this data reveals that the transcriptomes of E.coli and A.aceti are comparable in terms of gene length. In both organisms, the vast majority of coding sequences are around the 800-900 bp length, with both organisms having a smaller number of very large genes, in the approximate range of 2500-8000 bp long, as shown by the data points above the plots' upper whisker lines. The presence of these large genes are the reason the calculated average length is larger than the calculated median length, in both organisms. 

#### 2.4.1. Calculate the frequency of DNA bases in the total coding sequences for both organisms. Create a bar plot for nucleotide frequency. 

Using the 'unlist()' command, the cDNA sequences from each organism's coding DNA sequence (cds) data frame can be saved as a separate vector. The 'count()' command can then be used on this vector, with the instruction to count the number of nucleotides present in the vector in combinations of 1's (ie. as single nucleotides, a, c, g, and t), saving this result as another object.

The proportion of each nucleotide in each organism's total coding sequence can then be calculated, with the results graphed in a grouped bar graph, using the packages "ggplot2" and "ggpubr" (see Figure 5).

```{r Nucleotide composition and frequencies in E.coli and A.aceti}

# Determining DNA composition for each organism
dna_ecoli <- unlist(ecoli_cds) 
dna_composition_ecoli <- count(dna_ecoli,1)

dna_aaceti <- unlist(aaceti_cds) 
dna_composition_aaceti <- count(dna_aaceti,1) 

# Determining DNA proportion for each organism
dna_proportion_ecoli <- dna_composition_ecoli/sum(dna_composition_ecoli)
dna_proportion_ecoli <- round(dna_proportion_ecoli, digits = 3) # digits=3 means that the result is rounded to 3 decimal places.

dna_proportion_aaceti <- dna_composition_aaceti/sum(dna_composition_aaceti)
dna_proportion_aaceti <- round(dna_proportion_aaceti, digits = 3)
```

``` {r Creating the bar chart for nucleotide frequencies using ggplot2, fig.cap="*Figure 5. Grouped bar chat showing the frequencies of nucleotide composition in the two organisms, E.coli and A.aceti.*"}
# Create data frames from the nucleotide proportions, and add the column for Organism
dna_proportion_ecoli_df <- as.data.frame(dna_proportion_ecoli)
dna_proportion_ecoli_df$Organism <- c("E.coli")

dna_proportion_aaceti_df <- as.data.frame(dna_proportion_aaceti)
dna_proportion_aaceti_df$Organism <- c("A.aceti")

# Merge the two data frames together, to create one data frame for the bar plot
proportion_df <- rbind(dna_proportion_ecoli_df, dna_proportion_aaceti_df)
names(proportion_df) <- c("Nucleotide", "Proportion", "Organism") 

# Define the variables for the grouped bar chart
Organism <- proportion_df$Organism
Nucleotide <- proportion_df$Nucleotide
Proportion <- proportion_df$Proportion

# Create the grouped bar chart. 
# The code specifies the data frame to use, and then which variables to plot where. It defines the position is side-by-side with the code 'position=position_dodge', and ensures there is a space between the columns with the code '(0.8)'. A blue colour palette was chosen for the graph with the code scale_fill_brewer(palette = "Blues".
ggbarplot(data=proportion_df, x="Nucleotide", y="Proportion", fill="Organism", position=position_dodge (0.8))+
  scale_fill_brewer(palette = "Blues") 
```

#### 2.4.2. Describe any differences between the two organisms.

While nucleotide composition is comparable between E.coli and A.aceti, it is observed that the coding genome of the organism of interest, A.aceti, contains an decreased proportion of adenine (a) and thymine (t), and an increased proportion of cytosine (c) and guanine (g), as compared to the coding genome of E.coli. 

#### 2.4.3. Calculate the frequency of amino acids in the total protein sequences for both organisms. Create a bar plot for amino acid frequency. 

In order to analyse the protein sequences in these organisms, the coding sequences must first be translated to their respective protein sequences.
Here, the lapply function performs a loop, applying the 'translate' code to all elements of each organisms's cDNA data, and the output is saved as new object. 

Next, the 'unlist()' command is again used to extract the amino acid sequences from the newly created data objects, and the total protein sequence is saving as a new single vector. 

We then use the count() command again, but we first must define the list of possible amino acids which the count() command is to reference. To do this, we can determine the unique elements present in one of the organism's proteome vectors, remove all asterisks that represent the stop codons, and save it as a new vector 'aa'. Then, when the count() command is used here, we define the alphabet to be the newly created 'aa' vector, and the wordsize to be 1, such that the number of each of the single amino acids is calculated from the organisms entire proteome.

The proportion of each amino acid in the total proteome can then be calculated for each organism, with the results graphed in a grouped bar graph, using the packages "ggplot2" and "ggpubr" (see Figure 6).

```{r Determining the frequency of amino acids in the total protein sequences of each organism}
# Translate the cDNA for each organism
ecoli_prot <- lapply(ecoli_cds, translate)
aaceti_prot <- lapply(aaceti_cds, translate)

# Unlist the protein sequences and saving then as a single vector
ecoli_prot_total <- unlist(ecoli_prot)
aaceti_prot_total <- unlist(aaceti_prot) 

# Creates the list of possible amino acids to form our alphabet to use in the 'count' code.
aa <- unique(ecoli_prot_total) 

# Removes the stop codon asterisks.
aa <- aa[aa != "*"] 

# Count the number of each amino acid in each organisms proteome vector, in combinations of 1's (ie. single amino acids), using the aa alphabet previously determined, and save it as a new object. 
aa_composition_ecoli <- count(ecoli_prot_total, wordsize=1, alphabet=aa) 
aa_composition_aaceti <- count(aaceti_prot_total, wordsize=1, alphabet=aa)

# Calculate the proportion of each amino acid in each organisms' proteome
aa_proportion_ecoli <- aa_composition_ecoli/sum(aa_composition_ecoli)
aa_proportion_ecoli <- round(aa_proportion_ecoli, digits = 3)
aa_proportion_aaceti <- aa_composition_aaceti/sum(aa_composition_aaceti)
aa_proportion_aaceti <- round(aa_proportion_aaceti, digits = 3)
```

```{r Creating the bar chart for amino acid frequencies using ggplot2, fig.cap="*Figure 5. Grouped bar chat showing the frequencies of amino acid composition in the two organisms, E.coli and A.aceti.*"}
#Create data frames from the aa proportions, and add the column for Organism
aa_proportion_ecoli_df <- as.data.frame(aa_proportion_ecoli)
aa_proportion_ecoli_df$Organism <- c("E.coli")

aa_proportion_aaceti_df <- as.data.frame(aa_proportion_aaceti)
aa_proportion_aaceti_df$Organism <- c("A.aceti")

#Merge the two data frames together, to create one data frame for the bar plot
aa_proportion_df <- rbind(aa_proportion_ecoli_df, aa_proportion_aaceti_df)
names(aa_proportion_df) <- c("Amino_Acids", "Proportion_aa", "Organism_aa") 

#Defining the variables for the grouped bar chart
Organism_aa <- aa_proportion_df$Organism_aa
Amino_Acids <- aa_proportion_df$Amino_Acids
Proportion_aa <- aa_proportion_df$Proportion_aa

#Create the grouped bar chart
#The code specifies the data frame to use, and then which variables to plot where. It defines the position is side-by-side with the code 'position=position_dodge', and ensures there is a space between the columns with the code '(0.8)'. A green colour palette was chosen for the graph with the code scale_fill_brewer(palette = "Greens").
ggbarplot(data=aa_proportion_df, x="Amino_Acids", y="Proportion_aa", fill="Organism_aa", position=position_dodge (0.8))+
  scale_fill_brewer(palette = "Greens") 
```

#### 2.5.1. Create a codon usage table and quantify the codon usage bias among all coding sequences. 

Codon usage can be calculated by the 'uco()' command as part of the "seqinr" package. In order to calculate codon usage for all coding genes in each organism, the command is used on the vectors of 'unlisted' cDNA sequences, and the output is saved as a new data frame.

The code here includes instructions to calculate the RSCU value, which stands for Relative Synonymous Codon Usage. 
This is a measure that tells us whether there is any preferential usage of a particular codon for a particular amino acid, as compared to the relative usage of other codons for that same amino acid. 

The two codon usage data frames are then merged into a new data frame, for easy vewing in Table 7 below.

```{r Creating the codon usage table}
# Create the codon usage data frame for E.coli, showing RSCU values. 
codons_ecoli_df <- uco(dna_ecoli, index = "rscu", as.data.frame=TRUE)
codons_ecoli_table <- codons_ecoli_df
colnames(codons_ecoli_table) <- c("AA", "codon", "E.coli_eff", "E.coli_freq", "E.coli_RSCU")

# Create the codon usage data frame for A.aceti, showing RSCU values
codons_aaceti_df <- uco(dna_aaceti, index = "rscu", as.data.frame=TRUE)
codons_aaceti_table <- codons_aaceti_df
colnames(codons_aaceti_table) <- c("AA", "codon", "A.aceti_eff", "A.aceti_freq", "A.aceti_RSCU")

# Creating a new data frame from 'codons_ecoli_df' data frame, and adding the A.aceti columns from the 'codons_aaceti_df' data frame.
codons_table <- codons_ecoli_table
codons_table$A.aceti_eff <- codons_aaceti_table$A.aceti_eff
codons_table$A.aceti_freq <- codons_aaceti_table$A.aceti_freq
codons_table$A.aceti_RSCU <- codons_aaceti_table$A.aceti_RSCU
knitr::kable(codons_table,
             row.names = FALSE,
             caption = "Table 7. Codon usage statistics of E.coli and A.aceti. AA: Amino Acid three leter code, codon: encoding codon, eff: total codon counts, freq: relative frequencies, RSCU: Relative Synonymous Codon Usage.")
```

#### 2.5.2. Provide charts to support your observations.

The RSCU value is useful to focus on when comparing codon preference between organisms. 
An RSCU value of 1 means the codon is used as expected, values above 1 indicate the codon is being used more frequently than expected, and values below 1 indicate the codon is being used less frequently. 

The following code sets up a new data frame which is then the source for creating a grouped bar plot of Relative Synonymous Codon Usage in the organisms E.coli and A.aceti. By subtracting 1 from all RSCU values, the relative usage patters can be more easily discernible, as they fall either above or below a zero line, as seen in Figure 6(a) and Figure 6(b) below. 

```{r Creating the data frame which will be used to create the grouped bar chart for RSCU values}
# Using the original codon data frames, add a column for organism name
codons_ecoli_df$Organism <- c("E.coli")
codons_aaceti_df$Organism <- c("A.aceti")

# Merge these data frames into one
column_names <- colnames(codons_ecoli_df) # creates the vector of column names
codons_df <- rbind(codons_ecoli_df, codons_aaceti_df) # combines the rows of data
names(codons_df) <- column_names # assigns the column names of the new data frame to be the names from the set vector

# Create a new column that calculates the RSCU value minus 1, such that values are either above or below a zero line.
codons_df$"RSCU-1" <- codons_df$RSCU-1

# Create a new column that contains the amino acid name plus it's codon
codons_df$"AA_and_codon" <- paste(codons_df$AA, ":", codons_df$codon)

# Put the data in alphabetical order, according to the amino acid 3 letter code
ordered_codons_df <- codons_df[order(codons_df$AA), ]
```

```{r Grouped barplot (a) showing RSCU-1 values, fig.cap="*Figure 6(a). Grouped bar chat showing Relative Synonymous Codon Usage in E.coli and A.aceti, for the first subset of amino acids. RSCU values minus 1 are shown on the y axis, and on the x axis is each of the amino acids and their various respective codons.*"}
# Subsetting the data frame using [] after the data frame name, to select rows 1 to 70 to graph.
# In the theme() component, I've set the column text size to 10, have adjusted the angle to 90 degrees, have justified the labels to be straight against the x-axis, and have shifted them left such that they line up under the tick lines. 
ggbarplot(data=ordered_codons_df[1:70, ], x="AA_and_codon", y="RSCU-1", fill="Organism", position=position_dodge (0.8)) +
  scale_fill_brewer(palette = "Oranges") +
  ggtitle("Relative Synonymous Codon Usage between E.coli and A.aceti") +
  xlab("Amino Acid : Codon") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust=1, vjust=0.5))
```

```{r Grouped barplot (b) showing RSCU-1 values, fig.cap="*Figure 6 (b). Grouped bar chat showing Relative Synonymous Codon Usage in E.coli and A.aceti, for the second subset of amino acids. RSCU values minus 1 are shown on the y axis, and on the x axis is each of the amino acids and their various respective codons. Note that 'Stp' refers to stop codons, and not a specific amino acid.*"}
# Subsetting the data frame using [] after the data frame name, to select rows 71 to 128 to graph.
# In the theme() component, I've set the column text size to 10, have adjusted the angle to 90 degrees, have justified the labels to be straight against the x-axis, and have shifted them left such that they line up under the tick lines. 
ggbarplot(data=ordered_codons_df[71:128, ], x="AA_and_codon", y="RSCU-1", fill="Organism", position=position_dodge (0.8)) +
  scale_fill_brewer(palette = "Oranges") +
  ggtitle("Relative Synonymous Codon Usage between E.coli and A.aceti") +
  xlab("Amino Acid : Codon") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust=1, vjust=0.5))
```

#### 2.5.3. Describe any differences between the two organisms with respect to their codon usage bias. 

Figure 6 (a) and (b) reveals that while there are many codons being used to similar frequencies in E.coli and A.aceti, such as the Leucine codon 'ctg', there are also some interesting differences, such as the preferred codon usage for Serine, which in A.aceti is 'tcc' and in E.coli is 'agc'. Also interesting is the preferred different in Stop codon usgae, as denoted by 'Stp' in Figure 6 (b), with A.aceti preferencing the codon 'tga' and E.coli preferencing the codon 'taa'.

Such differences in codon preference provide interesting insights into the evolutionary histories of these two bacteria.

#### 2.6.1. In the organism of interest, identify 10 protein sequence k-mers of length 3-5 which are the most over- and under-represented k-mers in your organism of interest. Are these k-mers also over- and under-represented in E. coli to a similar extent? 

K-mers are short substrings of biological sequence of a fixed length, denoted by 'k'. 
Here, we'll investigate the presence of protein sequence k-mers in the total protein sequences of E.coli and A.aceti, of length 3, 4 and 5 amino acids.

To achieve this, we'll again use the 'count()' command with the 'aa' alphabet, but this time we'll instruct the word size to be either 3, 4, or 5, to return the counts for k-mers of 3, 4 and 5 amino acids in length, respectively. Adding the code 'freq=TRUE' ensures the frequencies are also calculated.

Once these calculations have been saved to a new data frame, the data can then be ordered by frequency value, and top and bottom 10 rows saved to new data frames, using the head() and tail() commands, thus subsetting the most- and least- frequently used k-mers of each size. 

These calculations can be seen in in Tables 8 to 13 below.

```{r Extracting K-mers frequency for sizes 3, 4 and 5, and saving as separate data frames}
# Determine all K-mers of size 3, and their frequency of use
aaceti_kmer3 <- data.frame(count(aaceti_prot_total, wordsize=3, alphabet=aa, freq=TRUE))
colnames(aaceti_kmer3) <- c("Kmer", "Frequency")

# Order this data by frequency, highest to lowest, then extract the top 10 and bottom 10 rows, to save in new data frames
aaceti_kmer3 <- aaceti_kmer3[order(-aaceti_kmer3$Frequency), ] # the minus sign before the data frame name makes the ordering highest to lowest.
aaceti_kmer3_top10 <- head(aaceti_kmer3, 10)
aaceti_kmer3_bottom10 <- tail(aaceti_kmer3, 10)

# Repeat this process for all K-mers of size 4 amino acids
aaceti_kmer4 <- data.frame(count(aaceti_prot_total, wordsize=4, alphabet=aa, freq=TRUE))
colnames(aaceti_kmer4) <- c("Kmer", "Frequency")

# Order this data by frequency, highest to lowest, then extract the top 10 and bottom 10 rows, to save in new data frames
aaceti_kmer4 <- aaceti_kmer4[order(-aaceti_kmer4$Frequency), ] # the minus sign before the data frame name makes the ordering highest to lowest.
aaceti_kmer4_top10 <- head(aaceti_kmer4, 10)
aaceti_kmer4_bottom10 <- tail(aaceti_kmer4, 10)

# Repeat this process for all K-mers of size 5 amino acids
aaceti_kmer5 <- data.frame(count(aaceti_prot_total, wordsize=5, alphabet=aa, freq=TRUE))
colnames(aaceti_kmer5) <- c("Kmer", "Frequency")

# Order this data by frequency, highest to lowest, then extract the top 10 and bottom 10 rows, to save in new data frames
aaceti_kmer5 <- aaceti_kmer5[order(-aaceti_kmer5$Frequency), ] # the minus sign before the data frame name makes the ordering highest to lowest.
aaceti_kmer5_top10 <- head(aaceti_kmer5, 10)
aaceti_kmer5_bottom10 <- tail(aaceti_kmer5, 10)

# Create all tables
knitr::kable(aaceti_kmer3_top10,
             row.names = FALSE,
             caption = "Table 8. Ten most frequently used amino acid 3-mers in A.aceti.")

knitr::kable(aaceti_kmer3_bottom10,
             row.names = FALSE,
             caption = "Table 9. Ten least frequently used amino acid 3-mers in A.aceti.")

knitr::kable(aaceti_kmer4_top10,
             row.names = FALSE,
             caption = "Table 10. Ten most frequently used amino acid 4-mers in A.aceti.")

knitr::kable(aaceti_kmer4_bottom10,
             row.names = FALSE,
             caption = "Table 11. Ten least frequently used amino acid 4-mers in A.aceti.")

knitr::kable(aaceti_kmer5_top10,
             row.names = FALSE,
             caption = "Table 12. Ten most frequently used amino acid 5-mers in A.aceti.")

knitr::kable(aaceti_kmer5_bottom10,
             row.names = FALSE,
             caption = "Table 13. Ten least frequently used amino acid 5-mers in A.aceti.")
```

This analysis reveals that there are many 4-mers and 5-mers which are not present at all in the proteome of A.aceti, as revealed by frequencies of zero in Tables 11 and 13 above. 
In addition, the frequency values returned for the top most used 3-mers are all larger than the frequency values for the top most used 4-mers and 5-mers. This is as expected, given that increased length of k-mer also brings increased specialty, and and hence is present at a reduced frequency within the same data set.

For these reasons, further comparisons between A.aceti and E.coli will focus just on 3-mers.

In order to ascertain in E.coli the frequency of the same 3-mers identified to be the most and least used in A.aceti, as shown in Tables 8 and 9 respectively, we first extract the list of these 20 3-mers, saving as separate vectors.

After the counts of all 3-mers has been performed on the E.coli total protein sequence, using the 'count()' command as with A.aceti, we can then subset the E.coli 3-mer data frame to select just for the rows containing the 10 most- and 10 least- frequently used 3-mers from A.aceti, saving the data as separate data frames.

```{r Extracting lists of the 10 most over- and under- represented 3 K-mers in A.aceti}
Top10_kmer3 <- as.vector(aaceti_kmer3_top10$Kmer)
Bottom10_kmer3 <- as.vector(aaceti_kmer3_bottom10$Kmer)
```

```{r Extracting frequency data of 3-mers in E.coli proteome}
ecoli_kmer3 <- data.frame(count(ecoli_prot_total, wordsize=3, alphabet=aa, freq=TRUE))
colnames(ecoli_kmer3) <- c("Kmer", "Frequency")
```

```{r Subsetting the E.coli lists of 3-mers}
# Subset the E.coli list, using the top 10 A.aceti list, and save as a new df
# This code tells R to save to a new data frame all columns of data from rows in the 'ecoli_kmer3' data frame, where the values in the 'Kmer' column match a value in the 'Top10_kmer3' vector. 
ecoli_kmer3_topsubset <- ecoli_kmer3[ecoli_kmer3$Kmer %in% Top10_kmer3, ] 
ecoli_kmer3_topsubset$Organism <- c("E.coli") # Add a column for Organism, and note it as E.coli

# Subset the E.coli list, using the bottom 10 A.aceti list, and save as a new df
ecoli_kmer3_bottomsubset <- ecoli_kmer3[ecoli_kmer3$Kmer %in% Bottom10_kmer3, ]
ecoli_kmer3_bottomsubset$Organism <- c("E.coli") # Add a column for Organism, and note it as E.coli
```

#### 2.6.2. Provide plots to support your observations of k-mer usage between E.coli and your organism of interest. 

Now that the 3-mer frequencies have been determined in E.coli for the same 3-mers identified in A.aceti, we can present this data as grouped bar plots.

```{r Preparing the data frame}
#Prepare the A.aceti K-mer 3 data frames, by adding a column for organism
aaceti_kmer3_bottom10$Organism <- c("A.aceti")
aaceti_kmer3_top10$Organism <- c("A.aceti")

# Merge the two organisms' data frames that contain the frequency data for the top 10 3-K-mers from A.aceti
column_names_kmer3 <- colnames(aaceti_kmer3_bottom10) # Sets what the column names will be
Top_kmer3_df <- rbind(aaceti_kmer3_top10, ecoli_kmer3_topsubset)
names(Top_kmer3_df) <- column_names_kmer3

# Merge the two organisms' data frames that contain the frequency data for the bottom 10 3-K-mers from A.aceti
Bottom_kmer3_df <- rbind(aaceti_kmer3_bottom10, ecoli_kmer3_bottomsubset)
names(Bottom_kmer3_df) <- column_names_kmer3
```

```{r Boxplot comparison of A.aceti top 10 k-mers, fig.cap="*Figure 7. Grouped bar chat showing frequency of top-10 3-mers in A.aceti, and their relative occurence in the E.coli proteome.*"}
# This line keeps the levels in the same order as they appear in the table, and not graphed in alphabetical order, as is default for ggplot. As there are duplicates in the Kmer column, the 'levels = ' code specifies for only unique factors.
Top_kmer3_df$Kmer <- factor(Top_kmer3_df$Kmer, levels = unique(Top_kmer3_df$Kmer))
ggbarplot(data=Top_kmer3_df, x="Kmer", y="Frequency", fill="Organism", position=position_dodge (0.8))+
  ggtitle("Ten most over-represented amino acid 3-mers in A.aceti, \n with comparison to frequency of presentation in E.coli") +
  scale_fill_brewer(palette = "Greens") 
```

```{r Boxplot comparison of A.aceti bottom 10 k-mers, fig.cap="*Figure 8. Grouped bar chat showing frequency of bottom-10 3-mers in A.aceti, and their relative occurence in the E.coli proteome.*"}
# This line keeps the levels in the same order as they appear in the table, and not graphed in alphabetical order, as is default for ggplot. As there are duplicates in the Kmer column, the 'levels = ' code specifies for only unique factors.
Bottom_kmer3_df$Kmer <- factor(Bottom_kmer3_df$Kmer, levels = unique(Bottom_kmer3_df$Kmer))
ggbarplot(data=Bottom_kmer3_df, x="Kmer", y="Frequency", fill="Organism", position=position_dodge (0.8))+
  ggtitle("Ten most under-represented amino acid 3-mers in A.aceti, \n with comparison to frequency of presentation in E.coli") +
  scale_fill_brewer(palette = "Reds") 
```

#### 2.6.3. Why do you think these sequences are present at different levels in the genomes of these organisms?

Figures 7 and 8 show the ten most over- and under-represented 3-mers in A.aceti, respectively, with their corresponding frequencies observed in E.coli as comparison. These bar plots clearly show differences in observed frequencies between the two organisms, with E.coli having lower representation of 7/10 of A.aceti's most over-represented 3-mers, and higher representation of 9/10 of A.aceti's most under-represented 3-mers.

However, in order to clearly visualize the relative difference of these differences in frequencies, it's helpful to see the data on the same bar plot, as shown below in Figure 9.

```{r Boxplot comparison of A.aceti bottom and top 10 k-mers, fig.cap="*Figure 9. Grouped bar chat showing frequency of bottom-10 3-mers and top-10 kmers in A.aceti, and their relative occurence in the E.coli proteome.*"}
# Combine the two data frames, and order first by organism, then by frequency
combined_kmer3_df <- rbind(Top_kmer3_df, Bottom_kmer3_df)
combined_kmer3_df <- combined_kmer3_df[order(combined_kmer3_df$Organism, combined_kmer3_df$Frequency), ]

# This first line again keeps the order as per the data frame, and not alphabetical.
combined_kmer3_df$Kmer <- factor(combined_kmer3_df$Kmer, levels = unique(combined_kmer3_df$Kmer))
ggbarplot(data=combined_kmer3_df, x="Kmer", y="Frequency", fill="Organism", position=position_dodge (0.8))+
  ggtitle("Ten most under- and over-represented amino acid 3-mers in A.aceti, \n with comparison to frequency of presentation in E.coli") +
  scale_fill_brewer(palette = "Purples")+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust=1, vjust=0.5))
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