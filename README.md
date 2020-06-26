# Fine-map QTL with NIL phenotypes
An R shiny web app was developed to visualize the results from the high-throughput drug-response assays

***Link to Shiny app: [here](https://katiesevans9.shinyapps.io/QTL_NIL/)***

## EXPLANATION OF FUNCTIONALITY
To begin analysis, the user can find all data controls in a panel on the left-hand side of the screen. A test dataset is provided (user should check “Use sample data” checkbox) or the user can upload a file from their local computer (see below for input file tips).  

<img src="https://github.com/katiesevans/finemap_NIL/blob/master/images/setup.png" width="250"/>

The user should select the control condition and the drug condition. In some cases, the user might want to further choose a specific assay to view if multiple options are available. The user also has the option to view one of many drug-response traits by selecting a trait from the drop-down menu. Finally, the user should choose a chromosome to view, generally the chromosome which contains the highlighted QTL for that particular assay.

Along the top of the main panel, the user can navigate several tabs including “Control”, “Condition”, “Regressed”, “NIL Genotypes”, and “Help!”. The “Control”, “Condition”, and “Regressed” tabs each show the NIL genotypes along the selected chromosome (left) and the NIL phenotypes for the selected trait (right) in the control condition (control), raw drug condition (condition), or regressed drug condition (regressed). 

<img src="https://github.com/katiesevans/finemap_NIL/blob/master/images/regressed.png" width="750"/>

Genotypes for strains that are not NILs generated from the N2 and CB4856 strains (such as NILs generated from different parents or CRISPR-generated deletion strains) are represented by a grey bar. The user can hover their mouse above the NIL genotype or phenotype plots to see more information or zoom in on a specific area of the plot. 

<img src="https://github.com/katiesevans/finemap_NIL/blob/master/images/plotly.png" width="750"/>

Below this plot is an interactive datatable containing the pairwise strain comparisons for this condition and trait. The “NIL Genotypes” tab shows a plot of the NIL genotypes across all chromosomes, not just the chromosome selected by the user (top) and an interactive datatable with the genotypes of each strain across all chromosomes (bottom). 

<img src="https://github.com/katiesevans/finemap_NIL/blob/master/images/genos.png" width="750"/>

The final tab, “Help!” provides the user with the instructions detailed here to help them use the Shiny App.

---

In addition to these basic controls, the user also has access to several advanced features. The user can choose a subset of strains to plot by checking the box labeled “Show a subset of strains?” and unchecking the boxes next to strains the user wishes to omit. Additionally, the user can plot the location of one or more QTL as a vertical line on the NIL genotype plot by checking the “Show QTL?” box. The user then chooses how many QTL to show and uses the appropriate slider input below to designate the genomic positions of each QTL. 

<img src="https://github.com/katiesevans/finemap_NIL/blob/master/images/advanced.png" width="250"/>

Finally, if the “Show genotype?” box is checked, the genotype of each strain at each QTL position will be shown on the phenotype plot as an orange “N” representing N2 and a blue “C” representing CB4856.

<img src="https://github.com/katiesevans/finemap_NIL/blob/master/images/qtlgeno.png" width="750"/>

---

## Input files
Input files should be in the R data file format (.Rda or .RData) or a CSV file and should be the pruned output from the [easysorter pipeline](github.com/andersenlab/easysorter). It is important that the “condition” column of the dataframe contains both a drug and a control and the “control” column contains either the control or “None” (for the control of the control). Control regression will be performed in the application.

**Example dataframe structure**

experiment | assay | plate | row | col | strain | condition | control | trait | phenotype
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
NIL | a | 1 | A | 1 | N2 | water | None | median.EXT | 456.5 | 
NIL | a | 2 | A | 1 | N2 | zinc | water | median.EXT | 200.5 | 
... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

