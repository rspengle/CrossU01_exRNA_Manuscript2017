# CrossU01_exRNA_Manuscript2017

The code and data files included in this repository are intended for use in reproducing the analysis and figures for our manuscript: "Accuracy, reproducibility and consequences of bias in small RNA-seq for quantitative microRNA profiling: a multiple protocol study across multiple laboratories". The code presented is not intended for general use, and works only with the provided data and metadata files. 

## Getting Started

The code was written and tested on a Windows 10 machine using RStudio version 1.0.143 and R version 3.4.0. Instructions are provided assuming one has already installed R and RStudio (https://www.rstudio.com/products/rstudio/download/), and is running the code from the RStudio Desktop environment. Follow the installation instructions on the RStudio website. Next, either download the CrossU01_exRNA_Manuscript2017 respository from github, or clone the repository using git, if installed. If you are using a Windows operating system, we suggest cloning the repository somewhere high in your directory tree (e.g. C:/CrossU01_exRNA_Manuscript2017). Some of the path names are quite long, and Windows may complain or fail to find some of the data files if they are placed in a subdirectory. Once downloaded, open the 201708_CROSSU01_NBT_RESUBMISSION_PROJECT_V3_USED.Rproj file using RStudio. This project uses the R package "packrat" as a package manager, and the packrat sources are included here. Loading the R project should initiate the packrat repository and install the packages in the project folder. See the Prerequesits section for instructions on installing packrat. The R script files "201708_CrossU01_Rebuttal_Code.R" and "editingAnalysis.R" files contain all the code. They can be run without packrat, but depend on several libraries being installed. See Prerequesites section for the list of libraries to install if packrat is not being used.

### Prerequisites

R -- tested using v3.4.0 

RStudio (recommended) -- tested using version 1.0.143

packrat (recommended)

git (recommended)

R packages (and their respective dependencies): 
limma, data.table, ggplot2, Biostrings, RColorBrewer, scales, tools, xlsx, stringr, edgeR, pheatmap, ggbeeswarm, Hmisc, cowplot, MASS, DESeq2, vegan, dendsort


### Installing

Install R (required) and RStudio Desktop (optional) from their respective sites. Follow the instructions provided.


Download the CrossU01_exRNA_Manuscript2017 repository from github to a top-level folder (e.g. C:/). If using linux/mac OS, the download location does not matter. 

or

Clone the repository using git.
```
cd /path/to/git/repo
git clone https://github.com/rspengle/CrossU01_exRNA_Manuscript2017.git
```

Decompress the downloaded file, if necessary.

Open Rstudio and install packrat, if using. Run the following lines in the R console

```
if (!require("devtools")) install.packages("devtools")
devtools::install_github("rstudio/packrat")
```

If using packrat, once installed, open the 201708_CROSSU01_NBT_RESUBMISSION_PROJECT_V3_USED.Rproj file in the CrossU01_exRNA_Manuscript2017 folder using the Rstudio File > Open Project...

If not using packrat, or if packrat installation of packages fails, the required packages (see prerequisites) can be installed using the R console: 

```
install.packages(c("limma", "data.table", "ggplot2" , "Biostrings", "RColorBrewer", "scales", "tools" , "xlsx" ,"stringr" , "edgeR","pheatmap", "ggbeeswarm" , "Hmisc" ,"cowplot" , "MASS", "DESeq2", "vegan" , "dendsort"))
```

## Running the scripts

The code to regenerate the analysis and output from the manuscript are included in two .R files: 201708_CrossU01_Rebuttal_Code.R and editingAnalysis.R. They should open automatically when the RProject is loaded (File > Open Project... > /path/to/CrossU01_exRNA_Manuscript2017/201708_CROSSU01_NBT_RESUBMISSION_PROJECT_V3_USED.Rproj), or can be opened manually in RStudio (File > Open File...). All the required data and metadata are included under the data/ subdirectory. Importing these files is done automatically in the execution of the scripts, and assumes that the code is being run from within the top-level of the project directory. To check the current working directory in RStudio run 

```
getwd()
[1] "C:/CrossU01_exRNA_Manuscript2017""
```

To change to the project directory use the command (substituting your path to the repository):

```
setwd("C:/CrossU01_exRNA_Manuscript2017")  
```

Running the scripts will import all the data files, perform all the wrangling/manipulations and calculations to generate the data for the manuscript. The R scripts are intended to be run line-by-line interactively, which will display any graphs in the RStudio viewer, and allow one to see summary statistics printed to the console. The scripts can also be run non-interactively. In either case, an "output" directory will be created with several subdirectories. The subdirectories will include all graphs and tables used in the manuscript. 

Note: running the scripts (particularly 201708_CrossU01_Rebuttal_Code.R) could take several minutes to run, depending on your machine. The workflow was also run on a machine with 32 Gb RAM. The RAM required is far less than 32 Gb; however, it is not recommended to run on a low-memory machine. 


## Authors (Code)

* **Ryan spengler, Ph.D.** [rspengle](https://github.com/rspengle)

## Authors (Study/Manuscript)
* **Maria D. Giraldez, Ryan M. Spengler, Alton Etheridge, Paula Maria Godoy, Andrea J. Barczak, Srimeenakshi Srinivasan, Peter L. De Hoff, Kahraman Tanriverdi, Amanda Courtright, Shulin Lu, Joseph Khoory, Renee Rubio, David Baxter, Tom A. P. Driedonks, Hank P. J. Buermans, Esther N. M. Nolte-‘t Hoen, Hui Jiang, Kai Wang, Ionita Ghiran, Yaoyu Wang, Kendall Van Keuren-Jensen, Jane E. Freedman, Prescott G. Woodruff, Louise C. Laurent, David J. Erle, David J. Galas, Muneesh Tewari** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The analysis presented here is the end-result of a project funded by the NIH Extracellular RNA Communication Common Fund. The authors of the study/manuscript generated the small RNA-seq libraries analyzed here. 
* M. Giraldez was the lead author of the study, coordinating all the sequencing efforts across 9 participating laboratories. 
* Corresponding authors were D. Galas and M. Tewari. 
* The code was generated as a member of the Tewari lab at the University of Michigan https://medicine.umich.edu/dept/dcmb/muneesh-tewari-phd

