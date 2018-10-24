# mirApp
Shiny app to plot SNPs modulating microRNA-gene interactions in tumors throughout the genome

### PREREQUISITES

* Must have R installed on your local machine. To install R, go to https://cran.r-project.org/mirrors.html and choose the mirror closest to your local. Follow the instructions for installation.

* Once installed, open R and install packages "shiny" and "ggplot2". To install packages, type in the R prompt:

 `>install.packages(c("shiny", "ggplot2"))`

  Allow R to install any package dependencies.

* Install RStudio to run the app. To install RStudio, go to https://www.rstudio.com/products/rstudio/download/ and follow the instructions for your operating system.

### GETTING STARTED

* Download mirApp.tar.gz and untar it into whichever local directory you wish. To untar, type `tar -xvzf mirApp.tar.gz`
in terminal or decompress using your software of choice.

* What lives here:
  * **app.R** Shiny app code including UX and server side functions
  * **data/** directory containing all data files for running Shiny app.  

* NOTE! Data files are not included in the repo due to their size. Data files can be downloaded from https://northwestern.box.com/s/co7kzy6qo7c5dlwbilapnocg34j7zjiv and should be placed in the data/ subdirectory.

### RUNNING APP

* Open RStudio and set your working directory to the mirApp path location. To change your working directory in R, type in the R prompt:
`>setwd("path/to/mirApp")`
 and replace `"path/to/mirApp/"` to your path.

* In RStudio, open **app.R**. In the top right corner of the code editor, press the "Run App" dropdown menu and run "Run in Window".

* When running the app, a window will appear. In the left panel, you must choose a cancer type in the dropdown menu, and a miRNA. miRNA names must follow TCGA convention. An explanation is provided in the "About" tab. Delete the instructions within
the window before typing in your miRNA name.

* The plots will appear in the "Plot" tab. Because the data is very large, the app may be slow and may take time to generate the figures. Please be patient!

* If the app returns an error after your inputs (cancer type and miRNA name), this means that the miRNA is unconsidered within the cancer type in the study.

### FULL RESULTS TABLES (without using Shiny)

* Full tables of the regQTL analysis results for all trios may be downloaded as tab-separated files from https://northwestern.box.com/s/qwxeww6rrwmdlp4cxkiyf3pvot4j8sqx
