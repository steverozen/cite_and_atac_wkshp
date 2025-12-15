# Preparation for R and Single Cell Multiomics Workshop

### Dec 17, 2025

## Install R, the Seurat R package and other R packages, and RStudio or other IDE

### IMPORTANT: This may take over an hour. You must do this before the workshop

Send email to sr110\@duke.edu if you have problems.

Please let me know if you find errors in these instructions.

#### 1 If not already installed, install R from https://cloud.r-project.org/

Select Download R for macOS or Windows (or Linux) depending on your computer's operating system.

#### 1.1 If R is installed, but version is \< 4.5 I suggest you upgrade now using the instructions above

To find your version of R type `R.version.string` at the R prompt.

#### 2 If not already installed, install RStudio from https://posit.co/download/rstudio-desktop/

Select choice "2: Install RStudio", assuming you already have R installed.

If you prefer you can use Positron: https://posit.co/products/ide/positron/ or VS Code instead

Make sure you can run RStudio. If when you run it it suggests an update then update it

#### 3 In RStudio (or other interactive development environment), in an R console, enter the following R commands

These install packages you will need and test that you can load and use the packages.

``` r
install.packages("BiocManager")
```

``` r
library(BiocManager)
BiocManager::version()
```

``` r
BiocManager::install("GenomicRanges")
```

``` r
library(GenomicRanges)
gr0 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               IRanges(1:10, width=10:1))
```

``` r
install.packages("Signac")
```

``` r
library(Signac)
x <- matrix(data = sample(c(0, 1), size = 25, replace = TRUE), ncol = 5)
Jaccard(x = x, y = x)
```

``` r
install.packages(c('Seurat'))
```

**Important** you need Seurat v5. If you already have Seurat installed, to check, do

``` r
packageVersion("Seurat")
```

``` r
library(Seurat)
pbmc_raw <- read.table(
  file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
  as.is = TRUE
)
```

``` r
install.packages(c('ggplot2', 'patchwork'))
```

``` r
library(patchwork)
library(ggplot2)
p1 <- ggplot(mtcars) + geom_point(aes(mpg, disp))
p2 <- ggplot(mtcars) + geom_boxplot(aes(gear, disp, group = gear))
p1 + p2
```

``` r
BiocManager::install("AnnotationHub")
```

``` r
library(AnnotationHub)
AnnotationHub()[1]
```

``` r
BiocManager::install("AnnotationDbi")
```
```r
BiocManager::install("org.Hs.eg.db")
```
```r
install.packages("igraph")
```
```r
BiocManager::install('enrichplot')
```

This next package `clusterProfiler` has a lot of dependencies.
If there issues installing dependencies you may 
have better success installing the dependencies
indiviually before trying to install
`clusterProfiler` again.

```r
BiocManager::install("clusterProfiler") 
```

###### If one or more of these failed you may need to update R or RStudio

###### 

###### 4.2 If you have other problems update R

-   See the instructions in point 1 above to install R.
-   In RStudio, you may have to point RStudio to the correct version by using the Tools drop down, then Global Options
-   In any case, you will need to restart RStudionand confirm that it using the expected new verions. Type `R.version.string` at the R prompt.
-   Reinstall all packages above (they may need to be upgraded).

### Download data needed for the workshop