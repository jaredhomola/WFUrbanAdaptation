WFUrbanAdaptation is a research compendium that brings together the data
and analyses associated with Homola et al. (in review). R scripts are
provided in analysis/, human-readable data in extData/, and R-readable
data in an .rda file in data/.

Note: All analyses originally conducted and package built using R
version 3.5.1 in Windows 10

### To get started using the package

    options(repos=structure(c(CRAN="http://cran.r-project.org")))
    install.packages("devtools")
    library(devtools)
    install_github("jaredhomola/WFUrbanAdaptation")

Scripts in /analysis
--------------------

1.  Build RDA file (buildRDA.R)
2.  Install required packages (installPackages.R)
3.  PCA of environmental variables (envPCA.R)
4.  Genetic diversity measures (genDiversity.R)
5.  Genetic diversity empirical cumulative density function
    (diversityECDF.R)
6.  Genetic differentiation (neiGst)
7.  Multivariate analyses including DAPCs and PCoAs
    (multivariateAnalyses.R)
8.  Outlier analysis using percentile-based approach across replicates
    (outlierAnalysis)
9.  Genotype-environment association analysis using logistic regression
    (logisticRegressionGEA.R)

### Contact

Jared J. Homola  
<jaredhomola20@gmail.com>  
www.jaredhomola.com

#### Copyright (c) 2018 Jared J. Homola

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
