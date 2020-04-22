---
title: "R code for bayesian analysis in Lightfoot et al 2019"
author: "Danny Buss"
date: "14/09/2019"
output: html_document
---
  
## R Markdown
This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

Load required packages:
```{r}
library(parallel)
library(coda)
library(MCMCglmm)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(dplyr)
library(reshape2)
library(lme4)
library(nlme)
```
Load dataframe:
```{r}
df <-read.csv("Data.csv")
```

Load functions:
```{r}
plot.estimates <- function(x) {
    if (class(x) != "summary.mcmc")
      x <- summary(x)
    n <- dim(x$statistics)[1]
    par(mar=c(2,13,4,1))
    plot(x$statistics[,1], n:1,
         yaxt = 'n', ylab='',
         xlim=range(x$quantiles)*1.2,
         pch=19,
         main='Posterior means and 95% credible intervals')
    grid()
    axis(2, at=n:1, rownames(x$statistics), las=2)
    arrows(x$quantiles[,1], n:1, x$quantiles[,5], n:1, code=0)
    abline(v=0, lty=2)
}

trace.plots<-function(x) {
  n <- dim(x)[2]
  par(mfrow=c(ceiling(n/2),2), mar=c(0,0.5,1,0.5))
  for (i in 1:n) {
    plot(as.numeric(x[,i]), t="l", main=colnames(x)[i], xaxt="n", yaxt="n")
  }
}

mu.link<-function(weight) post$a + post$b*weight

plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow=c(ceiling(n/2),2), mar=c(3,2,3,0))
  for (i in 1:n) {
    acf(x[,i], lag.max=100, main=colnames(x)[i])
    grid()
  }
}
```

#Bayesian Mixed-Effect Models
Table 1
```{r}

```
#Paired-site model
Table 1
```{r}

```
#TimePeriod model
Table 1
```{r}

```
#Site model
Table 1
```{r}

```
#Plots
Table 1
```{r}

```

#Variances - Levene's test

```{r}

```