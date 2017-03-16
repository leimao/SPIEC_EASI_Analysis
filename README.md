# SPIEC_EASI_Analysis

Author: Lei Mao

Date: 3/16/2017

## Introduction

This is a R and Python script to analyze the association between different microbiome species using SPIEC-EASI.

Sparse InversE Covariance estimation for Ecological Association and Statistical Inference (SPIEC-EASI) is a statistical inference method developed by Richard Bonneau at NYU school of medicine. This package will be useful to anybody who wants to infer graphical models for all sorts of compositional data, though primarily intended for microbiome relative abundance data (generated from 16S amplicon sequence data). It also includes a generator for [overdispersed, zero inflated] multivariate, correlated count data. Please see the paper published in [PLoS Comp Bio](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). It has been employed to analyze the relationships between different microbiome taxonomies. One of such analysis was published in [Science](http://science.sciencemag.org/content/early/2016/04/13/science.aaf3229).

Here I am trying to reproduce such analysis on several datasets including American Gut Project dataset and Human Microbiome Project dataset.

## Run the Analysis

AG: American Gut Project

HMP: Human Microbiome Project

The input of the analysis was ".biom" files.

To Perform the analysis:

Run "spiec-easi_analysis.R" and then run "data_postprocessing.py".

The association between OTUs and number of associations were stored in "network_Order.csv" and "num_edge_Order.csv", respectively.
