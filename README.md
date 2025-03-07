# Evolution_of_gene_expression_in_seasonal_envrionment
This repository contains R scripts for the paper "Evolution of gene expression in seasonal environments." The scripts process gene expression data, perform clustering and comparative analyses, detect rhythmic patterns, and examine evolutionary rates.

# Overview
## 1. Data processing
This script processes gene expression data through multiple steps to prepare it for comparative analysis. The steps include:
- Filtering low-expression genes
- Converting count data to GeTMM values
- Extracting one-to-one orthologous genes

## 2. Hierarchical Clustering
This script performs hierarchical clustering on both samples and genes based on gene expression and climatic data. Additionally, it evaluates the relationship between sample clustering and climatic conditions. The results are used to generate Figure 2 in the manuscript.

## 3. Principal Component Analysis (PCA)
This script performs Principal Component Analysis (PCA) to analyze gene expression patterns.

## 4. Rhythm Detection
This script detects rhythmic gene expression patterns using the RAIN algorithm (Thaben and Westermark, 2014) and calculates the period and peak in seasonal gene expression.

## 5. Species Comparison
This script performs a comparative analysis of gene expression patterns across species. It includes:
- Correlation analysis
- Calculation of mean absolute difference
- Calculation of seasonal peak difference
- Multiple comparisons of correlation (Nemenyi test)
- Gene Ontology enrichment analysis of genes with conserved/diverged expression pattern

## 6. Evolutionary Rate Analysis
This script examines the relationship between gene expression patterns and evolutionary rates. It involves:
- Outlier filtering
- Linear regression analysis (peak difference, mean expression levels vs dN/dS)
- Multiple comparisons (Steel-Dwass test) of dN/dS

# Usage
Each script is designed to be run in R. Ensure that all required packages are installed before executing the scripts. Detailed instructions on data input formats and dependencies are provided within each script.

