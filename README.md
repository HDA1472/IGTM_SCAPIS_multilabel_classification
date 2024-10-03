# Multi-Label Classification in SCAPIS & IGTM Cohorts

## Description

This repo contains the code generated in the context of the study “Plasma Proteomics-Based 
Multi-Label Classification of Co-occurring Metabolic Diseases”, where we built and evaluated
different machine learning models to predict the presence of multiple metabolic and cardiovascular
diseases in two independent general Swedish population cohorts.

## Study

While current studies focus mostly on case-control comparisons, the complex interplay and 
co-occurrence of cardiovascular and metabolic diseases necessitate the development of robust 
and sophisticated multi-label classification methods. This study aims to identify novel 
biomarker signatures that will facilitate early diagnosis and improved stratification of 
patients with multiple metabolic and cardiovascular conditions using advanced machine 
learning techniques.

In total, 3,000 plasma samples were analyzed using the Olink Explore 1536 platform, a 
highly sensitive and multiplexed antibody-based technology. This dataset consists of 
two independent cohorts of patients aged 50-65 from the general Swedish population that 
will be used as discovery (n=2,000) and validation (n=1,000) cohorts. We employed and 
evaluated the performance of a range of machine learning methods, starting with binary 
lasso classifiers as a baseline and expanding to more complex approaches including 
Classifier Chains, Neural Networks, and unsupervised techniques.

By applying our multifaceted machine learning approach to a comprehensive cardiometabolic 
plasma proteomics dataset, we identified novel biomarker signatures for multiple co-occurring 
metabolic and cardiovascular conditions. By integrating these biomarker signatures with 
clinical metadata, we defined phenotypic subtypes within our multi-diseased cohort. This 
provides a deeper understanding of co-occurring disease patterns and uncovers new insights 
into disease mechanisms and interactions. Ultimately, our findings offer insights into disease 
mechanisms and pave the way for personalized therapeutic strategies.

## Content
This repository includes the code to generate the results described above. It is organized as follows:

- .qmd files: the Quarto files contain the code to generate the results.

- Code.Rproj: R project file.

