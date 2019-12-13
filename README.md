# The upper respiratory tract microbiome of patients with chronic rhinosinusitis

The goal of this project is to compare the microbiome of the human upper respiratory tract between patients with chronic rhinosinusitis and healthy controls. This will allow us to assess whether the microbiome plays a role in this disease and if yes, what bacteria are potentially associated with health or disease.

## The data

The data consists of 16S rRNA amplicon sequenced samples (V4 region), sequenced in seven different MiSeq runs.

## Analyses

### General analysis of the data

The RMarkdown document src/data_analysis.Rmd contains the main overall analysis of the data:

1) A comparison between sampling locations (nose, nasopharynx, maxillary sinus and ethmoid sinus)
2) A comparison between samples from CRS patients and healthy controls, for the nose and nasopharynx samples
3) Association tests between the microbiome data of the CRS patients and some patient/disease related covariates

This analysis has been published in mSphere:

[De Boeck I, Wittouck S, Martens K, Claes J, Jorissen M, Steelant B, van den Broek MFL, Seys SF, Hellings PW, Vanderveken OM, Lebeer S. 2019. Anterior nares diversity and pathobionts represent sinus microbiome in chronic rhinosinusitis. mSphere 4:e00532-19. https://doi.org/10.1128/mSphere.00532-19.](https://doi.org/10.1128/mSphere.00532-19)

### Family-level analyses (CODA and Lactobacillus)

These analyses can be found in the folder src/family-level. 

We performed two types of analyses with the read counts aggregated on the family level: 

* A compositional data analysis (CODA) to test for differential abundance of families between CRS and healthy control samples. 
* A differential abundance test for the Lactobacillus Genus Complex (our favourite bacteria) on the relative abundances. 

These analyses are reported in a manuscript that is currently being peer reviewed. 
