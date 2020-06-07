# The microbiome of chronic rhinosinusitis

The goal of this project is to compare the microbiome of the human upper respiratory tract between patients with chronic rhinosinusitis (CRS) and healthy controls. This will allow us to assess whether the microbiome plays a role in this disease and if yes, what bacteria are potentially associated with health or disease.

## The data

The data consists of 16S rRNA amplicon sequenced samples (V4 region), sequenced in seven different MiSeq runs.

## Analyses

### General analysis of the data

The RMarkdown document src/data_analysis.Rmd contains the main overall analysis of the data:

1) A comparison between sampling sites (nose, nasopharynx, maxillary sinus and ethmoid sinus) for the CRS samples
2) A comparison between samples from CRS patients and healthy controls, for the nose and nasopharynx sampling sites
3) Association tests between the microbiome data of the CRS patients and some patient/disease related covariates

The results of this analysis have been published in mSphere:

[De Boeck I, Wittouck S, Martens K, Claes J, Jorissen M, Steelant B, van den Broek MFL, Seys SF, Hellings PW, Vanderveken OM, Lebeer S. 2019. Anterior nares diversity and pathobionts represent sinus microbiome in chronic rhinosinusitis. mSphere 4:e00532-19. https://doi.org/10.1128/mSphere.00532-19.](https://doi.org/10.1128/mSphere.00532-19)

### Family-level analyses (CODA and Lactobacillus)

These analyses can be found in the folder src/family-level. 

We performed two types of analyses with the read counts aggregated on the family level: 

* A compositional data analysis (CODA) to test for differential abundance of families between CRS and healthy control samples. 
* A differential abundance test for the Lactobacillus Genus Complex (our favourite bacteria) on the relative abundances. 

The results of this analysis are described in a publication in Cell Reports: 

[De Boeck, I., van den Broek, M. F. L., Allonsius, C. N., Spacova, I., Wittouck, S., Martens, K., Wuyts, S., Cauwenberghs, E., Jokicevic, K., Vandenheuvel, D., Eilers, T., Lemarcq, M., De Rudder, C., Thys, S., Timmermans, J.-P., Vroegop, A. V, Verplaetse, A., Van de Wiele, T., Kiekens, F., … Lebeer, S. (2020). Lactobacilli Have a Niche in the Human Nose. Cell Reports, 31(8), 107674. https://doi.org/https://doi.org/10.1016/j.celrep.2020.107674](https://doi.org/https://doi.org/10.1016/j.celrep.2020.107674)
