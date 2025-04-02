install.packages("Biobacmanager")
BiocManager::install("metabolomicsWorkbenchR")
BiocManager::install("SEtools")
library(SummarizedExperiment)
library(SEtools)

#S'obté l'objecte SummarizedExperiment a partir del paquet "metabolomicsWorkbenchR". 
#Permet obtenir directament aquest objecte però també podria crear-se a partir de les
#dades.
se=do_query(
  context='study',
  input_item="study_id",
  input_value="ST003789",
  output_item="SummarizedExperiment"
)

# Existeixen dos "summarizedExperiment" associats a la referència ST003789. 
#Corresponen a dos objectes que contenen dades de metaboits obtinguts per detecció de 
#metabolits mitjançant canal iònic positiu o negatiu.

m1<-se$AN006228
m1
m2<-se$AN006229
m2