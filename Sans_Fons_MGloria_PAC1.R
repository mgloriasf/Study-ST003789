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

#Observació de les dades

ncol(m1)
class(m1)
head(m1)
tail(m1)
dim(m1)
colData(m1)
dplyr::glimpse(m1)

#Anàlisis de les dades

BiocManager::install("POMA",force=TRUE)
library(POMA)
library(ggtext)
library(patchwork)
library(magrittr)
library(tidyverse)


is.na(se)


#treure valors na

class(m1)
dim(m1)
imputed<-PomaImpute(m1,method="knn", zeros_as_na=TRUE, remove_na=TRUE,cutoff=20)
imputed
normalized<-PomaNorm(imputed,method="log_pareto")
normalized
PomaBoxplots(normalized,x="samples")+
  ggtitle("normalized")+
  theme(legend.position="none")
PomaBoxplots(m1,x="samples")
PomaDensity(m1,x="samples")


PomaOutliers(normalized,do=clean)
do<-PomaOutliers(normalized,
             method="euclidean",
             type="median",
             outcome="Diet",
             coef=2,
             labels=FALSE
)
do$polygon_plot
PomaOutliers(normalized,
                              method="euclidean",
                              type="median",
                              outcome="Treatment",
                              coef=2,
                              labels=FALSE
)


PomaOutliers(outlier_results,do="clean")

PomaUnivariate(normalized,method="ttest",covs= Null,error=NULL, paired=FALSE,adjust="fdr")


PomaPCA(
  normalized,
  outcome="Treatment",
  center=TRUE,
  scale=TRUE,
  ncomp=5,
  labels=FALSE,
  ellipse=FALSE,
  load_length=1
)


PomaPCA(
  normalized,
  outcome="Diet",
  center=TRUE,
  scale=TRUE,
  ncomp=5,
  labels=FALSE,
  ellipse=FALSE,
  load_length=1
)




poma_cor <- PomaCorr(normalized)
poma_cor$correlations
poma_cor$corrplot

BiocManager::install("PRONE")
library(PRONE)
plot_condition_overview(normalized, condition="Diet")
plot_condition_overview(normalized, condition="Treatment")

plot_NA_heatmap(normalized, color_by=NULL,label_by=NULL, cluster_samples=TRUE,cluster_proteins= TRUE)

