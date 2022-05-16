
# Epigenetic regulation of transcription factor binding motifs promotes Th1 response in Chagas disease Cardiomyopathy

## Article information

Title : Specific methylation marks in promoter regions are associated to the pathogenic process of Chronic Chagas disease Cardiomyopathy by modifying transcription factor binding patterns

Authors : Pauline Brochet1, Barbara Ianni2, Laurie Laugier3, Amanda Farage Frade2,4,5, João Paulo Silva Nunes1,2,4,5, Priscila Camillo Teixeira2,4,5, Charles Mady6, Ludmila Rodrigues Pinto Ferreira7, Quentin Ferré1, Ronaldo Honorato Barros Santos8, Andreia Kuramoto2 , Sandrine Cabantous3, Samuel Steffen8,9, Antonio Noedir Stolf9, Pablo Pomerantzeff10, Alfredo Inacio Fiorelli9, Edimar A Bocchi9, Cristina Wide Pissetti11, Bruno Saba12, Darlan da Silva Cândido2,4,5, Fabrício Dias13, Marcelo Sampaio12, Fabio Antônio Gaiotto8,9, José Antonio Marin-Neto13, Abílio Fragata12, Ricardo Costa Fernandes Zaniratto2, Sergio Siqueira14, Giselle de lima Peixoto14, Vagner Oliveira-Carvalho Rigaud2,15, Fernando Bacal8, Paula Buck10, Rafael Almeida Ribeiro2,4,5, Hui Tzu Lin-Wang12, José Antonio Marin-Neto13, André Schmidt13, Martino Martinelli14, Mario Hiroyuki Hirata16, Eduardo Donadi13, Alexandre Costa Pereira10, Virmondes Rodrigues Junior11, Denis Puthier1, Benoit Ballester1, Pierre Milpied17, Jorge Kalil2,4,5, Lionel Spinelli1,17,c, Edecio Cunha-Neto2,4,5,*,c, Christophe Chevillard1,*,c. 

* Authors had an equal contribution. 
c corresponding authors 

 
Affiliations: 

1. INSERM, UMR_1090, Aix Marseille Université, TAGC Theories and Approaches of Genomic Complexity, Institut MarMaRa, Marseille, France. 

2. Laboratory of Immunology, Heart Institute (InCor), University of São Paulo, School of Medicine, São Paulo, Brazil. 

3. Aix Marseille Université, Génétique et Immunologie des Maladies Parasitaires, Inserm, UMR_906, Marseille, France. 

4. Division of Clinical Immunology and Allergy, University of São Paulo, School of Medicine, São Paulo, Brazil. 

5. Instituto Nacional de Ciência e Tecnologia, INCT, iii- Institute for Investigation in Immunology, São Paulo, Brazil. 

6.  Myocardiopathies Unit, Heart Institute (InCor), School of Medicine, University of São Paulo, São Paulo, Brazil.  

7. RNA Systems Biology Laboratory (RSBL), Departamento de Morfologia, Instituto de Ciências Biológicas, Universidade Federal de Minas Gerais, Belo Horizonte, Minas Gerais, Brazil. 

8. Division of Transplantation, Heart Institute (InCor), University of São Paulo, School of Medicine, São Paulo, Brazil. 

9.  Division of Surgery, Heart Institute (InCor), University of São Paulo, School of Medicine, São Paulo, Brazil. 

10. Heart Institute (InCor), School of Medicine, University of São Paulo, São Paulo, São Paulo, Brazil. 

11.  Laboratory of Immunology, Universidade Federal Do Triângulo Mineiro (UFTM), Uberaba, Brazil. 

12. Laboratório de Investigação Molecular em Cardiologia, Instituto de Cardiologia Dante Pazzanese (IDPC), São Paulo, Brazil. 

13.  School of Medicine of Ribeirão Preto (FMRP), University of São Paulo, Ribeirão Preto, Brazil. 

14. Pacemaker Clinic, Heart Institute (InCor), School of Medicine, University of São Paulo, São Paulo, Brazil. 

15. Heart Failure Unit, Heart Institute (InCor) School of Medicine, University of Sao Paulo, Sao Paulo, Brazil. 

16. Department of Clinical and Toxicological Analyses, Faculty of Pharmaceutical Sciences, University of São Paulo (USP), São Paulo, Brazil. 

17. Aix Marseille Univ, CNRS, Inserm, Centre d'Immunologie de Marseille-Luminy, Marseille, France. 

 

Name and complete address for correspondence: 

Edecio Cunha Neto, MD, PhD 

Laboratório de Imunologia, Instituto do Coração, Hospital das Clínicas, Faculdade de Medicina da Universidade de São Paulo. 

Av. Dr. Eneas de Carvalho Aguiar, 44 - Bloco II, 9º andar - 05403-000 São Paulo, SP, Brasil 

Fax number: +55 11 2661-5953 / Telephone number: +55 11 2661-5914/5906 / E-mail address: edecunha@gmail.com 

Christophe Chevillard, PhD 


Abstract :


Chagas disease, caused by the protozoan Trypanosoma cruzi, is an endemic parasitic disease of Latin America, affecting 7 million people. Although most patients are asymptomatic, 30% develop complications, including the often-fatal Chronic Chagasic Cardiomyopathy (CCC). The pathogenic process remains poorly understood. 

Based on bulk RNA-seq and whole genome DNA methylation data, we investigated the genetic and epigenetic deregulations present in the moderate and severe stages of CCC. Based on heart tissue gene expression profile, we had identified 1407 differentially expressed transcripts (DEGs) specific from CCC patients. A tissue DNA methylation analysis done on the same tissue has permitted the identification of 92 regulatory Differentially Methylated Regions (DMR) localized in the promoter of DEGs. An in-depth study of the transcription factors binding sites (TFBS) in the DMRs corroborated the importance of TFBS’s DNA methylation for gene expression in CCC myocardium. TBX21, RUNX3 and EBF1 are the transcription factors whose binding motif appears to be affected by DNA methylation in the largest number of genes. 

By combining both transcriptomic and methylomic analysis on heart tissue, and methylomic analysis on blood, 4 biological processes affected by severe CCC have been identified, including immune response, ion transport, cardiac muscle processes and nervous system. An additional study on blood methylation of moderate CCC samples put forward the importance of ion transport and nervous system in the development of the disease.



## Goal of the repository

This GitHub repository contains the instructions and material to reproduce the analysis performed on *Specific methylation marks in promoter regions are associated to the pathogenic process of Chronic Chagas disease Cardiomyopathy by modifying transcription factor binding patterns* paper [DOI]. Builded Docker/Singularity images are available on download and required data may be downloaded from their original sources. Instructions necessary to reproduce the analysis are provided below.

To build the database, you first need to prepare the environments and then follow the steps described there.



## Description of the datasets

Three datasets are used in this study :

- Heart tissue RNA-seq data, from 6 control, 8 dilated cardiomyopathies and 8 severe CCC samples.

- Heart tissue methylation, from 6 control and 8 severe CCC samples.

- Heart tissue RNA-seq data, from 48 asymptomatic, 48 moderate and 48 severe CCC dilated cardiomyopathies samples.


All those data are available at [GEO].



## System requirement and dependencies

All the analysis has been done on a Linux system using Docker and Singularity containers. We recommend to use systems with at least 24 threads and 192 GB of RAM available. A stable connection to the Internet is also required as some information are queried from different databases accessible online.



## Environment preparation

In order to prepare the environment for analysis execution, it is required to:

1. Clone the current GitHub repository and set the `WORKING_DIR` environment variable
2. Download the Docker image (.tar.gz) and Singularity image (.img) files
3. Install [Docker](https://www.docker.com/), [Docker-compose](https://docs.docker.com/compose/ and [Singularity](https://singularity.lbl.gov/)
4. Load the Docker images on your system and start the containers
 
This section provides additional information for each of these steps.
 

### Clone the GitHub repository

Use you favorite method to clone this repository in a chosen folder (see [GitHub documentation](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) for more information). This will create a folder called `EpiChagas` containing all the code and documentation. 

Then, set an environment variable called `WORKING_DIR` with a value set to the path to this folder. For instance, if you cloned the repository in `/home/user/workspace`, then the `WORKING_DIR` variable needs be set to `/home/user/workspace/EpiChagas`.

**On Linux:**

```
export WORKING_DIR=/home/user/workspace/EpiChagas
```
 

### Download the Docker and Singularity images

### Load Docker images on the system

## Additional information

### Documentation

### Dates of download

### Tree view


```
.
├── 01.Container
│   ├── DockerToSingularity.sh
│   ├── Methylation
│   │   └── Dockerfile
│   ├── Ologram
│   │   ├── Dockerfile
│   │   └── pygtftk
│   └── RNAseq
│       ├── Dockerfile
│       └── qualimap-build-31-08-20
├── 02.Config
│   ├── clusterConfig.json
│   └── config.yaml
├── 03.Rules
│   ├── 01.RNAsequencing_analysis.smk
│   ├── 02.Tissue_methylation_analysis.smk
│   └── 03.Blood_methylation_analysis.smk
├── 04.Scripts
│   ├── Blood_methylation_analysis.R
│   ├── buildDB.R
│   ├── cellProportion.R
│   ├── GeneOntology.R
│   ├── Heatmap.R
│   ├── Kegg_pathways.R
│   ├── ncRNA_analysis.R
│   ├── parseOlogram.R
│   ├── parseOlogramResults.R
│   ├── parseReMAP.sh
│   ├── PCA.R
│   ├── RNAseq_normalization_DEGs.R
│   ├── subsetTF.R
│   ├── Tissue_methylation_analysis.R
│   └── volcanoPlot.R
├── README.md
├── runMeso.sh
├── run.sh
├── Snakefile

```
