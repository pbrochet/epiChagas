# Configuration file
configfile: "02.Config/config.yaml"

# Rules
include: "03.Rules/01.RNAsequencing_analysis.smk"
include: "03.Rules/02.Tissue_methylation_analysis.smk"
include: "03.Rules/03.Blood_methylation_analysis.smk"

# Variables list
SAMPLES=["sevCCC1", "sevCCC2", "sevCCC3", "sevCCC4", "sevCCC5", "sevCCC6", "sevCCC7", "sevCCC8", "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5", "CTRL6", "DCM1", "DCM2", "DCM3", "DCM4", "DCM5", "DCM6", "DCM7", "DCM8"]
TEST=["CTRL_sevCCC", "DCM_sevCCC", "CTRL_DCM"]
STRAND=["R1", "R2"]
TRIMMO=["paired", "unpaired"]

# All steps
rule all:
    input:
        ### RNA-seq steps
            # Include all fastQC analysis on raw data
        multiqc_raw = expand("{path}/01.QualityControl/01.fastQC/multiqc/multiqc_report.html", path = config["path"]["result"]["RNAseq"])
            # Include quality control of all preprocessed steps : trimming, alignment, gene quantification
        multiqc = expand("{path}/01.QualityControl/Total_QC/multiqc_report.html", path = config["path"]["result"]["RNAseq"]),
            # Not necessary, include gene normalization and differential expression test, for each phenotype of interest
        degs = expand("{path}/04.Differential_expression/{test}/DEGs_test.tab", path = config["path"]["result"]["RNAseq"], test=TEST),
            # Deconvolution analysis to identify cell proportion in heart tissue samples
        cellProportion = expand("{path}/06.cellTypeProportion/{test}/All_enrichment.tab", path = config["path"]["result"]["RNAseq"], test = TEST),

        ### Tissue methylation steps
            # Include all tissue methylation steps : quality control, normalization, batch effect correction, differential methylation position test, differential methylation region test (necessite DEGs results)
        dmr = expand("{path}/05.DMR/myDMR_DEG.bed", path = config["path"]["result"]["Tissue_methylation"]),
       
        ### Transcription factor analysis steps
            # Include ReMap database subsetting, OLOGRAM pairwise and OLOGRAM N-wise analysis 
        ologramGraph = expand("{path}/06.TF/05.OlogramNwise/Ologram_finalGraph.pdf", path = config["path"]["result"]["Tissue_methylation"]),

        ### Need to add JASPAR analysis ###

        ### Blood methylation steps
            # Include all blood methylation steps : quality control, normalization, batch effect correction, differential methylation position test
        dmp_blood = expand("{path}/04.DMP/modCCC_sevCCC.tab", path = config["path"]["result"]["Blood"]["Blood_methylation"])

        ### Need to add biomarker analysis ###
