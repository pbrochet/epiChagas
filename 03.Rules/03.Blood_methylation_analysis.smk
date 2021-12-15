# This script concern all the steps involved in blood methylation and biomarker analysis
# It begin with the methylation (quality control, normalization, batch effect correction, 
# differentially methylated position test) and continues
# with the biomarker research

rule Blood_methylation:
  input:
    file1 = os.path.join(
      config["path"]["data"]["methylation"]["blood"],
      "SampleSheet.csv"),
    file2 = config["ref"]["phenoTable"],
    file3 = config["ref"]["convTable"]
  output:
    os.path.join(
      config["path"]["result"]["Blood"]["Blood_methylation"],
      "04.DMP",
      "modCCC_sevCCC.tab")
  singularity:
    config["container"]["Methylation"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    Rscript 04.Scripts/Blood_methylation_analysis.R {config[path][data][methylation][blood]} {config[path][result][Blood][Blood_methylation]} {input.file2} {input.file3}
    """

