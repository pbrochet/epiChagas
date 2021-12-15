# This script concern all the steps involved in tissue methylation and transcription factor analysis
# It begin with the methylation (quality control, normalization, batch effect correction, 
# differentially methylated position/region test) and continues
# with all the analysis related to transcription factor analysis (ReMap constrution,
# OLOGRAM pairwise and Nwise analysis, TFBS enrichment with JASPAR...)


rule Tissue_methylation:
  input:
    sampleSheet = os.path.join(
      config["path"]["data"]["methylation"]["tissue"],
      "SampleSheet.csv"),
    phenoTable = config["ref"]["phenoTable"],
    conversionTable = config["ref"]["convTable"],
    DEGs_test = os.path.join(
      config["path"]["result"]["RNAseq"],
      "04.Differential_expression",
      "CTRL_sevCCC",
      "DEGs_test.tab"),
    ReMAP_db = config["ref"]["ReMAP"],
    GenePromoter = config["ref"]["genePromoter"]
  output:
    output1 = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "DifferentiallyMethylated_RegRegion.bed"),
    output2 = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion.bed"),
    output4=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "05.DMR",
      "myDMR_DEG.bed")
  singularity:
    config["container"]["Methylation"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    Rscript 04.Scripts/Tissue_methylation_analysis.R {config[path][data][methylation][tissue]} {config[path][result][Tissue_methylation]} {input.phenoTable} {input.conversionTable} {config[Methylation_tissue][FCcutoff]} {input.DEGs_test} {config[DESeq2][FCcutoff]} {input.ReMAP_db} {input.GenePromoter}
    """

rule parseReMAP:
  input:
    config["ref"]["ReMAP"],
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  singularity:
    config["container"]["Methylation"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    # 1. Select the unique TF present in ReMap database selection
    cat {input} | awk -F"\t" '{{print $4}}' | sort -u > {config[path][result][Tissue_methylation]}/06.TF/listeTF_subset.txt
    # 2. Create a bed file for each TF
    mkdir -p {config[path][result][Tissue_methylation]}/06.TF/All_ReMap_TF
    sh 04.Scripts/parseReMAP.sh {config[path][result][Tissue_methylation]}/06.TF/listeTF_subset.txt {config[path][result][Tissue_methylation]}/06.TF/All_ReMap_TF {config[ref][ReMAP]}
    touch {output}
    """

rule OlogramPairwise1:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "DifferentiallyMethylated_RegRegion.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "01.OlogramPairwise1",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/All_ReMap_TF/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/01.OlogramPairwise1 --force-chrom-peak --force-chrom-more-bed \
      -k {threads} -mn 10 -ms 10 -V 3 \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/01.OlogramPairwise1/*.tsv {config[path][result][Tissue_methylation]}/06.TF/01.OlogramPairwise1/Ologram.tsv
    """

rule OlogramPairwise2:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "05.DMR",
      "myDMR_DEG.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "02.OlogramPairwise2",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/All_ReMap_TF/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/02.OlogramPairwise2 --force-chrom-peak --force-chrom-more-bed \
      -k {threads} -mn 10 -ms 10 -V 3 \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/02.OlogramPairwise2/*.tsv {config[path][result][Tissue_methylation]}/06.TF/02.OlogramPairwise2/Ologram.tsv
    """

rule OlogramPairwise3:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "05.DMR",
      "myDMR_DEG.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "DifferentiallyMethylated_RegRegion.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "03.OlogramPairwise3",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/All_ReMap_TF/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/03.OlogramPairwise3 --force-chrom-peak --force-chrom-more-bed \
      -k {threads} -mn 10 -ms 10 -V 3 \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/03.OlogramPairwise3/*.tsv {config[path][result][Tissue_methylation]}/06.TF/03.OlogramPairwise3/Ologram.tsv
    """

rule OlogramPairwise4:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "05.DMR",
      "myDMR_DEG.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion2.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "04.OlogramPairwise4",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/All_ReMap_TF/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/04.OlogramPairwise4 --force-chrom-peak --force-chrom-more-bed \
      -k {threads} -mn 10 -ms 10 -V 3 \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/04.OlogramPairwise4/*.tsv {config[path][result][Tissue_methylation]}/06.TF/04.OlogramPairwise4/Ologram.tsv
    """ 
# rule parseMethylation:
#   input:
#     listeTF = os.path.join(
#       config["path"]["result"]["Tissue_methylation"],
#       "06.TF",
#       "listeTF_subset2.txt")
#   output:
#     os.path.join(
#       config["path"]["result"]["Tissue_methylation"],
#       "06.TF",
#       "myTFs",
#       "All_done.txt")
#   singularity:
#     config["container"]["Methylation"]
#   threads:
#     config["threads"]["medium"]
#   shell:
#     """
#     # 1. We subset the TFs binding at least 50% of our genes having a DMR
#     # bedtools intersect -wa -wb -b {input.DMRs} -a {config[ref][ReMAP]} > {config[path][result][Tissue_methylation]}/06.TF/TFs_Genes.tab
#     # Rscript 04.Scripts/subsetTF.R {config[path][result][Tissue_methylation]}/06.TF/TFs_Genes.tab {input.listeTF} {input.DMRs} {config[path][result][Tissue_methylation]}/06.TF
#     # 2. Parse the final TF for ologram analysis
#     mkdir -p {config[path][result][Tissue_methylation]}/06.TF/myTFs
#     sh 04.Scripts/parseReMAP.sh {input.listeTF} {config[path][result][Tissue_methylation]}/06.TF/myTFs {config[ref][ReMAP]}
#     touch {output}
#     """


rule OlogramPairwise5:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "DifferentiallyMethylated_RegRegion.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "05.OlogramPairwise5",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/testsRegions/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/05.OlogramPairwise5 --force-chrom-peak --force-chrom-more-bed \
      -k {threads} -mn 10 -ms 10 -V 3 \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/05.OlogramPairwise5/*.tsv {config[path][result][Tissue_methylation]}/06.TF/05.OlogramPairwise5/Ologram.tsv
    """ 

rule OlogramPairwise6:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "05.DMR",
      "myDMR_DEG.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "All_ReMap_TF",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "06.OlogramPairwise6",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/testsRegions/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/06.OlogramPairwise6 --force-chrom-peak --force-chrom-more-bed \
      -k {threads} -mn 10 -ms 10 -V 3 \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/06.OlogramPairwise6/*.tsv {config[path][result][Tissue_methylation]}/06.TF/06.OlogramPairwise6/Ologram.tsv
    """ 

rule parseMethylation:
  input:
    listeTF = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "listeTF_subset2.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "myTFs",
      "All_done.txt")
  singularity:
    config["container"]["Methylation"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    # 2. Parse the final TF for ologram analysis
    mkdir -p {config[path][result][Tissue_methylation]}/06.TF/myTFs
    sh 04.Scripts/parseReMAP.sh {input.listeTF} {config[path][result][Tissue_methylation]}/06.TF/myTFs {config[ref][ReMAP]}
    touch {output}
    """

rule OlogramNwise:
  input:
    region = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "05.DMR",
      "myDMR_DEG.bed"),
    control = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "otherRegRegion.bed"),
    All_TFs=os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "myTFs",
      "All_done.txt")
  output:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "05.OlogramNwise",
      "Ologram.tsv")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    sed -i 's/"//g' {input.region}
    sed -i 's/"//g' {input.control}
    gtftk ologram -z -c hg19 -p {input.region} --more-bed `ls {config[path][result][Tissue_methylation]}/06.TF/myTFs/*.bed` \
      -o {config[path][result][Tissue_methylation]}/06.TF/05.OlogramNwise --force-chrom-peak --force-chrom-more-bed \
      -k 24 -mn 10 -ms 10 -V 3  --more-bed-multiple-overlap \
      --bed-incl {input.control}
    mv {config[path][result][Tissue_methylation]}/06.TF/05.OlogramNwise/*.tsv {config[path][result][Tissue_methylation]}/06.TF/05.OlogramNwise/Ologram.tsv
    """

rule Ologram_graph:
  input:
    os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "05.OlogramNwise",
      "Ologram.tsv")
  output:
    total = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "05.OlogramNwise",
      "Ologram_totalGraph.pdf"),
    final = os.path.join(
      config["path"]["result"]["Tissue_methylation"],
      "06.TF",
      "05.OlogramNwise",
      "Ologram_finalGraph.pdf")
  singularity:
    config["container"]["Ologram"]
  threads:
    config["threads"]["low"]
  shell:
    """
    # 1. All complex tree graph
    gtftk ologram_modl_treeify -i {input} -o {output.total}
    # 2. Subset tree graph
    Rscript 04.Scripts/parseOlogram.R {input} {config[path][result][Tissue_methylation]}/06.TF/05.OlogramNwise {config[TFs][Ologram_complexSize]}
    gtftk ologram_modl_treeify -i {config[path][result][Tissue_methylation]}/06.TF/05.OlogramNwise/Ologram_subset.tsv -o {output.final} -mh {config[TFs][Ologram_heritance]}
    """