# This script concern all the steps involved in the RNAseq analysis.
# It begin with the standard analysis (quality control, alignment) and continues
# with all the analysis related to differentially expressed genes (normalization,
# statisical test, functionnal enrichment...)


######################
# 1. Quality control #
######################


# 1. Quality control for raw data
rule fastQC:
  input:
    file1=os.path.join(
        config["path"]["data"]["rna"],
        "{all_sample}_R1_001.fastq.gz"),
    file2=os.path.join(
        config["path"]["data"]["rna"],
        "{all_sample}_R2_001.fastq.gz")
  output:
    os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "01.fastQC",
        "{all_sample}_{brin}_001_fastqc.zip")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    fastqc -t {threads} -o {config[path][result][RNAseq]}/01.QualityControl/01.fastQC {input.file1} {input.file2}
    """ 

rule multiQC_raw:
  input:
    expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "01.fastQC",
        "{all_sample}_{brin}_001_fastqc.zip"), all_sample=TOTAL_RNASEQ, brin=STRAND)
  output:
    os.path.join(
            config["path"]["result"]["RNAseq"],
            "01.QualityControl",
            "01.fastQC",
            "multiqc",
            "multiqc_report.html")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    multiqc -f -o {config[path][result][RNAseq]}/01.QualityControl/01.fastQC/multiqc/ {input}
    """

# 2.Trimming. This step remove low-quality bases and sequencing adaptators
rule Trimmomatic:
  input:
    file1=os.path.join(
        config["path"]["data"]["rna"],
        "{all_sample}_R1_001.fastq.gz"),
    file2=os.path.join(
        config["path"]["data"]["rna"],
        "{all_sample}_R2_001.fastq.gz")
  output:
    file1=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R1.fq.gz"),
    file2=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R1_unpaired.fq.gz"),
    file3=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R2.fq.gz"),
    file4=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R2_unpaired.fq.gz"),
    file5 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_trim.log")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["high"]
  shell:
    """
    TrimmomaticPE -threads {threads} -phred33 <(unpigz -c {input.file1}) <(unpigz -c {input.file2}) {output.file1} {output.file2} {output.file3} {output.file4} ILLUMINACLIP:{config[trimmomatic][adaptateurs]} TRAILING:{config[trimmomatic][TRAILING]} SLIDINGWINDOW:{config[trimmomatic][SLIDINGWINDOW]} MINLEN:{config[trimmomatic][MINLEN]} 2> {output.file5}
    """ 

# 3. Quality control for trimmed data 
rule fastQC2:
  input:
    file1=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R1.fq.gz"),
    file2=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R2.fq.gz")
  output:
    os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "03.fastQC",
        "{all_sample}_{brin}_fastqc.html")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    fastqc -t {threads} -o {config[path][result][RNAseq]}/01.QualityControl/03.fastQC {input.file1} {input.file2}
    """ 


#############################
# 2. Alignment against hg19 #
#############################

# 1. Index for STAR alignment
rule index:
  output:
    os.path.join(
          config["path"]["result"]["RNAseq"],
          "02.Alignment",
          "01.Index",
          "genomeParameters.txt")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {config[path][result][RNAseq]}/02.Alignment/01.Index --genomeFastaFiles {config[ref][fasta]} \
    --sjdbGTFfile {config[ref][gtf]} --sjdbOverhang {config[STARpass][Overhang]} --limitGenomeGenerateRAM=124000000000
    """

# 2. STAR alignment
rule STAR:
  input:
    file1=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R1.fq.gz"),
    file2=os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_R2.fq.gz"),
    index =  os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "01.Index",
        "genomeParameters.txt")
  output:
    bam = os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "02.Alignment",
      "{all_sample}/Aligned.sortedByCoord.out.bam"),
    log = os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "02.Alignment",
      "{all_sample}/Log.final.out")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["high"]
  shell:
    """
    STAR --runThreadN {threads} --genomeDir {config[path][result][RNAseq]}/02.Alignment/01.Index --readFilesIn <(unpigz -c {input.file1}) <(unpigz -c {input.file2}) --outSAMtype {config[STARpass][outSAMtype]} --outReadsUnmapped {config[STARpass][outReadsUnmapped]} --outFileNamePrefix {config[path][result][RNAseq]}/02.Alignment/02.Alignment/{wildcards.all_sample}/
    """

# Alignment quality control
rule bamQC:
  input: 
    os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "02.Alignment",
      "{all_sample}/Aligned.sortedByCoord.out.bam")
  output: 
    bamQC1 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "03.QualityControl",
      "{all_sample}/raw_data_qualimapReport",
      "coverage_profile_along_genes_(high).txt"),
    bamQC2 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "03.QualityControl",
      "{all_sample}/raw_data_qualimapReport",
      "coverage_profile_along_genes_(low).txt"),
    bamQC3 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "03.QualityControl",
      "{all_sample}/raw_data_qualimapReport",
      "coverage_profile_along_genes_(total).txt"),
    bamQC4 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "03.QualityControl",
      "{all_sample}/rnaseq_qc_results.txt")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    ./{config[container][qualimap]} rnaseq -bam {input} --java-mem-size=128G -p strand-specific-forward -gtf {config[ref][gtf]} -pe -outdir {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}
    """

# Alignment quality control
rule rename_BamQC:
  input: 
    bamQC1 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport",
        "coverage_profile_along_genes_(high).txt"),
    bamQC2 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport",
        "coverage_profile_along_genes_(low).txt"),
    bamQC3 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport",
        "coverage_profile_along_genes_(total).txt")
  output: 
    bamQC1 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport_RENAME",
        "coverage_profile_along_genes_high.txt"),
    bamQC2 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport_RENAME",
        "coverage_profile_along_genes_low.txt"),
    bamQC3 = os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport_RENAME",
        "coverage_profile_along_genes_total.txt")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    cp {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}/raw_data_qualimapReport/coverage_profile_along_genes_\(high\).txt {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}/raw_data_qualimapReport_RENAME/coverage_profile_along_genes_high.txt
    cp {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}/raw_data_qualimapReport/coverage_profile_along_genes_\(low\).txt {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}/raw_data_qualimapReport_RENAME/coverage_profile_along_genes_low.txt
    cp {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}/raw_data_qualimapReport/coverage_profile_along_genes_\(total\).txt {config[path][result][RNAseq]}/02.Alignment/03.QualityControl/{wildcards.all_sample}/raw_data_qualimapReport_RENAME/coverage_profile_along_genes_total.txt
    """

###########################
# 3. Reads quantification #
###########################

rule FeatureCounts:
  input:
    expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "02.Alignment",
        "{all_sample}/Aligned.sortedByCoord.out.bam"), all_sample=SAMPLES)
  output:
    file1 = os.path.join(
         config["path"]["result"]["RNAseq"],
        "03.Quantification",
        "Reads_counts.tab"),
    file2 = os.path.join(
         config["path"]["result"]["RNAseq"],
        "03.Quantification",
        "Reads_counts.tab.summary")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["high"]
  shell:
    """
    featureCounts -G {config[ref][fasta]} -a {config[ref][gtf]} -t exon -g gene_id -J -o {output.file1} -T {threads} {input}
    """


# Final quality control

rule multiQC:
  input:
    fastqc = expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "03.fastQC",
        "{all_sample}_{brin}_fastqc.zip"), all_sample=SAMPLES, brin = STRAND),
    trimmomatic = expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "01.QualityControl",
        "02.Trimmomatic",
        "{all_sample}_trim.log"),all_sample=SAMPLES),
    bamQC1 = expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport_RENAME",
        "coverage_profile_along_genes_high.txt"),all_sample=SAMPLES),
    bamQC2 = expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport_RENAME",
        "coverage_profile_along_genes_low.txt"),all_sample=SAMPLES),
    bamQC3 = expand(os.path.join(
        config["path"]["result"]["RNAseq"],
        "02.Alignment",
        "03.QualityControl",
        "{all_sample}/raw_data_qualimapReport_RENAME",
        "coverage_profile_along_genes_total.txt"),all_sample=SAMPLES),
    bamQC4 = expand(os.path.join(
          config["path"]["result"]["RNAseq"],
          "02.Alignment",
          "03.QualityControl",
          "{all_sample}/rnaseq_qc_results.txt"),all_sample=SAMPLES),
    star = expand(os.path.join(
      config["path"]["result"]["RNAseq"],
      "02.Alignment",
      "02.Alignment",
      "{all_sample}/Log.final.out"),all_sample=SAMPLES),
    featureCounts = os.path.join(
         config["path"]["result"]["RNAseq"],
        "03.Quantification",
        "Reads_counts.tab.summary")
  output:
    os.path.join(
            config["path"]["result"]["RNAseq"],
            "01.QualityControl",
            "Total_QC",
            "multiqc_report.html")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    multiqc -f -o {config[path][result][RNAseq]}/01.QualityControl/Total_QC/ {input.fastqc} {input.trimmomatic} {input.bamQC1} {input.bamQC2} {input.bamQC3} {input.bamQC4} {input.star} {input.featureCounts}
    """


##############################################################
# 4. Reads normalization and test of differential expression #
##############################################################

rule DEG_analysis:
  input:
    file1 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "03.Quantification",
      "Reads_counts.tab"),
    file2 = config["ref"]["phenoTable"],
    file3 = config["ref"]["convTable"]
  output:
    output1 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "04.Differential_expression",
      "{test}",
      "Normalized_data.tab"),
    output2 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "04.Differential_expression",
      "{test}",
      "DEGs_test.tab")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    Rscript 04.Scripts/RNAseq_normalization_DEGs.R {config[path][result][RNAseq]} {input.file1} {input.file2} {input.file3} {config[DESeq2][FCcutoff]} {wildcards.test}
    """



######################
# 5. lncRNA analysis #
######################

rule DL_databases:
  output:
    lncRNA2Target = os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "00.Databases",
      "lncRNA2Target.csv"),
    lncTarD = os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "00.Databases",
      "lncTarD.txt"),
    lncRNAdisease = os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "00.Databases",
      "lncRNAdisease.csv")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    # 1. lncRNA2Target database
    wget -P {config[path][result][RNAseq]}/05.lncRNA/00.Databases http://bio-annotation.cn/lncrna2target/data/lncRNA_target_from_low_throughput_experiments.xlsx
    xlsx2csv -d 'tab' {config[path][result][RNAseq]}/05.lncRNA/00.Databases/lncRNA_target_from_low_throughput_experiments.xlsx > {output.lncRNA2Target}
    sed -i 's/"//g' {output.lncRNA2Target}
    sed -i "s/'//g" {output.lncRNA2Target}

    # 2. lncTarD
    wget -P {config[path][result][RNAseq]}/05.lncRNA/00.Databases http://bio-bigdata.hrbmu.edu.cn/LncTarD/download/lncTarD.txt
    sed -i 's/"//g' {output.lncTarD}
    sed -i "s/'//g" {output.lncTarD}

    # 3. lncRNAdisease
    wget -P {config[path][result][RNAseq]}/05.lncRNA/00.Databases -O {config[path][result][RNAseq]}/05.lncRNA/00.Databases/lncRNAdisease.xlsx http://www.rnanut.net/lncrnadisease/static/download/experimental%20lncRNA-disease%20information.xlsx
    xlsx2csv -d 'tab' {config[path][result][RNAseq]}/05.lncRNA/00.Databases/lncRNAdisease.xlsx > {output.lncRNAdisease}
    sed -i 's/"//g' {output.lncRNAdisease}
    sed -i "s/'//g" {output.lncRNAdisease}
    """

rule ncRNA_analysis:
  input:
    normData = os.path.join(
      config["path"]["result"]["RNAseq"],
      "04.Differential_expression",
      "{test}",
      "Normalized_data.tab"),
    DEGs = os.path.join(
      config["path"]["result"]["RNAseq"],
      "04.Differential_expression",
      "{test}",
      "DEGs_test.tab"),
    phenoTable = config["ref"]["phenoTable"],
    lncRNA2Target = os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "00.Databases",
      "lncRNA2Target.csv"),
    lncTarD = os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "00.Databases",
      "lncTarD.txt"),
    lncRNAdisease = os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "00.Databases",
      "lncRNAdisease.csv")
  output:
    os.path.join(
      config["path"]["result"]["RNAseq"],
      "05.lncRNA",
      "{test}",
      "lncRNAdisease_enrichment.tab")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    Rscript 04.Scripts/ncRNA_analysis.R {config[path][result][RNAseq]}/05.lncRNA/{wildcards.test} {input.normData} {input.DEGs} {input.phenoTable} {config[DESeq2][FCcutoff]} {input.lncRNA2Target} {input.lncTarD} {input.lncRNAdisease}
    """


###########################
# 6. Cell type proportion #
###########################

rule dlDB:
  output:
    file1 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "GSE109816_normal_heart_umi_matrix.csv"),
    file2 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "GSE109816_normal_heart_cell_info.txt"),
    file3 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "GSE109816_normal_heart_cell_cluster_info.txt")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    # 1. Download files for heart cells dabatase
    mkdir -p {config[path][result][RNAseq]}/06.cellTypeProportion/00.Database
    wget -P {config[path][result][RNAseq]}/06.cellTypeProportion/00.Database https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109816/suppl/GSE109816_normal_heart_umi_matrix.csv.gz
    wget -P {config[path][result][RNAseq]}/06.cellTypeProportion/00.Database https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109816/suppl/GSE109816_normal_heart_cell_info.txt.gz
    wget -P {config[path][result][RNAseq]}/06.cellTypeProportion/00.Database https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109816/suppl/GSE109816_normal_heart_cell_cluster_info.txt.gz

    # 2. Decompress all files
    unpigz {config[path][result][RNAseq]}/06.cellTypeProportion/00.Database/*
    """

rule buildDB:
  input:
    file1 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "GSE109816_normal_heart_umi_matrix.csv"),
    file2 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "GSE109816_normal_heart_cell_info.txt"),
    file3 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "GSE109816_normal_heart_cell_cluster_info.txt")
  output:
    os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "heartDB.tab")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["medium"]
  shell:
    """
    Rscript 04.Scripts/buildDB.R {config[path][result][RNAseq]}/06.cellTypeProportion
    """

rule cellType:
  input:
    file1 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "03.Quantification",
      "Reads_counts.tab"),
    file2 = config["ref"]["phenoTable"],
    file3 = config["ref"]["convTable"],
    file4 = os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "00.Database",
      "heartDB.tab")
  output:
    os.path.join(
      config["path"]["result"]["RNAseq"],
      "06.cellTypeProportion",
      "{test}",
      "All_enrichment.tab")
  singularity:
    config["container"]["RNAseq"]
  threads:
    config["threads"]["low"]
  shell:
    """
    Rscript 04.Scripts/cellProportion.R {config[path][result][RNAseq]}/06.cellTypeProportion {input.file1} {input.file2} {input.file3} {wildcards.test}
    """