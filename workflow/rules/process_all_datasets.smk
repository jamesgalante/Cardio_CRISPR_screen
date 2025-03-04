
# Main rule file for processing all datasets (raw data -> matrix)

# Sean M. Wu Lab Dataset(s) (smu)
# Combined lineage tracing and scRNA-seq reveals unexpected first heart field predominance of human iPSC differentiation

# Loading in files from GEO
rule smu_download_geo_files:
  output:
    # Files from GEO
    metadata = "results/sean_m_wu_lab/data/GSE202398_scRNA-seq-Run_Sample_Descriptions.xlsx",
    data = "results/sean_m_wu_lab/data/GSE202398_RAW.tar"
  shell:
    """
    # Make directory to put files in
    mkdir -p results/sean_m_wu_lab/data/
    
    # Grab information file from GEO
    wget -O {output.metadata} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202398/suppl/GSE202398%5FscRNA%2Dseq%2DRun%5FSample%5FDescriptions.xlsx
    
    # Grab tar file with data from GEO
    wget -O {output.data} "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE202398&format=file"
    """

# Untar the RAW file
rule smu_untar_GSE_file:
  input:
    tar_file = "results/sean_m_wu_lab/data/GSE202398_RAW.tar",
    sample_metadata = "results/sean_m_wu_lab/data/GSE202398_scRNA-seq-Run_Sample_Descriptions.xlsx"
  output:
    processed_data = "results/sean_m_wu_lab/raw_counts.rds",
    day30_cells = "results/sean_m_wu_lab/day30_cells.rds"
  log: "results/logs/smu_untar_GSE_file.log"
  conda:
    "../envs/sean_m_wu.yml"
  resources:
    mem = "96G",
    time = "2:00:00"
  script:
    "../scripts/process_all_datasets/smu_untar_GSE_file.R"
    

### ========================================================================

rule download_and_index_hg38:
    output:
        genome = "results/hg38/hg38.fa",
        hisat2_index = expand("results/hg38/hg38.{idx}.ht2", idx=range(1, 9))
    params:
        rsync_path = "rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
    resources:
        mem = "32G",
        time = "4:00:00"
    threads: 8
    shell:
        """
        # Load necessary modules
        ml system
        ml biology
        ml hisat2/2.1.0

        # Ensure the output directory exists
        mkdir -p results/hg38/

        # Download the reference genome using rsync
        rsync -avz --progress {params.rsync_path} results/hg38/hg38.fa.gz

        # Uncompress the downloaded genome
        gunzip -f results/hg38/hg38.fa.gz

        # Build HISAT2 index
        hisat2-build -p {threads} results/hg38/hg38.fa results/hg38/hg38
        """

rule download_hg38_v29_annotation:
    output:
        annotation = "results/hg38/gencode.v29.annotation.gtf"
    params:
        url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"
    resources:
        mem = "4G",
        time = "1:00:00"
    shell:
        """
        # Load necessary modules
        ml load system
        ml load biology
        
        # Ensure the output directory exists
        mkdir -p results/hg38/
        
        # Download the annotation file
        wget {params.url} -O {output.annotation}.gz
    
        # Uncompress the downloaded annotation
        gunzip -f {output.annotation}.gz
        """

rule download_individual_sra_file:
    output:
        fastq_1 = "results/joseph_c_wu_lab/fastq/SRR{sra_id}_1.fastq.gz",
        fastq_2 = "results/joseph_c_wu_lab/fastq/SRR{sra_id}_2.fastq.gz"
    params:
        sra_id = "{sra_id}"
    resources:
        mem = "32G",
        time = "12:00:00"
    threads: 4
    log: "results/logs/download_sra_file_{sra_id}.log"
    run:
        # Strip the '.gz' extension to get the uncompressed filenames
        uncompressed_fastq_1 = output.fastq_1.rstrip('.gz')
        uncompressed_fastq_2 = output.fastq_2.rstrip('.gz')

        shell("""
            # Load required modules
            ml system
            ml biology
            ml sra-tools/3.0.7

            # Create output directory
            mkdir -p $(dirname {output.fastq_1})

            # Download and convert SRA file directly to FASTQ
            fasterq-dump --threads {threads} --progress \
                --outdir $(dirname {output.fastq_1}) \
                --split-files \
                SRR{params.sra_id}

            # Compress the FASTQ files
            gzip -f {uncompressed_fastq_1} {uncompressed_fastq_2}
            """)


rule align_to_hg38:
    input:
        reads=[
            "results/joseph_c_wu_lab/fastq/SRR{sra_id}_1.fastq.gz",
            "results/joseph_c_wu_lab/fastq/SRR{sra_id}_2.fastq.gz"
        ],
        idx=expand("results/hg38/hg38.{idx}.ht2", idx=range(1, 9))
    output:
        "results/joseph_c_wu_lab/bam/SRR{sra_id}.bam"
    log:
        "results/logs/align_sra_file_{sra_id}.log"
    params:
        extra=""  # No additional parameters
    threads: 8
    resources:
        mem="32G",
        time="12:00:00"
    wrapper:
        "v5.0.1/bio/hisat2/align"


rule count_genes:
    input:
        sam="results/joseph_c_wu_lab/bam/SRR{sra_id}.bam",
        annotation="results/hg38/gencode.v29.annotation.gtf"
    output:
        counts="results/joseph_c_wu_lab/counts/SRR{sra_id}.featureCounts",
        summary="results/joseph_c_wu_lab/counts/SRR{sra_id}.featureCounts.summary",
        jcounts="results/joseph_c_wu_lab/counts/SRR{sra_id}.featureCounts.jcounts"
    threads: 4
    params:
        extra="-O --fracOverlap 0.2 -p"
    log: "results/logs/count_genes_{sra_id}.log"
    resources:
      mem = "32G",
      time = "3:00:00"
    wrapper:
        "0.72.0/bio/subread/featurecounts"
        
# List with all SRA IDs
sra_ids = list(range(3618755, 3618794))
rule complete_joseph_c_wu_counts:
    input:
      featureCounts = expand("results/joseph_c_wu_lab/counts/SRR{sra_id}.featureCounts", sra_id=sra_ids), 
      metadata = "resources/joseph_c_wu_lab/SraRunTable.csv", 
      annot = "results/hg38/gencode.v29.annotation.gtf", 
      cardiomyocyte_genes = "results/extra_gene_lists/PanglaoDB_markers_27_Mar_2020.tsv.gz"
    output:
      top_n_genes_w_tpm = "results/joseph_c_wu_lab/plots/top_n_genes_w_tpm.pdf", 
      tpm_distr_per_day = "results/joseph_c_wu_lab/plots/tpm_distr_per_day.pdf",
      expression_variance = "results/joseph_c_wu_lab/plots/expression_variance.pdf",
      threshold_plot = "results/joseph_c_wu_lab/plots/threshold_plot.pdf",
      cardio_genes_that_decrease_expr = "results/joseph_c_wu_lab/plots/cardio_genes_that_decrease_expr.pdf",
      failed_thresh_at_least_once = "results/joseph_c_wu_lab/plots/failed_thresh_at_least_once.pdf",
      gene_list_df = "results/joseph_c_wu_lab/final_gene_list.rds"
    log: "results/logs/complete_joseph_c_wu_counts.log"
    conda:
      "../envs/analyze_crispr_screen.yml"
    resources:
      mem = "32G",
      time = "12:00:00"
    script:
      "../scripts/process_all_datasets/combine_feature_counts.R"
      
      
### ========================================================================

# - the human-kinome-info-brunello files can be found in `https://www.addgene.org/pooled-library/broadgpp-human-kinome/` and was downloaded with `https://media.addgene.org/cms/filer_public/24/f3/24f3f127-7257-463b-90c4-1718c52e48f9/human-kinome-info-brunello.zip`
# - bassik_human_library_composition.xlsx was taken from `https://www.addgene.org/pooled-library/bassik-human-crispr-knockout/` and downloaded from `https://media.addgene.org/cms/filer_public/42/29/4229381e-3796-4922-950f-dc19e3f6afab/bassik_human_library_composition.xlsx`

# Download additional gene lists
rule download_additional_gene_lists:
  output:
    kinase_genes = "results/extra_gene_lists/broadgpp_brunello_Km_pool_14.xlsx",
    membrane_genes = "results/extra_gene_lists/bassik_human_library_composition.xlsx",
    cardiomyocyte_genes = "results/extra_gene_lists/PanglaoDB_markers_27_Mar_2020.tsv.gz"
  shell:
    """
    # Make directory to put files in
    mkdir -p results/extra_gene_lists/
    
    # Download the kinase genes
    wget -O results/extra_gene_lists/human-kinome-info-brunello.zip https://media.addgene.org/cms/filer_public/24/f3/24f3f127-7257-463b-90c4-1718c52e48f9/human-kinome-info-brunello.zip
    unzip -o results/extra_gene_lists/human-kinome-info-brunello.zip -d results/extra_gene_lists/
    
    # Download the membrane genes
    wget -O {output.membrane_genes} https://media.addgene.org/cms/filer_public/42/29/4229381e-3796-4922-950f-dc19e3f6afab/bassik_human_library_composition.xlsx
    
    # Download the cardiomyocyte genes
    # Panglao doesn't like when you try to use wget on the cluster, so we have to change the user-agent
    wget --user-agent="Mozilla/5.0" -O {output.cardiomyocyte_genes} https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz
    """

# Compare the processed datasets by TPM in "mature" cardiomyocytes
rule compare_datasets:
    input:
      sean_m_wu = "results/sean_m_wu_lab/day30_cells.rds",
      joseph_c_wu = "results/joseph_c_wu_lab/final_gene_list.rds",
      liangfu = "results/liangfu_sc/downstream_analysis/cpm_df.rds",
      cardiomyocyte_genes = "results/extra_gene_lists/PanglaoDB_markers_27_Mar_2020.tsv.gz",
      kinase_genes = "results/extra_gene_lists/broadgpp_brunello_Km_pool_14.xlsx",
      membrane_genes = "results/extra_gene_lists/bassik_human_library_composition.xlsx",
    output:
      final_top_13k = "results/final_top_13k.txt",
      plot_joe_sean = "results/compare_datasets/plot_joe_sean.pdf" # There are more plots, but we only need to ensure one is made, so that the "compare_datasets" directory is dynamically created by snakemake
    log: "results/logs/compare_datasets.log"
    conda:
      "../envs/analyze_crispr_screen.yml"
    resources:
      mem = "32G",
      time = "12:00:00"
    script:
      "../scripts/process_all_datasets/compare_datasets.R"
      
      
### =========================================================================
  






