
# Data: scRNA-seq
# Sample: iPSCs derived CM (2D differentiation)-Day30
# Lab work by: Liangfu

# USE STARSOLO TO ALIGN READS TO GENOME ======================================

SAMPLES = ["P1_24s001312-1-1_Xie_lane124s001312", "P2_24s001313-1-1_Xie_lane124s001313"]

rule star_index:
    input:
        fasta="results/hg38/hg38.fa",
        gtf="results/hg38/gencode.v29.annotation.gtf"
    output:
        directory("results/liangfu_sc/star_index")
    threads: 8
    resources:
      mem = "64G",
      time = "4:00:00"
    shell:
        """
        ml load biology
        ml load star/2.7.10b
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf}
        """

rule run_starsolo:
    input:
        whitelist="resources/liangfu_dataset/3M-february-2018.txt.gz",
        r2="resources/liangfu_dataset/2024-11-19-AAG52FLM5/AAG52FLM5_CM-{sample}_2_sequence.txt.gz",    # cDNA
        r1="resources/liangfu_dataset/2024-11-19-AAG52FLM5/AAG52FLM5_CM-{sample}_1_sequence.txt.gz",    # Barcode
        index="results/liangfu_sc/star_index"
    output:
        matrix="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/matrix.mtx",
        features="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/features.tsv",
        barcodes="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/barcodes.tsv",
        filter_matrix="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/matrix.mtx",
        filtered_features="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/features.tsv",
        filtered_barcodes="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/barcodes.tsv"
    resources:
      mem = "64G",
      time = "8:00:00"
    threads: 12
    shell:
       """
        ml load biology
        ml load star/2.7.10b
        mkdir -p $(dirname {output.matrix})
        # Decompress whitelist into a temporary file
        tmp_whitelist=$(mktemp)
        gunzip -c {input.whitelist} > $tmp_whitelist
        STAR --genomeDir {input.index} \
             --readFilesIn {input.r2} {input.r1} \
             --readFilesCommand zcat \
             --runThreadN {threads} \
             --soloType CB_UMI_Simple \
             --soloCBwhitelist {input.whitelist} \
             --soloUMIlen 12 \
             --outFileNamePrefix $(dirname $(dirname $(dirname $(dirname {output.matrix}))))/
        """
        
# To get the Barcode rank plot colored like cellranger, i used this github thread: https://github.com/alexdobin/STAR/issues/1158
# The following is run: 
  # curl -O https://raw.githubusercontent.com/alexdobin/STAR/master/extras/scripts/soloUMIperCell.awk
  # chmod +x soloUMIperCell.awk
rule generate_umi_per_cell:
    input:
        raw_matrix="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/matrix.mtx",
        raw_barcodes="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/barcodes.tsv",
        filtered_barcodes="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/barcodes.tsv"
    output:
        # "UMIperCell.txt" will have 2 columns: totalUMIs, isFiltered(0/1)
        "results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/UMIperCell.txt"
    resources:
        mem = "8G",
        time = "1:00:00"
    shell:
        """
        wget -O workflow/scripts/liangfu/soloUMIperCell.awk https://raw.githubusercontent.com/alexdobin/STAR/master/extras/scripts/soloUMIperCell.awk 
        chmod +x workflow/scripts/liangfu/soloUMIperCell.awk
        
        # Make sure AWK is available or use a conda env if needed
        # Use the AWK script from your scripts folder:
        awk -f workflow/scripts/liangfu/soloUMIperCell.awk \
            {input.raw_matrix} \
            {input.raw_barcodes} \
            {input.filtered_barcodes} \
        | sort -k1,1rn \
        > {output}
        """

rule compress_raw_outputs:
    input:
        mat="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/matrix.mtx",
        bar="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/barcodes.tsv",
        feat="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/features.tsv"
    output:
        mat_gz="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/matrix.mtx.gz",
        bar_gz="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/barcodes.tsv.gz",
        feat_gz="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/features.tsv.gz"
    shell:
        """
        gzip -c {input.mat}  > {output.mat_gz}
        gzip -c {input.bar}  > {output.bar_gz}
        gzip -c {input.feat} > {output.feat_gz}
        """
        
rule compress_filtered_outputs:
    input:
        mat="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/matrix.mtx",
        bar="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/barcodes.tsv",
        feat="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/features.tsv",
    output:
        mat_gz="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/matrix.mtx.gz",
        bar_gz="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/barcodes.tsv.gz",
        feat_gz="results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/features.tsv.gz",
    shell:
        """
        gzip -c {input.mat}  > {output.mat_gz}
        gzip -c {input.bar}  > {output.bar_gz}
        gzip -c {input.feat} > {output.feat_gz}
        """

rule downstream_analysis:
    input:
        raw_matrices = expand("results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/raw/matrix.mtx.gz", sample=SAMPLES),
        filt_matrices = expand("results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/filtered/matrix.mtx.gz", sample=SAMPLES),
        umi_files = expand("results/liangfu_sc/run_starsolo/{sample}/Solo.out/Gene/UMIperCell.txt", sample=SAMPLES)
    output:
        combined_seurat = "results/liangfu_sc/downstream_analysis/combined_seurat.rds",
        cpm_df = "results/liangfu_sc/downstream_analysis/cpm_df.rds"
    log: "results/liangfu_sc/logs/downstream_analysis.log"
    conda:
        "../envs/sean_m_wu.yml"
    resources:
      mem = "64G",
      time = "6:00:00"
    script:
        "../scripts/liangfu/processing_raw_counts_data.R"
