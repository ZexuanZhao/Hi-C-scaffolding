rule copy_assembly:
    input:
        config["asm"]
    output:
        os.path.join(out_dir,"assembly","assembly.fasta")
    shell:
        """
        cp {input} {output}
        """

rule index_assembly:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"assembly","assembly.fasta")
    output:
        os.path.join(out_dir,"assembly","assembly.fasta.bwt"),
        os.path.join(out_dir,"assembly","assembly.fasta.fai")
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        """

rule hic_mapping:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        bwt = os.path.join(out_dir,"assembly","assembly.fasta.bwt"),
        assembly = os.path.join(out_dir,"assembly","assembly.fasta"),
        R1 = config["hiC_read1"],
        R2 = config["hiC_read2"]
    output:
        R1_mapped = os.path.join(out_dir,"hic_mapping", "R1toasm.bam"),
        R2_mapped = os.path.join(out_dir,"hic_mapping", "R2toasm.bam")
    threads:
        threads
    shell:
        """
        bwa mem \
            -t {threads}\
            {input.assembly} \
            {input.R1}| \
            samtools view -@ {threads} -Sb - \
            > {output.R1_mapped}
        bwa mem \
            -t {threads}\
            {input.assembly} \
            {input.R2}| \
            samtools view -@ {threads} -Sb - \
            > {output.R2_mapped}
        """

rule filter5end:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        R1_mapped = os.path.join(out_dir,"hic_mapping", "R1toasm.bam"),
        R2_mapped = os.path.join(out_dir,"hic_mapping", "R2toasm.bam")
    output:
        R1_5endFiltered = os.path.join(out_dir,"hic_mapping", "R1toasm.5endFiltered.bam"),
        R2_5endFiltered = os.path.join(out_dir,"hic_mapping", "R2toasm.5endFiltered.bam")
    threads:
        threads
    shell:
        """
        samtools view -h -@ {threads} {input.R1_mapped} | \
            filter_five_end.pl | \
            samtools view -Sb -@ {threads} - \
            > {output.R1_5endFiltered}
        samtools view -h -@ {threads} {input.R2_mapped} | \
            filter_five_end.pl | \
            samtools view -Sb -@ {threads} - \
            > {output.R2_5endFiltered}
        """

rule conbine_and_filter_bams:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        R1_5endFiltered = os.path.join(out_dir,"hic_mapping", "R1toasm.5endFiltered.bam"),
        R2_5endFiltered = os.path.join(out_dir,"hic_mapping", "R2toasm.5endFiltered.bam"),
        fai = os.path.join(out_dir,"assembly","assembly.fasta.fai")
    output:
         os.path.join(out_dir,"hic_mapping", "hic2asm.combined.filtered.bam")
    threads:
        threads
    params:
        mapq_filter=config["mapq_filter"]
    shell:
        """
        two_read_bam_combiner.pl {input.R1_5endFiltered} {input.R2_5endFiltered} samtools {params.mapq_filter} | \
            samtools view -bS -@ {threads} -t {input.fai} - | \
            samtools sort -@ {threads} -o {output}
        """

rule mark_duplicate:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "hic2asm.combined.filtered.bam")
    output:
        metric = os.path.join(out_dir,"qc", "markDuplicate.metric.txt"),
        bam = os.path.join(out_dir,"hic_mapping", "hic2asm.combined.filtered.purged.bam")
    params:
        mem = config["mem"]
    shell:
        """
        picard MarkDuplicates \
            -Xmx{params.mem} -XX:-UseGCOverheadLimit \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metric} \
            ASSUME_SORTED=TRUE \
            VALIDATION_STRINGENCY=LENIENT\
            REMOVE_DUPLICATES=TRUE
        """

rule sort_by_name_bam:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "hic2asm.combined.filtered.purged.bam")
    output:
        os.path.join(out_dir,"hic_mapping", "hic2asm.combined.filtered.purged.sorted.bam")
    threads:
        threads
    shell:
        """
        samtools sort -@ {threads} -o {output} -n {input}
        """
