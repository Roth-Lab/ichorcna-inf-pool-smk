from snakemake.utils import min_version, validate

min_version("9.12")

conda: "envs/global.yaml"

validate(config, "schemas/config.schema.yaml")

from utils import ConfigManager

config = ConfigManager(config)

onsuccess:
    config.notification(
        on="success",
        workflow="ichorcna-inf-pool-tf-smk",
        configfile=workflow.configfiles[0],
        imgs=[config.copied_config, config.tfs_plot_file]
    )


onerror:
    config.notification(
        on="error", 
        workflow="ichorcna-inf-pool-tf-smk",
        configfile=workflow.configfiles[0],
    )

rule all:
    input:
        config.pipeline_files


rule build_config_file:
    input:
        workflow.configfiles[0]
    output:
        config.copied_config
    shell:
        "cp {input} {output}"


rule build_cfclone_clone_cn_files:
    input:
        c=config.clone_filter_file,
        i=config.hapclone_results_file,
    output:
        config.cfclone_clone_cn_template,
    params:
        config.get_num_bins_arg
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.cfclone_clone_cn_template),
    shell:
        "(python scripts/build_clone_cn_file.py "
        "-c {input.c} "
        "-i {input.i} "
        "-n {params} "
        "-o {output} ) >{log} 2>&1"


rule build_cfclone_ctdna_file:
    input:
        c=config.cfclone_clone_cn_template,
        d=config.hapclone_data_file,
        r=config.hapclone_results_file,
        s=config.snp_file,
    output:
        config.cfclone_ctdna_template,
    params:
        r=config.read_length,
        p=config.clone_prevalences,
        t=config.get_tumour_content_arg,
        c=config.get_coverage_arg,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.cfclone_ctdna_template),
    benchmark:
        config.get_benchmark_file(config.cfclone_ctdna_template),
    shell:
        "(python scripts/build_ctdna_file.py "
        "-c {input.c} "
        "-d {input.d} "
        "-r {input.r} "
        "-s {input.s} "
        "--coverage {params.c} "
        "--read-length {params.r} "
        "--seed {wildcards.seed} "
        "--tumour-content {params.t} "
        "--clone-prevalence-file {params.p} "
        "-o {output}) >{log} 2>&1"


rule build_ctdna_wig_file:
    input:
        config.cfclone_ctdna_template,
    output:
        config.ctdna_wig_template,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.ctdna_wig_template)
    params:
        config.wig_template_file
    shell:
        "(python scripts/build_ctdna_wig.py "
        "-i {input} "
        "-o {output} "
        "--wig-template-file {params} ) >{log} 2>&1"


rule run_ichorcna:
    input:
        config.ctdna_wig_template,
    output:
        config.replicate_params_txt_file_template,
    log:
        config.get_log_file(config.replicate_params_txt_file_template)
    benchmark:
        config.get_benchmark_file(config.replicate_params_txt_file_template)
    conda:
        "envs/ichorcna.yaml"
    params:
        o=lambda wildcards: config.get_replicate_out_dir(wildcards),
        s=config.ichorcna_settings_file,
    shell:
        """
        myichorcna run -o {params.o} -s {params.s} -c {input} >{log} 2>&1
        """


rule build_replicate_summary_file:
    input:
        config.replicate_params_txt_file_template
    output:
        config.replicate_summary_file_template
    log:
        config.get_log_file(config.replicate_summary_file_template)
    params:
        cp=config.clone_prevalences,
        n=config.get_num_bins_arg,
        cov=config.get_coverage_arg,
        tc=config.get_tumour_content_arg,
    shell:
        "(python scripts/write_summary_file.py "
        "-i {input} "
        "-o {output} "
        "--coverage {params.cov} "
        "--tumour-content {params.tc} "
        "--clone-prevalences {params.cp} "
        "--num-bins {params.n} ) >{log} 2>&1"


rule merge_summary_files:
    input:
        config.get_summary_files
    output:
        config.summary_file_tfs
    params:
        config.out_dir
    log:
        config.get_log_file(config.summary_file_tfs)
    conda:
        "envs/python.yaml"
    shell:
        "(python scripts/merge_tables.py -i {input} -o {output}) >{log} 2>&1"


rule plot_tfs_summary:
    input:
        config.summary_file_tfs
    output:
        config.tfs_plot_file
    log:
        config.get_log_file(config.tfs_plot_file)
    conda:
        "envs/plot.yaml"
    shell:
        "(python scripts/plot_summary_tfs.py -i {input} -o {output}) >{log} 2>&1"