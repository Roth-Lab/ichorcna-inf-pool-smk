import pathlib

from itertools import product

from snakemake.shell import shell

class ConfigManager(object):
    def __init__(self, config):
        self.config = config
        
    # WILDCARDS 
   
    @property
    def coverage(self):
        return self.config["coverage"]
    
    @property
    def coverage_ids(self):
        return list(range(len(self.coverage)))

    @property
    def tumour_content(self):
        return self.config["tumour_content"]

    @property
    def tumour_content_ids(self):
        return list(range(len(self.tumour_content)))
    
    @property
    def num_replicates(self):
        return self.config.get("num_replicates", 1)
    
    # DATA GENERATION SETTINGS 
    
    @property
    def clone_prevalences(self):
        return self.config['clone_prevalences']
    
    @property
    def num_bins(self):
        return self.config['num_bins']
    
    @property
    def read_length(self):
        return self.config.get("read_length", 150)

    # INPUT FILES FOR DATA GEN
    
    @property
    def clone_filter_file(self):
        return pathlib.Path(self.config["clone_filter_file"]).resolve()

    @property
    def hapclone_data_file(self):
        return pathlib.Path(self.config["hapclone_data_file"]).resolve()

    @property
    def hapclone_results_file(self):
        return pathlib.Path(self.config["hapclone_results_file"]).resolve()

    @property
    def snp_file(self):
        return pathlib.Path(self.config["snp_file"]).resolve()
    
    # ICHORCHA MODEL SETTINGS 
    
    @property
    def wig_template_file(self):
        return pathlib.Path(self.config["wig_template_file"]).resolve()
    
    @property
    def ichorcna_settings_file(self):
        return pathlib.Path(self.config["ichorcna_settings_file"]).resolve()

    # OUTPUT DIRECTORIES FILES 
    
    @property
    def out_dir(self):
        return pathlib.Path(self.config["out_dir"]).resolve()

    @property
    def copied_config(self):
        return self.out_dir.joinpath("config.yaml")

    @property
    def pipeline_dir(self):
        return pathlib.Path(self.config["pipeline_dir"]).resolve()

    @property
    def log_dir(self):
        return self.pipeline_dir.joinpath("log")
    
    @property
    def benchmark_dir(self):
        return self.pipeline_dir.joinpath("benchmark")

    @property
    def tmp_dir(self):
        return self.pipeline_dir.joinpath("tmp")

    # DATA GENERATION OUTPUT TEMPLATES
    
    @property
    def cfclone_clone_cn_dir(self):
        return self.out_dir.joinpath("input", "clone_cn")
    
    @property
    def cfclone_clone_cn_template(self):
        return self.cfclone_clone_cn_dir.joinpath("clone_cn.tsv.gz")

    @property
    def cfclone_ctdna_dir(self):
        return self.out_dir.joinpath(
            "input",
            "ctdna",
            "coverage_{coverage_id}",
            "tc_{tumour_content_id}",
        )
    
    @property
    def cfclone_ctdna_template(self):
        return self.cfclone_ctdna_dir.joinpath("replicate_{seed}.tsv.gz")
    
    @property
    def ctdna_wig_template(self):
        return self.cfclone_ctdna_dir.joinpath("replicate_{seed}.wig")
        
    # ICHORCNA OUTPUT TEMPLATES
    
    @property
    def replicate_dir(self):
        return self.tmp_dir.joinpath(
            "coverage_{coverage_id}",
            "tc_{tumour_content_id}", 
            "replicate_{seed}",
        )

    @property
    def replicate_out_dir(self):
        return self.replicate_dir.joinpath("results")

    @property
    def replicate_pipeline_dir(self):
        return self.replicate_dir.joinpath("tmp")
    
    @property
    def replicate_params_txt_file_template(self):
        return self.replicate_out_dir.joinpath("DummySampleID.params.txt")
    
    @property
    def replicate_summary_file_template(self):
        return self.replicate_out_dir.joinpath("summary_tfs.tsv")
    
    # FINAL OUTPUT FILES
    
    @property
    def summary_file_tfs(self):
        return self.out_dir.joinpath("summary_tfs.tsv")
    
    @property
    def tfs_plot_file(self):
        return self.out_dir.joinpath("tfs.png")
    
    @property
    def pipeline_files(self):
        
        files = []
        
        files.append(self.copied_config)
        
        files.append(self.tfs_plot_file)
        
        return files

    # HELPER FUNCTIONS FOR RULES 

    def get_coverage_arg(self, wildcards):
        return float(self.coverage[int(wildcards.coverage_id)])

    def get_tumour_content_arg(self, wildcards):
        return float(self.tumour_content[int(wildcards.tumour_content_id)])
    
    def get_clone_prevalences_file_arg(self, wildcards):
        return str(self.clone_prevalences[int(wildcards.clone_prevalences_id)])
    
    @property
    def get_num_bins_arg(self):
        num_bins = self.num_bins
        if num_bins == "all":
            return num_bins
        else:
            return int(num_bins)
    
    def get_replicate_out_dir(self, wildcards):
        return str(
            self.replicate_out_dir
        ).format(
            coverage_id=wildcards.coverage_id,
            tumour_content_id=wildcards.tumour_content_id,
            seed=wildcards.seed,
        )
    
    @property
    def combos(self):
        return [self.coverage_ids, self.tumour_content_ids, range(self.num_replicates)]
    
    @property
    def get_summary_files(self):
        
        files = []
        
        for cov, tc, seed in product(*self.combos):
            
            files.append(str(self.replicate_summary_file_template).format(
                coverage_id=cov,
                tumour_content_id=tc,
                seed=seed
                )
            )
                        
        return files 
    
    
    # HELPERS FOR LOG AND BENCHMARK FILES
    
    def get_log_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.log_dir.joinpath(parent, rel_path)
    
    def get_benchmark_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.benchmark_dir.joinpath(parent, rel_path)

    def _get_relative_path(self, template):
        try:
            rel_path = template.relative_to(self.pipeline_dir)
            parent = "working"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path
    
    # HELPERS FOR EMAIL NOTIFICATIONS
    
    @property
    def email(self) -> str:
        return self.config.get("email", "lepurmatteo@gmail.com")

    def notification(self, config_file: str, on: str) -> None:

        msg_template = "configfile: {config}"

        subj_template = "ichorcna-inf-pool-power-calc-smk: {on}"
        
        cmd_template = "echo {msg} | mail -s {sub} {email}"
        
        msg = msg_template.format(config=config_file)
        
        subj = subj_template.format(on=on)

        cmd = cmd_template.format(msg=msg, sub=subj, email=self.email)

        shell(cmd)
