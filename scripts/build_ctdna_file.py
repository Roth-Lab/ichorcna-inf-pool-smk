from statsmodels.nonparametric.smoothers_lowess import lowess

import h5py
import math
import numpy as np
import pandas as pd
import pysam
import scipy
import statsmodels.formula.api as smf


def main(args):
    rng = np.random.RandomState(args.seed)

    cell_profiles, cell_to_clone, data = load_data(
        args.clone_cn_file, args.hapclone_data_file, args.hapclone_results_file
    )

    ploidy = cell_profiles.sum(axis=-1).mean(axis=1)

    target_normal_reads, target_tumour_reads = compute_target_reads(
        args.coverage,
        data,
        ploidy,
        args.read_length,
        args.tumour_content,
    )

    normal_reads = generate_normal_reads(
        data,
        target_normal_reads,
        rng,
    )

    cell_coverage = generate_cell_coverage(
        cell_to_clone,
        data,
        rng,
        target_tumour_reads,
        clone_prevalence_file=args.clone_prevalence_file,
    )

    cell_reads = generate_cell_reads(
        cell_coverage,
        data,
        rng,
    )

    df = data.bin_df

    df["reads"] = normal_reads + cell_reads.sum(axis=0)

    snp_density = compute_snp_density(data, args.snp_file)

    df["a"], df["b"] = generate_baf(
        cell_profiles,
        cell_reads,
        data,
        normal_reads,
        args.read_length,
        snp_density,
        rng,
    )

    df["gc"] = data.gc

    df["map"] = data.map

    df["valid"] = data.valid

    df = gc_correction(df)

    df = df.drop(columns="rdr")

    df = df.rename(columns={"beg": "start", "rdr_cor": "rdr"})

    df = df[df["valid"]]
    
    # df_clone_cn = pd.read_csv(args.clone_cn_file, sep="\t")
    
    # df_bins = df_clone_cn.drop_duplicates(subset=["chrom", "start", "end"])
    
    # df = df.merge(df_bins, on=["chrom", "start", "end"], how="inner")

    df.to_csv(args.out_file, index=False, sep="\t")


def load_data(clone_cn_file, data_file, results_file):
    clones = pd.read_csv(clone_cn_file, sep="\t")["clone"].astype(int).unique()

    data = DataSet.from_file(data_file)

    results_df = pd.read_csv(results_file, sep="\t")

    cell_to_clone = results_df[["cluster_id", "cell_id"]].drop_duplicates().set_index("cell_id")["cluster_id"].to_dict()

    bin_df = results_df[["chrom", "beg", "end"]].drop_duplicates()

    bin_df["bin_id"] = bin_df.apply(lambda row: "{chrom}:{beg}:{end}".format(**row.to_dict()), axis=1)

    results_df = pd.merge(bin_df, results_df, on=["chrom", "beg", "end"])

    data.filter_bins(results_df["bin_id"].unique())

    cells = results_df[results_df["cluster_id"].isin(clones)]["cell_id"].unique()

    data.filter_cells(cells)

    cell_cn_a = results_df.pivot(index="cell_id", columns="bin_id", values="cn_A_cell").loc[data.cells, data.bins]

    cell_cn_b = results_df.pivot(index="cell_id", columns="bin_id", values="cn_B_cell").loc[data.cells, data.bins]

    cell_profiles = np.stack([cell_cn_a.values, cell_cn_b.values], axis=-1)

    return cell_profiles, cell_to_clone, data


def compute_snp_density(data, snp_file):
    snp_reader = pysam.VariantFile(snp_file, "r")

    df = data.bin_df

    snp_density = np.zeros(data.num_bins, dtype=float)

    for i, (_, row) in enumerate(df.iterrows()):
        count = 0

        for record in snp_reader.fetch(row["chrom"], row["beg"], row["end"]):
            # Skip indels
            if len(record.alleles[0]) != 1 or len(record.alleles[1]) != 1:
                continue

            sample_id = list(record.samples.keys())[0]

            sample_info = record.samples[sample_id]

            # Count het SNPs
            if sample_info["GT"] in set([(0, 1), (1, 0)]):
                count += 1

        bin_len = row["end"] - row["beg"]

        snp_density[i] = count / bin_len

    return snp_density


def compute_target_reads(coverage, data, ploidy, read_length, tumour_content):
    target_reads = (coverage * data.genome_size) / read_length

    m = ploidy.mean()

    target_normal_reads = math.ceil((2 / (2 + m)) * (1 - tumour_content) * target_reads)

    target_tumour_reads = math.ceil((m / (2 + m)) * tumour_content * target_reads)

    return target_normal_reads, target_tumour_reads


def generate_cell_coverage(cell_to_clone, data, rng, target_reads, clone_prevalence_file=None):
    """
    Generate the total number of reads coming from each cell
    """
    if clone_prevalence_file is None:
        p = rng.dirichlet(np.ones(data.num_cells))

        # Sample total number DNA fragments from cells
        cell_coverage = rng.multinomial(target_reads, p)

    else:
        df = pd.read_csv(clone_prevalence_file, converters={"clone_id": int}, sep="\t")

        # Allows the output of `cfclone write-prevalence-stats` to be passed without modification
        df = df.rename(columns={"mean_prevalence": "prevalence"})

        clone_prev = df.set_index("clone_id")["prevalence"].to_dict()

        cell_coverage = np.zeros(data.num_cells)

        for c in clone_prev:
            # Number of reads for the clone
            clone_reads = np.ceil(target_reads * clone_prev[c])

            clone_cells = [k for k, v in cell_to_clone.items() if v == c]

            idxs = [data.cells.index(x) for x in clone_cells]

            # Randomly allocate reads to cells belonging to clone
            cell_props = rng.dirichlet(np.ones(len(clone_cells)))

            cell_coverage[idxs] = rng.multinomial(clone_reads, cell_props)

    return cell_coverage


def generate_cell_reads(cell_coverage, data, rng):
    """
    Generate reads coverage profile for each cell
    """
    cell_reads = np.zeros((data.num_cells, data.num_bins), dtype=int)

    # Sample coverage by bin for each cell
    for i in range(data.num_cells):
        p = rng.dirichlet(data.reads[i].sum(axis=-1) + 1e-6)

        cell_reads[i] = rng.multinomial(cell_coverage[i], p)

    return cell_reads


def generate_normal_reads(data, target_reads, rng):
    p = rng.dirichlet(data.normal_reads + 1e-6)
    return rng.multinomial(target_reads, p)


def generate_baf(
    cell_profiles,
    cell_reads,
    data,
    normal_reads,
    read_length,
    snp_density,
    rng,
):
    tumour_a, tumour_b = generate_tumour_baf(
        cell_profiles,
        cell_reads,
        data,
        read_length,
        snp_density,
        rng,
    )

    normal_a, normal_b = generate_normal_baf(
        normal_reads,
        read_length,
        snp_density,
        rng,
    )

    a = tumour_a + normal_a

    b = tumour_b + normal_b

    return a, b


def generate_normal_baf(reads, read_length, snp_density, rng):
    e_snp = reads * snp_density * read_length

    d = rng.poisson(e_snp)

    normal_a = rng.binomial(d, 0.5)

    normal_b = d - normal_a

    return normal_a, normal_b


def generate_tumour_baf(cell_profiles, cell_reads, data, read_length, snp_density, rng):
    tumour_a = np.zeros(data.num_bins, dtype=int)

    tumour_b = np.zeros(data.num_bins, dtype=int)

    for i in range(data.num_cells):
        for j in range(data.num_bins):
            # Expected number of SNPs covered by a read in a bin
            e_snp = snp_density[j] * cell_reads[i, j] * read_length

            d = rng.poisson(e_snp)

            m = cell_profiles[i, j, 0] / cell_profiles[i, j].sum()

            a = rng.binomial(d, m)

            b = d - a

            tumour_a[j] += a

            tumour_b[j] += b

    return tumour_a, tumour_b


def gc_correction(df):
    df["rdr"] = df["reads"]

    # Filtering and sorting
    df_regression = df[df["valid"]].copy()

    df_regression.sort_values(by="gc", inplace=True)

    df_regression = modal_quantile_regression(df_regression, lowess_frac=0.2)

    df["gc_correction"] = np.nan

    df["rdr_cor"] = np.nan

    df.loc[df_regression.index, "gc_correction"] = df_regression["modal_curve"]

    df["rdr_cor"] = df["rdr"] / df["gc_correction"]

    return df


def modal_quantile_regression(df_regression, lowess_frac=0.2, degree=2, knots=(0.38,)):
    """
    Fits a B-spline polynomial curve through the "modal" quantile of the data:
    * Runs quantile regression to fit a B-spline curve for each percentile 10-90
    * Estimates the modal quantile as the quantile where difference in AUC is minimized
    * Uses the curve fit to this modal quantile for normalization

    Parameters:
        df_regression: pandas.DataFrame with at least columns [chr, start, end, rdr, gc]
        lowess_frac: float, fraction of data used to estimate each y-value in Lowess smoothing of AUC curve
        degree: int, degree of polynomial to fit to each section of the B-spline curve
        knots: list of floats, GC values where B-spline polynomial is allowed to change

    Returns:
        pandas.DataFrame with additional columns
            modal_curve: modal curve's predicted # rdr for GC value in this row
            modal_quantile: quantile selected as the mode (should be the same for all bins)
            modal_corrected: corrected read count (i.e., rdr / modal_curve)
    """

    q_range = range(10, 91, 1)
    quantiles = np.array(q_range) / 100
    quantile_names = [str(x) for x in q_range]

    # need at least 3 values to compute the quantiles
    if len(df_regression) < 10 or sum(df_regression["rdr"]) < 100:
        df_regression["modal_quantile"] = None
        df_regression["modal_curve"] = None
        df_regression["modal_corrected"] = None
        return df_regression

    poly_quantile_model = smf.quantreg(
        f"rdr ~ bs(gc, degree={degree}, knots={knots}, include_intercept = True)",
        data=df_regression,
    )

    poly_quantile_fit = [poly_quantile_model.fit(q=q, max_iter=10000) for q in quantiles]

    poly_quantile_predict = [poly_quantile_fit[i].predict(df_regression) for i in range(len(quantiles))]

    poly_quantile_params = pd.DataFrame()

    for i in range(len(quantiles)):
        df_regression[quantile_names[i]] = poly_quantile_predict[i]
        poly_quantile_params[quantile_names[i]] = poly_quantile_fit[i].params

    # integration and mode selection
    gc_min = df_regression["gc"].quantile(q=0.10)
    gc_max = df_regression["gc"].quantile(q=0.90)

    true_min = df_regression["gc"].min()
    true_max = df_regression["gc"].max()

    poly_quantile_integration = np.zeros(len(quantiles) + 1)

    # form (k+1)-regular knot vector
    repeats = degree + 1

    my_t = np.r_[[true_min] * repeats, knots, [true_max] * repeats]

    for i in range(len(quantiles)):
        # compose params into piecewise polynomial
        params = poly_quantile_params[quantile_names[i]].to_numpy()
        pp = scipy.interpolate.PPoly.from_spline((my_t, params[1:] + params[0], degree))

        # compute integral
        poly_quantile_integration[i + 1] = pp.integrate(gc_min, gc_max)

    # find the modal quantile
    distances = poly_quantile_integration[1:] - poly_quantile_integration[:-1]

    df_dist = pd.DataFrame(
        {
            "quantiles": quantiles,
            "quantile_names": quantile_names,
            "distances": distances,
        }
    )
    dist_max = df_dist["distances"].quantile(q=0.95)
    df_dist_filter = df_dist[df_dist["distances"] < dist_max].copy()
    df_dist_filter["lowess"] = lowess(
        df_dist_filter["distances"],
        df_dist_filter["quantiles"],
        frac=lowess_frac,
        return_sorted=False,
    )

    modal_quantile = df_dist_filter.set_index("quantile_names")["lowess"].idxmin()

    # add values to table
    df_regression["modal_quantile"] = modal_quantile

    df_regression["modal_curve"] = df_regression[modal_quantile]

    return df_regression


class DataSet(object):
    @staticmethod
    def from_file(file_name):
        with h5py.File(file_name, "r") as fh:
            baf = fh["baf"][()]

            reads = fh["reads"][()]

            normal_reads = fh["normal"][()]

            bins = [x.decode() for x in fh["bins"]]

            cells = [x.decode() for x in fh["cells"]]

            gc = fh["gc"][()]

            mappability = fh["map"][()]

            regions = fh["region"][()]

            valid = fh["valid"][()]

        return DataSet(bins, cells, baf, reads, normal_reads, gc, mappability, regions, valid)

    def __init__(self, bins, cells, baf, reads, normal_reads, gc, mappability, regions, valid):
        self.bins = bins

        self.cells = cells

        self.baf = baf

        self.reads = reads

        self.normal_reads = normal_reads

        self.gc = gc

        self.map = mappability

        self.regions = regions

        self.valid = valid

    @property
    def bin_df(self):
        df = pd.DataFrame([x.split(":") for x in self.bins], columns=["chrom", "beg", "end"])

        df["chrom"] = df["chrom"].astype(str)

        df["beg"] = df["beg"].astype(int)

        df["end"] = df["end"].astype(int)

        return df

    @property
    def genome_size(self):
        df = self.bin_df

        df["len"] = df["end"] - df["beg"]

        return df["len"].sum()

    @property
    def num_bins(self):
        return self.baf.shape[1]

    @property
    def num_blocks(self):
        return self.baf.shape[2]

    @property
    def num_cells(self):
        return self.baf.shape[0]

    def filter_bins(self, bins):
        idxs = [self.bins.index(x) for x in bins]

        self.bins = [self.bins[i] for i in idxs]

        self.baf = self.baf[:, idxs]

        self.reads = self.reads[:, idxs]

        self.normal_reads = self.normal_reads[idxs]

        self.gc = self.gc[idxs]

        self.map = self.map[idxs]

        self.regions = self.regions[idxs]

        self.valid = self.valid[idxs]

    def filter_cells(self, cells):
        idxs = [self.cells.index(x) for x in cells]

        self.cells = [self.cells[i] for i in idxs]

        self.baf = self.baf[idxs]

        self.reads = self.reads[idxs]

    def remove_centromeres(self):
        bins = [self.bins[i] for i, r in enumerate(self.regions) if r != -1]

        self.filter_bins(bins)


if __name__ == "__main__":
    import argparse
    
    default0 = "/home/matteo/projects/cfdna/wfs/results/cfclone-inf-pool-cal-smk/TFRI004/out_dir/input/clone_cn.tsv.gz"
    
    default1 = "/home/matteo/projects/cfdna/wfs/data/cfclone-inf-pool-cal-smk/TFRI004/hapclone/data.h5"
    
    default2 = "/home/matteo/projects/cfdna/wfs/data/cfclone-inf-pool-cal-smk/TFRI004/hapclone_refit/merged_results.tsv.gz"
    
    default3 = "/home/matteo/projects/cfdna/wfs/data/cfclone-inf-pool-cal-smk/TFRI004/hapclone_refit/rephased_snps.bcf"
    
    default4 = "/home/matteo/projects/cfdna/wfs/results/cfclone-inf-pool-cal-smk/TFRI004/out_dir/input/ctdna/coverage_0/tc_0/cp_0/replicate_0.tsv.gz"
    
    default5 = "/home/matteo/projects/cfdna/wfs/configs/cfclone-inf-pool-cal-smk/TFRI004/clone-prevs/clone_prevs_00.tsv"
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--clone-cn-file", default=default0)

    parser.add_argument("-d", "--hapclone-data-file", default=default1)

    parser.add_argument("-r", "--hapclone-results-file", default=default2)

    parser.add_argument("-s", "--snp-file", default=default3)

    parser.add_argument("-o", "--out-file", default=default4)

    parser.add_argument("--coverage", default=1, type=float)

    parser.add_argument("--read-length", default=150, type=int)

    parser.add_argument("--seed", default=None, type=int)

    parser.add_argument("--tumour-content", default=0.1, type=float)
    
    parser.add_argument("--clone-prevalence-file", default=default5)

    cli_args = parser.parse_args()

    main(cli_args)