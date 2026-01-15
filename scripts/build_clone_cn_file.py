import h5py
import numpy as np
import pandas as pd


def main(args):
    df = pd.read_csv(args.in_file, sep="\t")

    df = df[["chrom", "beg", "end", "cluster_id", "cn_A", "cn_B"]].drop_duplicates()

    df = df.rename(
        columns={"cn_A": "cn_a", "cn_B": "cn_b", "beg": "start", "cluster_id": "clone"}
    )

    df["clone"] = df["clone"].astype(int)

    if args.clone_filter_file is not None:
        
        clone_df = pd.read_csv(args.clone_filter_file, sep="\t")

        clone_df = clone_df[clone_df["keep"]]

        clones = clone_df["clone_id"].astype(int).unique()

        df = df[df["clone"].isin(clones)]

    bin_df = df[["chrom", "start", "end"]].drop_duplicates()
    
    if args.num_bins != 'all':

        if 0 < args.num_bins < bin_df.shape[0]:
            
            idxs = np.random.choice(bin_df.shape[0], replace=False, size=args.num_bins)

            idxs = sorted(idxs)

            bin_df = bin_df.iloc[idxs]

            df = pd.merge(df, bin_df, on=["chrom", "start", "end"])

    df.to_csv(args.out_file, index=False, sep="\t")


def convert_cn_state_to_profiles(A, cell_states):
    state_map = get_state_map(A)

    cell_profiles = np.zeros((cell_states.shape[0], cell_states.shape[1], 2), dtype=int)

    for i in range(cell_states.shape[0]):
        for j in range(cell_states.shape[1]):
            cell_profiles[i, j] = state_map[cell_states[i, j]]

    return cell_profiles


def get_bin_df(bins):
    df = pd.DataFrame([x.split(":") for x in bins], columns=["chrom", "start", "end"])

    df["chrom"] = df["chrom"].astype(str)

    df["start"] = df["start"].astype(int)

    df["end"] = df["end"].astype(int)

    return df


def get_state_map(A):
    s_map = np.empty((A * A, 2), dtype=np.int32)

    s = 0

    for u in range(A):
        for v in range(A):
            s_map[s, 0] = u

            s_map[s, 1] = v

            s += 1

    return s_map


if __name__ == "__main__":
    import argparse
    
    def num_bins_parser(arg):
        if arg == "all":
            return arg
        else:
            return int(arg)

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-c", "--clone-filter-file", default=None)

    parser.add_argument("-n", "--num-bins", type=num_bins_parser)

    cli_args = parser.parse_args()

    main(cli_args)
