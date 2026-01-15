import pandas as pd


def main(args):
    
    df = []

    for file_name in args.in_files:
        
        df.append(pd.read_csv(file_name, sep="\t"))

    df = pd.concat(df)

    df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-files", nargs="+", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
