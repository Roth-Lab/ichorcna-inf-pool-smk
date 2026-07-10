import pandas as pd

from pandas.io.common import StringIO

from pathlib import Path


def load_tfs(params_txt_file: str | Path) -> pd.DataFrame:
    """
    Args:
        file_path (str | Path): Path to the text file containing ichorCNA output.

    Returns:
        pd.DataFrame: DataFrame with columns such as:
            [
                'init',
                'n_est',
                'phi_est',
                'BIC',
                'Frac_genome_subclonal',
                'Frac_CNA_subclonal',
                'loglik'
            ]
    """
    
    # LOAD RESULLTS FROM TXT FILE 
    
    file_path = Path(params_txt_file)
    
    with file_path.open("r") as f:
        
        lines = f.readlines()
        
    # ICHORCNA FAILED TO RUN 
        
    if lines[0].startswith("Failed"):
        
        return pd.DataFrame({'status': ['Failed']})
    
    # FIND START OF TABLE WITH PARAMETERS FOR EACH EM INITIALISATION 
    
    start_idx = None
    
    for i, line in enumerate(lines):
        
        if line.strip().startswith("init"):
            
            start_idx = i
            
            break
    
    if start_idx is None:
        
        raise ValueError(f"No 'init' table header found in {file_path.name}")
    
    # CREATE DATAFRAME WITH RESULTS 
    
    table_text = "".join(lines[start_idx:])
    
    data = StringIO(table_text)
    
    df = pd.read_csv(data, sep="\t")
    
    df['status'] = 'Completed'

    return df


def main(args):
    
    df = load_tfs(params_txt_file=args.in_file)
    
    df.to_csv(args.out_file, sep="\t", index=False)
    
    df.insert(0, 'coverage', args.coverage)
    
    df.insert(1, 'tumour_content', args.tumour_content)
    
    df.insert(2, 'clone_prevalence', Path(args.clone_prevalences).stem)
    
    df.insert(3, 'num_bins', args.num_bins)
    
    df.insert(4, 'ichorcna_settings_file', args.ichorcna_settings_file)
    
    df.insert(5, 'replicate', args.replicate)
    
    df.to_csv(args.out_file, sep="\t", index=False)
   
 
if __name__ == "__main__":
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    
    parser.add_argument('-i', '--in-file', type=str, required=True)
    
    parser.add_argument('-o', '--out-file', type=str, required=True)
    
    parser.add_argument('--coverage', type=float, required=True)
    
    parser.add_argument('--tumour-content', type=float, required=True)
    
    parser.add_argument('--clone-prevalences', type=str, required=True)
    
    parser.add_argument('--num-bins', type=str, required=True)
    
    parser.add_argument('--replicate', type=int, required=True)
    
    parser.add_argument('--ichorcna-settings-file', type=str, required=True)
    
    
    cli_args = parser.parse_args()
    
    main(cli_args)