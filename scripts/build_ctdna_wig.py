import pandas as pd 


def df_to_wig(df: pd.DataFrame, genome_build: str, output_path: str) -> None:
    
    with open(output_path, "w") as f:
        
        for c in df['chrom'].unique():
            
            df_c = (
                
                df
                
                .loc[df['chrom'] == c]
                
                .sort_values("start")
                
            )
            
            start = df_c["start"].iloc[0]
            
            end = df_c['end'].iloc[0]
            
            step = int(end - start)
            
            if start == 0:
                
                start = 1
            
            if genome_build == 'hg19':
                chrom = c.replace('chr')
            else:
                chrom = c
             
            # write header    
            
            s = "fixedStep chrom={chrom} start={start} step={step} span={step}\n".format(
                start=start,
                step=step,
                chrom=chrom
            )
            
            f.write(s)
            
            # write reads
            
            for read in df_c['reads'].values:
                
                f.write("{reads}\n".format(reads=int(read)))


def wig_to_df(wig_path):
    
    chroms, starts, ends, chrms, values = [], [], [], [], []

    with open(wig_path, "r") as f:
        
        chrom = None
        
        start = None
        
        step = None
        
        span = None
        
        pos = None

        for line in f:
            
            line = line.strip()
            
            if not line:
                
                continue

            if line.startswith("fixedStep"):    # header line 
                
                parts = dict(item.split("=") for item in line.split()[1:])
                
                chrom = parts["chrom"]
                
                start = int(parts["start"])
                
                step = int(parts["step"])
                
                span = int(parts.get("span", step))
                
                pos = start
                
            else:
                
                val = float(line)
                
                chroms.append(chrom)
                
                starts.append(pos)
                
                ends.append(pos + span)
                
                chrms.append(f"{chrom}")
                
                values.append(val)
                
                pos += step
                
    return pd.DataFrame(
        data={
        'chrom': chrms,
        'start': starts,
        'end': ends,
        'values': values
        }
    )

 
def main(args):

    # load wig to dataframe 
    
    df_wig = wig_to_df(args.wig_template_file)
    
    # load cfdna data
    
    df_cfdna = pd.read_csv(args.in_file, sep='\t')
    
    df = (
        
        df_cfdna
        
        # 1 based index bins 
        
        .assign(
            start=lambda df: df['start'] + 1,     
            end=lambda df: df['end'] + 1
        )
        
        # get missing bins 
        
        .merge(df_wig, on=['chrom', 'start'], how='outer', suffixes=('_cfdna', '')) 
        
        # fill missing bins
        
        .fillna({'reads': 0})
        
        .astype({'reads': int})
        
        .filter(items=['chrom', 'start', 'end', 'reads'])
        
    )
    
    # convert dataframe to wig
    
    df_to_wig(df=df, genome_build='hg38', output_path=args.out_file)
    
    pass

 
if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    
    default0 = "/home/matteo/projects/cfdna/wfs/results/ichorcna-inf-pool-tf-smk/TFRI004/out_dir/input/ctdna/coverage_0/tc_0/cp_0/num_bins_0/replicate_0.tsv.gz"
    
    default1 = "/home/matteo/projects/cfdna/wfs/src/ichorcna-inf-pool-tf-smk/scripts/gc_hg38_500kb.wig"

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)
    
    parser.add_argument("-t", "--wig-template-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
