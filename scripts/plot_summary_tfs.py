import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error


plot_settings = {
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "axes.titleweight": "bold",
    "font.size": 5,
}

def plot(
    df_plot: pd.DataFrame,
    xlims: tuple[float, float] | None = None,
    ylims: tuple[float, float] | None = None,
    log_scale: bool = False,
    figsize: tuple[float, float] | None = None,
    show: bool = False,
    output_path: str | None = None,
    add_inits: bool = False,
) -> None:
    with plt.rc_context(plot_settings):
        _plot(
            df_plot=df_plot,
            xlims=xlims,
            ylims=ylims,
            log_scale=log_scale,
            figsize=figsize,
            output_path=output_path,
            add_inits=add_inits,
            show=show
        )

def _plot(
    df_plot: pd.DataFrame,
    xlims: tuple[float, float] | None = None,
    ylims: tuple[float, float] | None = None,
    log_scale: bool = False,
    figsize: tuple[float, float] | None = None,
    show: bool = False,
    add_inits: bool = False,
    output_path: str | None = None,
) -> None:
    """
    Each subplot is a different coverage.
    
    Args:
        df_plot (pd.DataFrame): with the following columns:
            [
                'coverage',
                'tumour_content',
                'clone_prevalence',
                'num_bins',
                'replicate',
                'init',                 # ichor
                'tf_est',               # ichor
                'phi_est',              # ichor
                'loglik',               # ichor
            ]
    """
    coverages = df_plot['coverage'].unique().tolist()
    
    num_covs = len(coverages)
    
    num_cols = int(np.ceil(np.sqrt(num_covs)))
    
    num_rows = int(np.ceil(num_covs / num_cols))
    
    if figsize is None:
        
        figsize = (num_cols * 3, num_rows * 3)
    
    fig = plt.figure(figsize=figsize, constrained_layout = True)
    
    gs = fig.add_gridspec(nrows=num_rows, ncols=num_cols)
    
    for i in range(num_rows):
        
        for j in range(num_cols):
            
            ax = fig.add_subplot(gs[i, j])
            
            idx = i * num_cols + j
            
            if idx < num_covs:
                
                df = df_plot[df_plot['coverage'] == coverages[idx]]
                
                plot_tf_estimates(
                    df=df,
                    ax=ax,
                    xlims=xlims,
                    ylims=ylims,
                    log_scale=log_scale,
                    add_inits=add_inits
                )
                
            else:
                
                ax.axis('off')
                
    if output_path is not None:
        
        plt.savefig(output_path)
    
    if show:
        
        plt.show()
        
    plt.close()



def plot_tf_estimates(
    df: pd.DataFrame,
    ax: plt.Axes,
    xlims: tuple[float, float] | None = None,
    ylims: tuple[float, float] | None = None,
    log_scale: bool = False,
    add_inits: bool = False
):
    """
    Args:
        df (pd.DataFrame): dataframe with columns 
            [
                'coverage',
                'tumour_content',
                'replicate',
                'init',                 # ichor
                'tf_est',               # ichor
                'phi_est',              # ichor
                'loglik',               # ichor
                # 'mean',                 # cfclone
                # 'lower_hdi',            # cfclone
                # 'upper_hdi',            # cfclone
                # 'median',               # cfclone
            ]
    """
    # make sure dataframe is valid 
    
    covs = df['coverage'].unique().tolist()
    
    assert len(covs) == 1
    
    # ADD PERFECT INFERENCE LINE 
    
    add_y_equals_x(ax=ax)
    
    df_max = (
        df.
        
        sort_values(
            by=['coverage', 'tumour_content', 'replicate', 'loglik'], 
            ascending=True
        )
        
        .drop_duplicates(
            subset=['coverage', 'tumour_content', 'replicate'], 
            keep='last'
        )
        
        .reset_index(drop=True)
    )
    
    x = df_max['tumour_content'].values
    
    y = df_max['tf_est'].values
    
    # ADD OLS OF ICHORCNA
    
    add_ols(
        x=x,
        y=y,
        ax=ax,
        label='ichorcna',
        color='green',
        text_coords=(0.8, 0.05)
    )
    
    # IF TRUE ADD ALL EM INITS
    
    if add_inits:
        
        add_different_inits(df, ax)
    
    add_legend(ax=ax)

    ax.set_title('Coverage: {}X'.format(covs[0]))
    

        
    if xlims is not None:
        
        ax.set_xlim(xlims)
        
    if ylims is not None:
        
        ax.set_ylim(ylims)
        
        
    if log_scale:
        
        ax.set_xlabel('log expected tumour fraction')
        
        ax.set_ylabel('log estimated tumour fraction')
        
        ax.set_xscale('log')
        
        ax.set_yscale('log')
        
    else:
        
        ax.set_xlabel('expected tumour fraction')
        
        ax.set_ylabel('estimated Tumour Fraction')
        
        
def add_different_inits(df, ax):
    
    g = df.groupby(['coverage', 'tumour_content', 'replicate', 'init'])
    
    for (cov, tc, r, init), df in g:
        
        phi_est = df['phi_est']
        
        ax.scatter(
            x=tc,
            y=df['tf_est'],
            alpha=0.1,
            color='green'
        )

 
def add_y_equals_x(ax: plt.Axes):
    
    line_min, line_max = ax.get_xlim()
    
    ax.plot(
        [line_min, line_max],
        [line_min, line_max],
        color='red',
        linestyle='--',
        linewidth=1.5,
        alpha=0.8,
        label="y = x"
    )
    
    
def add_legend(ax: plt.Axes):
    
    handles, labels = ax.get_legend_handles_labels()
    
    keep = ["cfclone", "ichorcna"]
    
    filtered = [(h, l) for h, l in zip(handles, labels) if l in keep]
    
    if filtered:
        
        ax.legend(*zip(*filtered))
        

def add_ols(
    x: np.ndarray,
    y: np.ndarray,
    ax: plt.Axes,
    label: str | None = None,
    color: str | None = None,
    text_coords: tuple[float, float] = (0.5, 0.5),
) -> tuple[float, float, float]:
    colours = get_colours(y, color)

    y_scatter = get_y_scatter(y)
    # plot data 
    ax.scatter(
        x=x, 
        y=y_scatter,
        c=colours,
        s=40,
        alpha=0.7,
        linewidths=0.5,
        label=label,
        zorder=3
    )
    
    # plot ols
    nan_idxs = np.isnan(y)
    min_x, max_x = ax.get_xlim()
    x_ols = x[~nan_idxs]
    y_ols = y[~nan_idxs]
    
    try:
        coeffs = np.polyfit(x_ols, y_ols, 1)
        fit_y = np.polyval(coeffs, [min_x, max_x])
        ax.plot(
            [min_x, max_x],
            fit_y,
            color=color,
            linewidth=1.8,
            alpha=0.5,
            label=f"Least squares: y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}",
            zorder=2
        )
        # add text of regression metrics
        r2 = np.corrcoef(x_ols, y_ols)[0, 1] ** 2  # R^2 via correlation
        
        mae = mean_absolute_error(x_ols, y_ols)
        max_ae = np.max(np.abs(x_ols - y_ols))
        
        ax.text(
            x=text_coords[0],
            y=text_coords[1],  # relative position in axes coordinates
            s=f"$R^2$ = {r2:.3f} \n  MAE = {mae:.3f} \n  MaxAE = {max_ae:.3f}",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            color=color,
        )
    except:
        print("OLS didn't converge")
        
        
def get_colours(y: np.ndarray, colour: str, nan_colour: str = 'red') -> np.ndarray:
        nan_idxs = np.isnan(y)
        colours = np.full(y.shape[0], colour)
        colours[nan_idxs] = nan_colour
        return colours
        

def get_y_scatter(y: np.ndarray, fill_nan: float = -0.05) -> np.ndarray:
    nan_idxs = np.isnan(y)
    if np.any(nan_idxs):
        y_scatter = np.where(nan_idxs, fill_nan, y)
    else:
        y_scatter = y
    return y_scatter


def main(args):
    
    df = pd.read_csv(args.in_file, sep="\t")
    
    df['tf_est'] = 1. - df['n_est']
    
    plot(
        df_plot=df,
        output_path=args.out_file, 
        add_inits=args.add_inits,
        log_scale=args.log_scale,
    )
        



if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", type=str, required=True)

    parser.add_argument("-o", "--out-file", type=str, required=True)
    
    parser.add_argument("--add-inits", action="store_true")
    
    parser.add_argument("--log-scale", action="store_true")

    cli_args = parser.parse_args()

    main(cli_args)
