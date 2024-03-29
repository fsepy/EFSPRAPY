from concurrent.futures import ProcessPoolExecutor
from os import path

import matplotlib.pyplot as plt
import numpy as np
from sfeprapy.func.xlsx import xlsx_to_dict
from tqdm import tqdm


# plt.style.use('seaborn-v0_8-talk')


def format_ax(
        ax,

        xlabel: str = None,
        xlabel_fontsize=None,
        xticks=None,
        xticks_minor=None,
        xticklabels=None,
        xlim=None,
        xscale=None,

        ylabel: str = None,
        ylabel_fontsize=None,
        yticks=None,
        yticks_minor=None,
        yticklabels=None,
        ylim=None,
        yscale=None,

        ticks_labelsize=None,
        ticklabel_format_style=None,
        ticklabel_format_axis='y',
        ticklabel_format_scilimits=(0, 0),

        legend_title: str = None,
        legend_loc: int = 0,
        legend_ncol=1,
        legend_fontsize=None,
        legend_title_fontsize=None,
        legend_visible: bool = True,
        legend_borderpad=0.4,
        legend_labelspacing=0.5,

        grid_which='both',
        grid_ls='--',
        grid_lw=.5
):
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=xlabel_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=ylabel_fontsize)
    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)
    if xticks_minor is not None:
        ax.set_xticks(xticks_minor, minor=True)
    if yticks_minor is not None:
        ax.set_yticks(yticks_minor, minor=True)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if legend_visible is True:
        ax.legend(
            title=legend_title, loc=legend_loc, ncol=legend_ncol, frameon=True, fontsize=legend_fontsize,
            fancybox=False, title_fontsize=legend_title_fontsize, facecolor=(1, 1, 1, 0.5), edgecolor=(0, 0, 0),
            borderpad=legend_borderpad, labelspacing=legend_labelspacing
        ).set_visible(legend_visible)
    ax.grid(which=grid_which, ls=grid_ls, lw=grid_lw)
    ax.tick_params(labelsize=ticks_labelsize)

    if ticklabel_format_style is not None:
        ax.ticklabel_format(
            style=ticklabel_format_style,
            scilimits=ticklabel_format_scilimits,
            axis=ticklabel_format_axis,
            useMathText=True
        )

    if ticklabel_format_style is not None:
        ax.yaxis.offsetText.set_fontsize('x-small')


def plot_contour(
        ax,
        xx,
        yy,
        zz,

        xticks_minor=None,
        levels=None,
        clabel_fmt=lambda x: f'{x:.0f}',
        clabel_manual=False,
        cmap: str = 'Greys',
        **kwargs
):
    cmap = plt.get_cmap('viridis')
    norm = np.linspace(min(levels), max(levels), len(levels)) / max(levels)
    colors = list(cmap(norm))

    cs = ax.contour(xx, yy, zz, levels=levels, linewidths=0.5,
                    colors='k', linestyles='dotted', antialiased=True)
    cf = ax.contourf(xx, yy, zz, levels=levels, cmap=None,
                     alpha=0.6, extend='neither', colors=colors)
    ax.clabel(cs, cs.levels, inline=True, fmt=clabel_fmt, fontsize='x-small', manual=clabel_manual,
              use_clabeltext=True)
    format_ax(ax=ax, **kwargs)
    if xticks_minor is not None:
        ax.set_xticks(xticks_minor, minor=True)

    return cf


def read_output_and_gen_xyz_worker(args):
    print(args)
    fp_output, W, H, N = args
    try:
        mcs_output = np.genfromtxt(fp_output, delimiter=',', skip_header=1)
        t_ig_ftp = mcs_output[:, 1]
        return W, H, sum(np.logical_and(t_ig_ftp > 0, t_ig_ftp < np.inf)) / float(N)
    except:
        return np.array([W, H, np.nan])


def read_output_and_gen_xyz(fp_input, max_workers=8):
    mcs_input = xlsx_to_dict(fp_input)

    args = list()
    for k, v in mcs_input.items():
        W, H, _ = k.split('-')
        fp_output = path.join(path.dirname(
            fp_input), f'{path.splitext(path.basename(fp_input))[0]}.out', f'{k}.csv')
        args.append((fp_output, float(W), float(H), int(v['n_simulations'])))

    with ProcessPoolExecutor(max_workers=max_workers) as p:
        results = list(
            tqdm(p.map(read_output_and_gen_xyz_worker, args), total=len(args)))
    return np.array(results)
