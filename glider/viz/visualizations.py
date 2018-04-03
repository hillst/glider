import matplotlib as mpl

mpl.use('agg')
from matplotlib import pyplot
import seaborn as sb
import pandas as pd


def violins(df, variable, saveprefix):
    # filt = df.query('frag<870 & frag>120 & mrbq > 20 & vbq > 20 & mapq > 40')
    filt = df
    ax = pyplot.subplots(figsize=(19, 10))
    sb.violinplot(x="read_source", y=variable, hue="biotype", data=filt, palette="muted")
    pyplot.title("Fragment size (Controls interesected with all tumor VCFs)")
    return ax
    pyplot.savefig(saveprefix)


def biotype_hists(df, variable, saveprefix):
    # filt = df.query('frag<870 & frag>120 & mrbq > 20 & vbq > 20 & mapq > 40')
    filt = df
    ax = pyplot.subplots(figsize=(19, 10))
    sb.violinplot(x="read_source", y=variable, hue="biotype", data=filt, palette="muted")
    pyplot.title("Fragment size (Controls interesected with all tumor VCFs)")
    return ax
    pyplot.savefig(target)


def trinuc_boxes(df, saveprefix, biotype="control"):
    from glider.utils.trinucs import build_trinuc_lookup, trinuc_order, trinuc_colors
    import numpy as np
    sb.set_style("white")
    sb.despine()

    index2trinuc, trinuc2index = build_trinuc_lookup()

    # filtered set, into trinuc only
    trinuc_cols = [index2trinuc[index] for index in range(96)] + ["read_source", "mutation_source", "mutation_callset",
                                                                  "biotype"]

    trinuc_df = df[trinuc_cols]
    trinuc_df["read_source"]

    # organize
    trinuc_df = trinuc_df.groupby(["read_source", "mutation_source", "biotype"],
                                  squeeze=True).sum()  # gives us the sum of all the trinucleotides in these groups
    # normalize
    plot_view = trinuc_df.div(trinuc_df.sum(axis=1), axis=0)  # normed
    grouped = plot_view.reset_index()
    # this is basically the magic, we group by biotype here
    controls_only = grouped[grouped["biotype"] == biotype]
    melted = pd.melt(controls_only, id_vars=["read_source", "mutation_source", "biotype"])

    # this stuff is all specific to plotting regardless of grouping/subsets, note that boxes vs bars work different dependign on datatype! If you only have 1 example per feature (like in luad) then it doesnt work super great
    this_sample = []
    for trinuc in trinuc_order():
        this_sample.append(np.asarray(melted[melted["variable"] == trinuc]["value"]).flatten())

    fig, ax = pyplot.subplots(figsize=(30, 20))
    bplot = ax.boxplot(x=this_sample, labels=trinuc_order(), patch_artist=True, showmeans=True)
    for item in ax.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(18)

    for patch, color in zip(bplot['boxes'], trinuc_colors()):
        patch.set_facecolor(color)

    pyplot.ylim(0, 0.15)
    # pyplot.savefig(saveprefix)
    return ax
