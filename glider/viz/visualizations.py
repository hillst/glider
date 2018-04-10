import matplotlib as mpl
from matplotlib import pyplot
import seaborn as sb
from glider.utils.trinucs import build_trinuc_lookup, trinuc_order, trinuc_colors
import numpy as np
import pandas as pd

def violins(df, variable, saveprefix):
    import matplotlib as mpl
    from matplotlib import pyplot
    import seaborn as sb

    # filt = df.query('frag<870 & frag>120 & mrbq > 20 & vbq > 20 & mapq > 40')
    ax = pyplot.subplots(figsize=(19, 10))
    sb.violinplot(x="read_source", y=variable, hue="biotype", data=df, palette="muted")
    pyplot.title("Fragment size (Controls interesected with all tumor VCFs)")
    return ax


def biotype_hists(df, variable, hue="biotype", bins=None):
    import matplotlib as mpl
    from matplotlib import pyplot
    import seaborn as sb
    fig, ax = pyplot.subplots(figsize=(15, 10))
    for cat in df[hue].unique():
        sb.distplot(df[df[hue] == cat][variable], hist=True, bins=bins, label="{} == {}".format(hue, cat))
    return ax


def trinuc_odds(df,):
    import matplotlib as mpl
    from matplotlib import pyplot
    import seaborn as sb
    from glider.utils.trinucs import build_trinuc_lookup, trinuc_order, trinuc_colors
    import numpy as np
    import pandas as pd


    sb.set_style("white")
    sb.despine()
    fig, ax = pyplot.subplots(figsize=(30, 20))
    index2trinuc, trinuc2index = build_trinuc_lookup()

    # filtered set, into trinuc only
    trinuc_cols = [index2trinuc[index] for index in range(96)] + ["read_source", "mutation_source", "mutation_callset",
                                                                  "biotype"]
    trinuc_df = df[trinuc_cols]

    trinuc_df = trinuc_df.groupby(["read_source", "mutation_source", "biotype"],
                                  squeeze=True).sum()  # gives us the sum of all the trinucleotides in these groups
    trinuc_df = trinuc_df.div(trinuc_df.sum(axis=1), axis=0)  # normed
    trinuc_df = trinuc_df.reset_index()
    melted = pd.melt(trinuc_df, id_vars=["read_source", "mutation_source", "biotype"])

    luad = melted[melted["biotype"] == "luad"]
    control = melted[melted["biotype"] == "control"]
    #computes odds ratios
    tumor_odds_ratios = [luad[luad["variable"] == tn].mean().value / control[control["variable"] == tn].mean().value for tn in
                     trinuc_order()]
    ax.bar(range(len(trinuc_order())), tumor_odds_ratios, tick_label=trinuc_order(), edgecolor='black', color=trinuc_colors())

    control_odds_ratios = [-1 * control[control["variable"] == tn].mean().value / luad[luad["variable"] == tn].mean().value for tn
                     in trinuc_order()]
    ax.bar(range(len(trinuc_order())), control_odds_ratios, tick_label=trinuc_order(), edgecolor='black', color=trinuc_colors())
    ax.set_yticklabels([abs(tick) for tick in ax.get_yticks() ])


    for item in ax.get_xticklabels():
        item.set_rotation(90)

    ax.set_ylabel("Mutation Odds ratio (LUAD / Control)", fontsize=24)
    ax.set_ylim(-2.0, 2.0)


    # ax.title("MAPQ > 1 MRBQ > 1 VBQ > 1 FRAG > 1 PIR > 0")

def trinuc_boxes(df, saveprefix, biotype="control"):
    import matplotlib as mpl
    from matplotlib import pyplot
    import seaborn as sb
    from glider.utils.trinucs import build_trinuc_lookup, trinuc_order, trinuc_colors
    import numpy as np
    import pandas as pd

    #styling
    sb.set_style("white")
    sb.despine()
    trinuc2index, index2trinuc = build_trinuc_lookup()

    # filtered set, into trinuc only
    trinuc_cols = [index2trinuc[index] for index in range(96)] + ["read_source", "mutation_source", "mutation_callset",
                                                                  "biotype"]
    trinuc_df = df[trinuc_cols]

    # organize
    trinuc_df = trinuc_df.groupby(["read_source", "mutation_source", "biotype"],
                                  squeeze=True).sum()  # gives us the sum of all the trinucleotides in these groups
    # normalize
    plot_view = trinuc_df.div(trinuc_df.sum(axis=1), axis=0)  # normed
    grouped = plot_view.reset_index()
    # this is basically the magic, we group by biotype here
    controls_only = grouped[grouped["biotype"] == biotype]
    melted = pd.melt(controls_only, id_vars=["read_source", "mutation_source", "biotype"])

    # this stuff is all specific to plotting regardless of grouping/subsets,
    # note that boxes vs bars work different dependign on datatype! If you only have 1 example per feature
    # (like in luad) then it doesnt work super great
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
    return ax

def qc_visualization_routine(filt, variable='vbq', bins=None):
    """
    Plots the histograms by biotype and accompanies them with a violin plot.
    Assumes we are using matplotlib inline or something similar as it does not return any axes.


    :param filt:
    :return:
    """


    ax = biotype_hists(filt, variable=variable, hue="biotype", bins=bins)
    pyplot.title("{} ".format(variable))
    pyplot.legend()
    sb.despine()

    ax = biotype_hists(filt[filt["biotype"] == "luad"], variable=variable, hue="batch", bins=bins)
    pyplot.title("{} organized by batch, luad only.".format(variable))
    pyplot.legend()
    sb.despine()

    ax = biotype_hists(filt[filt["biotype"] == "control"], variable=variable, hue="batch", bins=bins)
    pyplot.title("{} organized by batch, control only.".format(variable))
    pyplot.legend()
    sb.despine()

    ax = violins(filt, variable=variable, saveprefix="source_violin.vbq.png")
    sb.despine()
