

import marimo

__generated_with = "0.13.2"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Compare kmer performance

        Check which of the kmers for the simulated datasets show the highest number of recovered sequences and the best accuracy
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""## Imports and Functions""")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import seaborn as sns
    import pandas as pd
    import re
    import matplotlib.pyplot as plt
    from pathlib import Path
    import matplotlib.patches as mpatches
    return Path, mpatches, pd, plt, re, sns


@app.cell
def _():
    from styles import set_style, get_palette

    set_style()
    return (get_palette,)


@app.cell
def _():
    family_dict={
        "Bovidae":"DQ124371.1",
        "Suidae":"KC469586.1",
        "Hyaenidae":"JF894377.1",
        "Elephantidae":"EU153448.1",
        "Hominidae":"KC879692.1",
    }
    return (family_dict,)


@app.cell
def _(pd):
    def compare_readwise(kmer):
        """
        The function returns the df with the summary of the family bam files.
        A table that shows read of a family bam file, if the sequence came from that family or not, then calculates a false-positive rate. 
        """
        readwise = pd.read_csv(
            f'assets/SI3_SI4/{kmer}/simadna_md_1000_Ziphiidae_{kmer}_readwise', 
            sep=','
        )

        summ = readwise.groupby(["RG","Family","Kmer"], as_index=False).sum()
        summ["FP_Rate"] = summ.False_Positive / summ.All

        return summ
    return (compare_readwise,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Main Loop""")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1. Load the data and prepare new dataframes

        We need to load all the data for the different kmer-runs and check for each read if it was a true-or false positive assignment

        - with the "pure" quicksand output
        - with the filtered output
        """
    )
    return


@app.cell
def _(compare_readwise, mo, pd):
    df = pd.DataFrame()
    rw = pd.DataFrame()

    for kmer in mo.status.progress_bar([f"kmer{x}" for x in [18,19,20,21,22,24,28,30]]):
        _tmp = pd.read_csv(
            f"assets/SI3_SI4/{kmer}/quicksand_v2.3/final_report.tsv", sep="\t")
        _tmp["Kmer"] = kmer
        _rwtmp = compare_readwise(kmer)

        df = pd.concat([df,_tmp], ignore_index=True)
        rw = pd.concat([rw,_rwtmp], ignore_index=True)
    return df, kmer, rw


@app.cell
def _(df, rw):
    # RG is dataset
    df_full = df.merge(rw, on=['RG', 'Family', 'Kmer'], validate="1:1")
    return (df_full,)


@app.cell
def _(df_full):
    df_full["FamPercentage"] = df_full.FamPercentage.astype(float)
    df_full["ProportionExpectedBreadth"] = df_full.ProportionExpectedBreadth.astype(float)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. quicksand accuracy 

        In this section we test, how the assignment accuracy of quicksand differs between the kmer-sizes.

        For all families in the final report, we count the number of bedfiltered sequences assigned to that family. While doing so, we check for each read, if it was assigned correctly or incorrectly.

        Eventually we group the results by RG (dataset) and kmer-size, compare the number of sequences vs. the number of true positive sequences and thus obtain only a single accuracy-value for each dataset and kmer. 

        This analysis is repeated after applying the 0.5% and 0.5 expected breadth filter (false-positive families filter) before grouping
        """
    )
    return


@app.cell
def _():
    damage_profiles_hr={
        'md':'Medium Damage',
        'hd':'High Damage'
    }
    return (damage_profiles_hr,)


@app.cell
def _(damage_profiles_hr, df_full, re):
    # lets go by the number of true positive and false positive sequences instead of families!
    kmer_data = df_full[[
            "RG", "Kmer", "All", "True_Positive",
            "False_Positive",
        ]].groupby(
            ["RG", "Kmer"], as_index=False
        ).sum()

    kmer_data["Accuracy"] = 100 * kmer_data["True_Positive"] / kmer_data["All"]

    kmer_data['Damage Profile'] = kmer_data.RG.apply(
        lambda x: damage_profiles_hr[re.search('[hm]d', x).group()] 
        if bool(re.search('[hm]d', x)) 
        else 'No Damage'
    )
    return (kmer_data,)


@app.cell
def _(damage_profiles_hr, df_full, re):
    # lets go by the number of true positive and false positive sequences instead of families! This time with the filters!
    kmer_data_filtered = df_full[
        (df_full.FamPercentage > 0.5)
        & (df_full.ProportionExpectedBreadth > 0.5)
    ][[
        "RG", "Kmer", "All", "True_Positive",
        "False_Positive","FamPercentage"
    ]].groupby(
        ["RG", "Kmer"], as_index=False
    ).sum()

    kmer_data_filtered["Accuracy"] = 100 * kmer_data_filtered["True_Positive"] / kmer_data_filtered["All"] 
    kmer_data_filtered['Damage Profile'] = kmer_data_filtered.RG.apply(
        lambda x: damage_profiles_hr[re.search('[hm]d', x).group()] 
        if bool(re.search('[hm]d', x)) 
        else 'No Damage'
    )
    return (kmer_data_filtered,)


@app.cell
def _(get_palette, kmer_data, kmer_data_filtered, plt, sns):
    import matplotlib as mpl

    def get_accuracy_plot(data, ax, **kwargs):
        sns.boxplot(
            data=data,
            x='Kmer', 
            y="Accuracy", 
            hue="Damage Profile",
            hue_order=['No Damage', 'Medium Damage', 'High Damage'],
            ax=ax, 
            **kwargs
        )
        return ax


    _fig = plt.figure(figsize=(12,5))
    _ax1 = _fig.add_subplot(1,2,1)
    _ax2 = _fig.add_subplot(2,2,(2,4), sharey=_ax1)
    #_ax3 = _fig.add_axes([0.6, 0.25, 0.3, 0.3])

    _ax1 = get_accuracy_plot(
        kmer_data, 
        _ax1, 
        legend=False, 
        palette=get_palette(3, r=True),
    )

    _ax2 = get_accuracy_plot(
        kmer_data_filtered, 
        _ax2, 
        palette=get_palette(3, r=True),
    )

    #_ax3 = get_accuracy_plot(
    #    kmer_data_filtered, 
    #    _ax3, 
    #    legend=False,
    #    palette=get_palette(3, r=True),
    #)

    plt.setp(_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.setp(_ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    #plt.setp(_ax3.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    #plt.setp(_ax3, xticks=[], )

    _ax1.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.2f}'))
    _ax2.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.2f}'))
    #_ax3.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.2f}'))

    _ax1.set_title("A")
    _ax1.set_xlabel("")
    _ax1.set_ylabel("Accuracy (%)")
    _ax2.set_title("B")
    _ax2.set_xlabel("")
    _ax2.set_ylabel("")
    #_ax3.set_ylim((99.5,100.1))
    #_ax3.set_xlabel("")
    #_ax3.set_ylabel("")


    _ax2.legend(
        prop={'size': 10}, 
        title='Damage', 
        facecolor='white', 
        edgecolor='white', 
        bbox_to_anchor=(1, 1, 0, 0)
    )

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(damage_profiles_hr, df_full, re):
    df_full['damage'] = df_full.RG.apply(
        lambda x: damage_profiles_hr[re.search('[hm]d', x).group()] 
        if bool(re.search('[hm]d', x)) 
        else 'No Damage'
    )

    df_full['percentage_mt'] = df_full.RG.apply(
        lambda x: re.search('[01]+', x).group() 
        if bool(re.search('[01]+', x)) 
        else '1'
    )
    return


@app.cell
def _(df_full, family_dict, get_palette, mpatches, plt, sns):
    _fig = plt.figure(figsize=(12, 5))

    _ax1 = _fig.add_subplot(131)
    _ax2 = _fig.add_subplot(132, sharey=_ax1)
    _ax3 = _fig.add_subplot(133, sharey=_ax1)

    # lets show the false positive families
    _subset = df_full[df_full.Family.isin(family_dict) == False].copy()

    _sizes = [
        '0.1', 
        '1', 
        '10', 
        '100', 
        '1000'
    ]

    _percent_mt = [
        '01', 
        '1', 
        '10', 
        '100', 
        '1000'
    ]

    _barcolors = { size:color for size,color in zip(_sizes, get_palette(5)) }
    _palette = { size:color for size,color in zip(_percent_mt, get_palette(5)) }


    sns.histplot(
        data=_subset[_subset.damage == 'No Damage'], 
        x='Kmer', 
        hue='percentage_mt',
        ax=_ax1,
        legend=False,
        palette = _palette
    )

    _ax1.set_xlabel('')
    _ax1.set_ylabel('Number of Families')
    _ax1.set_title('A', size=22)


    sns.histplot(
        data=_subset[_subset.damage == 'Medium Damage'], 
        x='Kmer', 
        hue='percentage_mt', 
        ax=_ax2,
        legend=False,
        palette = _palette
    )

    _ax2.set_xlabel('')
    _ax2.set_ylabel('')
    _ax2.set_title('B', size=22)


    sns.histplot(
        data=_subset[_subset.damage == 'High Damage'], 
        x='Kmer', 
        hue='percentage_mt', 
        ax=_ax3,
        legend=False,
        palette = _palette
    )
    _ax3.set_xlabel('')
    _ax3.set_ylabel('')
    _ax3.set_title('C', size=22)


    _handles = []

    for (_label, _color) in _barcolors.items():
        _patch = mpatches.Patch(color=_color, label=_label)
        _handles.append(_patch)

    plt.legend(
        handles=_handles, 
        prop={'size': 10}, 
        facecolor='white', 
        edgecolor='white', 
        title='Database Size', 
        bbox_to_anchor=(1, 1, 0, 0)
    )

    plt.setp(_ax1.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    plt.setp(_ax2.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    plt.setp(_ax3.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(Path, damage_profiles_hr, kmer, mo, pd, re):
    # I realized, that this is already filtered by the unique-kmer filter... So lets do it again, but parse the kraken-reports directly

    kraken = pd.DataFrame()

    for _kmer in mo.status.progress_bar([f"kmer{x}" for x in [18,19,20,21,22,24,28,30]]):
        for _report in Path(
            f"assets/SI3_SI4/{kmer}/quicksand_v2.3/stats/").glob("*.report"):
            _tmp = pd.read_csv(_report, sep="\t", skiprows=2)
            _tmp["Kmer"] = _kmer
            _tmp["dataset"] = _report.name.split(".")[0]

            kraken = pd.concat([kraken,_tmp], ignore_index=True)

    kraken = kraken[kraken['rank'] == 'family'].copy()
    kraken['damage'] = kraken.dataset.apply(
        lambda x: damage_profiles_hr[re.search('[hm]d', x).group()] 
        if bool(re.search('[hm]d', x)) 
        else 'No Damage'
    )
    kraken['Family'] = kraken.taxName.apply(lambda x: x.strip())
    kraken['percentage_mt'] = kraken.dataset.apply(
        lambda x: re.search('[01]+', x).group() 
        if bool(re.search('[01]+', x)) 
        else '1'
    )
    return (kraken,)


@app.cell
def _(family_dict, get_palette, kraken, mpatches, plt, sns):
    _fig = plt.figure(figsize=(12, 5))

    _ax1 = _fig.add_subplot(131)
    _ax2 = _fig.add_subplot(132, sharey=_ax1)
    _ax3 = _fig.add_subplot(133, sharey=_ax1)

    # lets show the false positive families
    _subset = kraken[kraken.Family.isin(family_dict) == False].copy()
    _subset.sort_values(['Kmer','percentage_mt'], inplace=True)

    _sizes = [
        '0.1', 
        '1', 
        '10', 
        '100', 
        '1000'
    ]

    _percent_mt = [
        '01', 
        '1', 
        '10', 
        '100', 
        '1000'
    ]

    _barcolors = { size:color for size,color in zip(_sizes, get_palette(5)) }
    _palette = { size:color for size,color in zip(_percent_mt, get_palette(5)) }


    sns.histplot(
        data=_subset[_subset.damage == 'No Damage'], 
        x='Kmer', 
        hue='percentage_mt',
        ax=_ax1,
        legend=False,
        palette = _palette
    )

    _ax1.set_xlabel('')
    _ax1.set_ylabel('Number of Families')
    _ax1.set_title('A', size=22)


    sns.histplot(
        data=_subset[_subset.damage == 'Medium Damage'], 
        x='Kmer', 
        hue='percentage_mt', 
        ax=_ax2,
        legend=False,
        palette = _palette
    )

    _ax2.set_xlabel('')
    _ax2.set_ylabel('')
    _ax2.set_title('B', size=22)


    sns.histplot(
        data=_subset[_subset.damage == 'High Damage'], 
        x='Kmer', 
        hue='percentage_mt', 
        ax=_ax3,
        legend=False,
        palette = _palette
    )
    _ax3.set_xlabel('')
    _ax3.set_ylabel('')
    _ax3.set_title('C', size=22)


    _handles = []

    for (_label, _color) in _barcolors.items():
        _patch = mpatches.Patch(color=_color, label=_label)
        _handles.append(_patch)

    plt.legend(
        handles=_handles, 
        prop={'size': 10}, 
        facecolor='white', 
        edgecolor='white', 
        title='Dataset Size', 
        bbox_to_anchor=(1, 1, 0, 0)
    )

    plt.setp(_ax1.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    plt.setp(_ax2.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    plt.setp(_ax3.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""## 2. quicksand sensitivity""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Check Kmer-performance

        Here we want to check, which kmer-size gives us most sequences back. We are interested in 

        1. the total number of reads per kmer and "real" family
        2. the proportion of false positives within the families

        We ignore false positive family assignments for now
        """
    )
    return


@app.cell
def _(plt, sns):
    hue_order = ["Hominidae","Suidae","Bovidae","Hyaenidae","Elephantidae"]


    def plot_kmerperformance(
        df, title=None, ax=None, colors=None, subplot="A", ylabel=True):

        fam = df.copy()

        def get_best(rg, f):
            tmp = fam[(fam.RG == rg)&(fam.Family == f)].copy()
            return max(tmp['ReadsDeduped'])

        fam['best'] = fam[['RG','Family']].apply(lambda x: get_best(x[0], x[1]), axis=1)
        fam["relative_best"] = fam['ReadsDeduped'] / fam['best']

        ax = sns.boxplot(
             data=fam, x="Kmer", y="relative_best", 
             hue="Family", 
             hue_order=hue_order, 
             ax=ax, 
             palette=colors,
        )

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        ax.set_title(title, size=22)
        ax.set_ylabel("Relative Sensitivity")
        ax.set_xlabel("")
        if not ylabel:
            ax.set_ylabel("")
        ax.get_legend().set_visible(False)

        return ax
    return (plot_kmerperformance,)


@app.cell
def _(df_full, mpatches, plot_kmerperformance, plt, sns):
    import matplotlib.colors

    _families = ['Hominidae', 'Suidae', 'Bovidae', 'Hyaenidae', 'Elephantidae']

    _colors = {fam:color for (fam, color) in zip(_families, [matplotlib.colors.to_rgb(x) for x in sns.color_palette('Spectral', len(_families)).as_hex()])}

    _fig = plt.figure(figsize=(12, 5))
    _ax1 = _fig.add_subplot(131)
    _ax2 = _fig.add_subplot(132, sharey=_ax1)
    _ax3 = _fig.add_subplot(133, sharey=_ax1)

    _ax1 = plot_kmerperformance(
        df_full[(df_full.damage == 'No Damage') 
        & df_full.Family.isin(_families)], 
        title='A', 
        ax=_ax1, 
        colors=_colors
    )

    _ax2 = plot_kmerperformance(
        df_full[(df_full.damage == 'Medium Damage') 
        & df_full.Family.isin(_families)], 
        title='B', 
        ax=_ax2, 
        colors=_colors, 
        ylabel=False
    )

    _ax3 = plot_kmerperformance(
        df_full[(df_full.damage == 'High Damage') 
        & df_full.Family.isin(_families)], 
        title='C', 
        ax=_ax3, 
        colors=_colors, 
        ylabel=False
    )

    _handles = []
    for (_label, _color) in _colors.items():
        _patch = mpatches.Patch(color=_color, label=_label)
        _handles.append(_patch)

    plt.legend(handles=_handles, prop={'size': 10}, facecolor='white', edgecolor='white', bbox_to_anchor=(1, 1, 0, 0))

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(df_full, get_palette, plt, sns):
    def plot_kmerperformance_merged(
        df, title=None, ax=None, subplot="A", ylabel=True, **kwargs
    ):

        fam = df.copy()

        def get_best(rg, f):
            tmp = fam[(fam.RG == rg)&(fam.Family == f)].copy()
            return max(tmp['ReadsDeduped'])

        fam['best'] = fam[['RG','Family']].apply(lambda x: get_best(x[0], x[1]), axis=1)
        fam["relative_best"] = fam['ReadsDeduped'] / fam['best']

        ax = sns.boxplot(
             data=fam, x="Kmer", y="relative_best", ax=ax, **kwargs 
        )

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        ax.set_title(title, size=15)
        ax.set_ylabel("Relative Sensitivity")
        ax.set_xlabel("")
        if not ylabel:
            ax.set_ylabel("")

        return ax

    _families = ['Hominidae', 'Suidae', 'Bovidae', 'Hyaenidae', 'Elephantidae']

    _fig = plt.figure(figsize=(12, 5))
    _ax1 = _fig.add_subplot(131)
    _ax2 = _fig.add_subplot(132, sharey=_ax1)
    _ax3 = _fig.add_subplot(133, sharey=_ax1)

    _ax1 = plot_kmerperformance_merged(
        df_full[(df_full.damage == 'No Damage') 
        & df_full.Family.isin(_families)], 
        title='A', 
        ax=_ax1, 
        color = get_palette(1)[0]
    )

    _ax2 = plot_kmerperformance_merged(
        df_full[(df_full.damage == 'Medium Damage') 
        & df_full.Family.isin(_families)], 
        title='B', 
        ax=_ax2, 
        ylabel=False,
        color = get_palette(1)[0]
    )

    _ax3 = plot_kmerperformance_merged(
        df_full[(df_full.damage == 'High Damage') 
        & df_full.Family.isin(_families)], 
        title='C', 
        ax=_ax3, 
        ylabel=False,
        color = get_palette(1)[0]
    )


    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(df_full, family_dict):
    # I want it also as a table
    df_summary = df_full[df_full.Family.isin(family_dict)].groupby(
        ["damage","RG","Kmer","Family"], as_index=False
    )[["All","False_Positive","True_Positive"]].sum()

    df_summary["Display"] = df_summary[["All","False_Positive","True_Positive"]].apply(
        lambda x: f"{x[0]} ({x[2]}/{x[1]})" ,axis=1
    )
    return (df_summary,)


@app.cell
def _(df_summary):
    df_res = df_summary.pivot(
        columns="Family", index=["damage","RG","Kmer"], values="Display" 
    )
    df_res = df_res.reset_index()
    df_res[df_res.Kmer=='kmer22']
    return


if __name__ == "__main__":
    app.run()
