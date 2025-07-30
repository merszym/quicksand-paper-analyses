

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from pathlib import Path
    return Path, pd, plt, sns


@app.cell
def _():
    from styles import set_style, get_palette

    set_style()
    return (get_palette,)


@app.cell
def _(pd):
    df = pd.read_csv("assets/SI6.1/quicksand_v2.3/final_report.tsv", sep='\t')
    return (df,)


@app.cell
def _(Path, pd):
    kmers = pd.DataFrame()

    for _report in Path("assets/SI6.1/quicksand_v2.3/stats/").glob('*.report'):
        _tmp = pd.read_csv(_report, sep="\t", skiprows=2)
        _tmp = _tmp[_tmp['rank']=='family']
        _tmp['RG'] = _report.name.replace(".kraken.report",'')
        kmers = pd.concat([kmers, _tmp], ignore_index=True)
    return (kmers,)


@app.cell
def _(kmers):
    kmers['Family'] = kmers.taxName.apply(lambda x: x.strip())
    return


@app.cell
def _(df, kmers):
    full = df.merge(kmers, on=['RG','Family'], how='left', validate="1:1")
    full = full.convert_dtypes()
    full['ReadsRaw'] = full["ReadsRaw"].astype(int)
    full = full[full.Family != "-"].copy()
    return (full,)


@app.cell
def _(full):
    full['3term'] = full['Deam3(95ci)'].apply(lambda x: float(x.split()[0]) if x.startswith('N')==False else None)
    full['5term'] = full['Deam5(95ci)'].apply(lambda x: float(x.split()[0]) if x.startswith('N')==False else None)
    return


@app.cell
def _(full):
    _subset = full[(full.FamPercentage >= 0.5)&(full.ProportionExpectedBreadth >= 0.5)]
    _subset[_subset.Family == 'Muridae']
    return


@app.cell
def _(full, get_palette, plt, sns):
    _subset = full[(full.FamPercentage >= 0.5)&(full.ProportionExpectedBreadth >= 0.5)]

    _fig = plt.figure(figsize=(10,8))

    _ax1 = _fig.add_subplot(221)
    _ax2 = _fig.add_subplot(222)
    _ax3 = _fig.add_subplot(223)
    _ax4 = _fig.add_subplot(224)
    #_ax5 = _fig.add_subplot(325)
    #_ax6 = _fig.add_subplot(326)

    # Plot 1
    sns.lineplot(
        data=_subset,
        x='ReadsRaw',
        y="ReadsDeduped",
        hue='Family',
        ax=_ax1,
        palette=get_palette(3, r=True)
    )
    sns.scatterplot(
        data=_subset[_subset.Family == 'Muridae'],
        x='ReadsRaw',
        y="ReadsDeduped",
        ax=_ax1,
        color=get_palette(3, r=True)[-1]
    )

    _ax1.set_title("A")
    _ax1.set_ylabel('Final Sequences')
    _ax1.set_xlabel('Number of Sequences')

    # Plot 2
    sns.lineplot(
        data=_subset,
        x='ReadsRaw',
        y="kmers",
        hue='Family',
        ax=_ax2,
        palette=get_palette(3, r=True)
    )
    sns.scatterplot(
        data=_subset[_subset.Family == 'Muridae'],
        x='ReadsRaw',
        y="kmers",
        ax=_ax2,
        color=get_palette(3, r=True)[-1]
    )

    _ax2.set_title("B")
    _ax2.set_ylabel('Unique Kmers')
    _ax2.set_xlabel('Number of Sequences')
    _ax2.axhline(129, ls='--', c='red')

    # Plot 3
    sns.lineplot(
        data=_subset[_subset.Family=='Hominidae'].groupby(['ReadsRaw','Family','Ancientness'], as_index=False).count(),
        x='ReadsRaw',
        y="ReadsMapped",
        hue='Ancientness',
        ax=_ax3,
        palette=get_palette(3, r=True)
    )

    _ax3.set_title("C")
    _ax3.set_ylabel('Number of Samples')
    _ax3.set_xlabel('Number of Sequences')

    # plot 4

    # plot 5
    sns.lineplot(
        data=_subset[_subset.Family=='Hominidae'],
        x='ReadsRaw',
        y="3term",
        ax=_ax4,
        color=get_palette(2, r=True)[0],
        label="3'C->T"
    )
    sns.lineplot(
        data=_subset[_subset.Family=='Hominidae'],
        x='ReadsRaw',
        y="5term",
        ax=_ax4,
        color=get_palette(2)[0],
        label="5'C->T"
    )
    _ax4.legend()
    _ax4.set_title("D")
    _ax4.set_ylabel('Term. Deamination Rate')
    _ax4.set_xlabel('Number of Sequences')

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(full):
    full[full.Ancientness=='++']
    return


if __name__ == "__main__":
    app.run()
