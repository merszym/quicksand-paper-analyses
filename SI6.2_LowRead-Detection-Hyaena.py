

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
def _(plt):
    from styles import get_palette

    plt.style.use("fast")
    return (get_palette,)


@app.cell
def _(pd):
    df = pd.read_csv("assets/SI6.2/quicksand_v2.3/final_report.tsv", sep="\t")
    return (df,)


@app.cell
def _(Path, pd):
    kmers = pd.DataFrame()

    for _report in Path("assets/SI6.2/quicksand_v2.3/stats/").glob('*.report'):
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
    full['Human Sequences'] = full['RG'].apply(lambda x: int(x.split('_')[-1]))
    full['Assignment'] = full['Family'].apply(lambda x: x if x in ['Hominidae','Hyaenidae'] else 'Other')
    return


@app.cell
def _(full, get_palette, plt, sns):
    _fig = plt.figure(figsize=(12,10))
    _ax1 = _fig.add_subplot(221)
    _ax2 = _fig.add_subplot(222)
    _ax3 = _fig.add_subplot(223)
    _ax4 = _fig.add_subplot(224)

    ## Plot 1

    sns.barplot(
        data=full.groupby('Human Sequences', as_index=False).count(),
        x='Human Sequences',
        y='FamPercentage',
        ax=_ax1
    )

    _ax1.set_title('A', size=21)
    _ax1.set_ylabel('Detected Families', size=15)
    _ax1.set_xlabel('Human Sequences', size=15)

    #_ax1.set_yscale('log')

    ## Plot 2

    sns.scatterplot(
        data=full[full.Family=='Hominidae'],
        x='Human Sequences',
        y='ReadsDeduped',
        ax=_ax2
    )

    _ax2.set_title('B', size=21)
    _ax2.set_ylabel('Final Sequences', size=15)
    _ax2.set_xlabel('Human Sequences', size=15)
    #_ax2.set_yscale('log')

    sns.scatterplot(
        data=full,
        x='FamPercentage',
        y='ProportionExpectedBreadth',
        ax=_ax3,
        hue='Assignment',
        palette=get_palette(3, r=True)
    )
    _ax3.set_title('C', size=21)
    _ax3.set_yscale('log')
    _ax3.set_xscale('log')
    _ax3.axhline(0.5, ls='--', c='red')
    _ax3.axvline(0.5, ls='--', c='red')
    _ax3.set_ylabel('PEB', size=15)
    _ax3.set_xlabel('PSF', size=15)

    sns.scatterplot(
        data=full,
        x='FamPercentage',
        y='kmers',
        ax=_ax4,
        hue='Assignment',
        palette=get_palette(3, r=True)
    )
    _ax4.set_title('D', size=21)
    _ax4.set_yscale('log')
    _ax4.set_xscale('log')
    _ax4.axhline(129, ls='--', c='red')
    _ax4.set_ylabel('Unique Kmers', size=15)
    _ax4.set_xlabel('PSF', size=15)

    _ax1.spines['top'].set_visible(False)
    _ax1.spines['right'].set_visible(False)
    _ax2.spines['top'].set_visible(False)
    _ax2.spines['right'].set_visible(False)
    _ax3.spines['top'].set_visible(False)
    _ax3.spines['right'].set_visible(False)
    _ax4.spines['top'].set_visible(False)
    _ax4.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig('out/Figure_S6.2.svg', dpi=600)
    plt.show()
    return


@app.cell
def _(full):
    full[full.Family=='Hyaenidae']
    return


if __name__ == "__main__":
    app.run()
