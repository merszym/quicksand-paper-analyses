

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path
    import re
    import matplotlib.ticker as ticker
    import marimo as mo
    return Path, mo, pd, plt, re, sns, ticker


@app.cell
def _():
    from styles import set_style, get_palette

    set_style()
    return (get_palette,)


@app.cell
def _():
    damage_profiles_hr={
        'md':'Medium Damage',
        'hd':'High Damage'
    }
    return (damage_profiles_hr,)


@app.cell
def _(Path, pd):
    df = pd.DataFrame()

    for kmer in [f"kmer{x}" for x in (18,19,20,21,22,24,28,30)]:
        _p = Path(f"assets/SI3_SI4/{kmer}/quicksand_v2.3/final_report.tsv")
        _tmp = pd.read_csv(_p, sep="\t")
        _tmp["kmer"] = kmer
        df = pd.concat([df, _tmp], ignore_index=True)
    return df, kmer


@app.cell
def _(Path, damage_profiles_hr, kmer, mo, pd, re):
    kraken = pd.DataFrame()

    for _kmer in mo.status.progress_bar([f"kmer{x}" for x in [18,19,20,21,22,24,28,30]]):
        for _report in Path(f"assets/SI3_SI4/{kmer}/quicksand_v2.3/stats/").glob("*.report"):
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
    kraken["Real"] = kraken["Family"].apply(lambda x: x in ["Hominidae","Elephantidae","Hyaenidae","Suidae","Bovidae"])
    kraken['Assignment'] = kraken['Real'].apply(lambda x: "True Positive" if x else "False Positive")
    return (kraken,)


@app.cell
def _(mo):
    mo.md(r"""## Uniq Kmer Filter""")
    return


@app.cell
def _(get_palette, kraken, plt, sns):
    _fig = plt.figure(figsize=(11,5))

    _ax1 = _fig.add_subplot(121)
    _ax2 = _fig.add_subplot(122)

    sns.scatterplot(
        data=kraken[kraken.Kmer == 'kmer22'],
        y="kmers",
        x="reads",
        hue="Assignment",
        palette=get_palette(2, r=True),
        hue_order=["True Positive","False Positive"],
        ax=_ax1
    )

    _ax1.set_title("A")
    _ax1.set_xscale("log")
    _ax1.set_yscale("log")
    _ax1.set_ylabel("Number of Unique Kmers")
    _ax1.set_xlabel("Number of Reads")
    #_ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}x'))

    _ax1.axhline(y=129, color='red', lw=.5, ls='-.')
    #_ax1.axvline(x=3, color='red', lw=.5, ls='-.')

    _g = sns.boxplot(
        data=kraken[kraken.Kmer == 'kmer22'],
        #y="reads",
        y="kmers",
        x="Assignment",
        hue="percentage_mt",
        palette=get_palette(5, r=True),
        hue_order=['01','1','10','100','1000'],
        order=['True Positive', 'False Positive'],
        ax=_ax2
    )

    _ax2.set_title("B")
    _ax2.set_yscale("log")
    _ax2.set_ylabel("Number of Unique Kmers")
    _ax2.set_xlabel("Family Assignments")


    _ax2.legend(
        prop={'size': 10}, 
        title='Dataset Size', 
        facecolor='white', 
        edgecolor='white', 
        bbox_to_anchor=(1, 1, 0, 0)
    )


    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(kraken):
    kraken[kraken.Assignment == "False Positive"].groupby("percentage_mt")['kmers'].describe()
    return


@app.cell
def _(mo):
    mo.md(r"""## Percentage Filter""")
    return


@app.cell
def _(df):
    filtered_df = df[df.ReadsBedfiltered > 0].copy()
    return (filtered_df,)


@app.cell
def _(filtered_df):
    filtered_df["Real"] = filtered_df["Family"].apply(lambda x: x in ["Hominidae","Elephantidae","Hyaenidae","Suidae","Bovidae"])
    return


@app.cell
def _(filtered_df):
    filtered_df[(filtered_df.Real==False)&(filtered_df.FamPercentage > 1)]
    return


@app.cell
def _(filtered_df, re):
    damage_dict={
        'md':'Medium Damage',
        'hd':'High Damage'
    }

    filtered_df['dataset'] = filtered_df.RG.apply(lambda x: "1" if not "1" in x else x.replace(".fq","").split("_")[-1])
    filtered_df['damage'] = filtered_df.RG.apply(lambda x: re.search("(md|hd)",x).group() if bool(re.search("(md|hd)",x)) else "nd")
    filtered_df['damage_hr'] = filtered_df['damage'].map(damage_dict)
    filtered_df['damage_hr'] = filtered_df['damage_hr'].fillna('No Damage')
    return


@app.cell
def _(filtered_df, get_palette, plt, sns, ticker):
    _fig = plt.figure(figsize=(10,5))

    _ax1 = _fig.add_subplot(121)
    _ax2 = _fig.add_subplot(122)

    filtered_df['Damage'] = filtered_df['damage_hr'].apply(lambda x: x.split()[0])
    filtered_df['Assignment'] = filtered_df['Real'].apply(lambda x: "True Positive" if x else "False Positive")

    sns.boxplot(
        data=filtered_df[filtered_df.kmer == 'kmer22'],
        y="FamPercentage",
        x="Assignment",
        order=['True Positive', 'False Positive'],
        ax=_ax1,
        color=get_palette(1)[0]
    )

    _ax1.set_yscale('log')
    _ax1.set_title("A")
    _ax1.set_ylabel("PSF (log)")
    _ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))


    sns.boxplot(
        data=filtered_df[(filtered_df.Real == False)&(filtered_df.kmer == 'kmer22')],
        y="FamPercentage",
        x="dataset",
        hue="Damage",
        hue_order=['No','Medium','High'],
        order=['01','1','10','100','1000'],
        ax=_ax2,
        palette=get_palette(3, r=True)
    )

    _ax2.set_ylabel("PSF")
    _ax2.set_xlabel("Dataset Size")
    _ax2.set_title("B")



    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""## Breadth of Coverage Filter""")
    return


@app.cell
def _(filtered_df, get_palette, plt, sns, ticker):
    _fig = plt.figure(figsize=(10,5))

    _ax1 = _fig.add_subplot(121)
    _ax2 = _fig.add_subplot(122)


    sns.scatterplot(
        data=filtered_df[filtered_df.kmer == 'kmer22'],
        y="Breadth",
        x="Coverage",
        hue="Assignment",
        palette=get_palette(2, r=True),
        hue_order=["True Positive","False Positive"],
        ax=_ax1
    )

    _ax1.set_title("A")
    _ax1.set_xscale("log")
    _ax1.set_ylabel("Breadth of Coverage")
    _ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}x'))


    sns.scatterplot(
        data=filtered_df[filtered_df.kmer == 'kmer22'],
        y="ProportionExpectedBreadth",
        x="Coverage",
        hue="Assignment",
        palette=get_palette(2, r=True),
        hue_order=["True Positive","False Positive"],
        ax=_ax2
    )

    _ax2.set_title("B")
    _ax2.set_xscale("log")
    _ax2.set_ylabel("PEB")
    _ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}x'))


    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""## Combined""")
    return


@app.cell
def _(filtered_df, get_palette, plt, sns, ticker):
    _fig = plt.figure(figsize=(5,5))

    _ax2 = _fig.add_subplot(111)


    _g = sns.scatterplot(
        data=filtered_df[filtered_df.kmer == 'kmer22'],
        y="ProportionExpectedBreadth",
        x="FamPercentage",
        hue="Assignment",
        palette=get_palette(2, r=True),
        hue_order=["True Positive","False Positive"],
        ax=_ax2
    )

    _ax2.set_xscale("log")
    _ax2.set_ylabel("Proportion of expected Breadth")
    _ax2.set_xlabel("Percentage of Sequences")
    _ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))

    _ax2.axhline(y=0.6, color='red', lw=.5, ls='-.')
    _ax2.axvline(x=0.5, color='red', lw=.5, ls='-.')

    plt.tight_layout()
    plt.show()
    return


if __name__ == "__main__":
    app.run()
