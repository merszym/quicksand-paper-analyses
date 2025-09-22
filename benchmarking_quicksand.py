

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from pathlib import Path
    import pysam
    import re
    return Path, pd, plt, re, sns, ticker


@app.cell
def _(mo):
    mo.md(
        r"""
        What is the plan?

        1. Load all the run-information for megan, quicksand and euka
        2. make some basic plots about the accuracy and sensitivity
        """
    )
    return


@app.cell
def _(plt):
    from styles import get_palette

    plt.style.use("fast")
    return (get_palette,)


@app.cell
def _(mo):
    mo.md(r"""# Accuracy""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Parse the Data""")
    return


@app.cell
def _():
    family_dict = {
        "Bovidae": "DQ124371.1",
        "Suidae": "KC469586.1",
        "Hyaenidae": "JF894377.1",
        "Elephantidae": "EU153448.1",
        "Hominidae": "KC879692.1",
    }
    acc_dict = {v: k for k, v in family_dict.items()}

    damage_dict = {"md": "Medium Damage", "hd": "High Damage"}
    return damage_dict, family_dict


@app.cell
def _(mo):
    mo.md(
        r"""
        ### Parse MEGAN + BLAST

        1. We need to reduce the dataset to 1 species per family and also calculate the FamPercentage
        2. Filter as usually (1% threshold)
        """
    )
    return


@app.cell
def _(Path, pd):
    # load the megan sim-data
    megan = pd.DataFrame()

    for _report in Path("assets/BENCH/blastmegan_nf/").glob("*/final_report.tsv"):
        _tmp = pd.read_csv(_report, sep="\t")
        megan = pd.concat([megan, _tmp], ignore_index=True)
    return (megan,)


@app.cell
def _(damage_dict, megan, re):
    # calculate the percentages and reduce the dataset to the highest species per family


    def megan_get_max(rg, fam, df):
        _tmp = df[(df.Family == fam) & (df.RG == rg)].copy()
        return max(_tmp["ReadsDeduped"])


    def get_fam_percentage(rg, fam, df):
        _tmp = df[df.RG == rg].copy()
        famreads = sum(
            _tmp[_tmp.Family == fam]["ReadsDeduped"]
        )  # its only one, but I prefer sum over iloc
        return famreads / sum(_tmp["ReadsDeduped"])


    # reduce to 1 species per family
    megan["damage"] = megan.RG.apply(
        lambda x: re.search("(md|hd)", x).group()
        if bool(re.search("(md|hd)", x))
        else "nd"
    )
    megan["damage_hr"] = megan["damage"].map(damage_dict)

    megan["MaxPerFamily"] = megan.apply(
        lambda x: megan_get_max(x["RG"], x["Family"], megan), axis=1
    )
    # calculate the percentage
    megan_v2 = (
        megan[megan.MaxPerFamily == megan.ReadsDeduped]
        .copy()
        .groupby(["RG", "Family"], as_index=False)
        .first()
    )
    megan_v2["FamPercentage"] = megan_v2.apply(
        lambda x: get_fam_percentage(x["RG"], x["Family"], megan_v2), axis=1
    )
    return (megan_v2,)


@app.cell
def _(pd):
    # read in each bamfile to get the number of correct reads
    def compare_readwise_megan(df):
        """
        In the df should be the columns RG, Family. The function returns another
        df with the summary of the deduped family bam file

        looks per read of a family, if the sequence comes from that family or not, then
        calculates a false-positive rate.
        """

        readwise = pd.read_csv("assets/BENCH/blastmegan_nf/MEGAN_readwise.csv")

        summ = readwise.groupby(["RG", "Family"], as_index=False).sum()
        summ["FP_Rate"] = summ.False_Positive / summ.All
        summ["Accuracy"] = 1 - summ.FP_Rate

        return summ
    return (compare_readwise_megan,)


@app.cell
def _(compare_readwise_megan, megan_v2):
    summ = compare_readwise_megan(megan_v2)
    return (summ,)


@app.cell
def _(megan_v2, summ):
    megan_v3 = megan_v2.merge(summ, on=["RG", "Family"], how="left", validate="1:1")
    return (megan_v3,)


@app.cell
def _(megan_v3):
    # summarize the data for box-plotting -> get accuracy value per dataset and damage profile

    megan_v3["dataset"] = megan_v3.RG.apply(
        lambda x: "1" if not "1" in x else x.replace(".fq", "").split("_")[-1]
    )
    megan_v3["damage_hr"] = megan_v3["damage_hr"].fillna("No Damage")

    megan_v3_grouped_filtered = (
        megan_v3[megan_v3.FamPercentage > 0.01]
        .groupby(["dataset", "damage_hr"], as_index=False)
        .sum(numeric_only=True)
    )
    megan_v3_grouped_filtered["Accuracy"] = (
        100
        * megan_v3_grouped_filtered["True_Positive"]
        / megan_v3_grouped_filtered["All"]
    )
    return (megan_v3_grouped_filtered,)


@app.cell
def _(megan_v3):
    megan_v2_sens = (
        megan_v3[megan_v3.FamPercentage > 0.01]
        .groupby(["RG", "Family"], as_index=False)
        .sum(numeric_only=True)
    )
    return


@app.cell
def _(mo):
    mo.md(r"""### Parse euka runs""")
    return


@app.cell
def _(Path, pd):
    euka = pd.DataFrame()

    for _report in Path("assets/BENCH/euka_nf/").glob("*/final_report.tsv"):
        _tmp = pd.read_csv(_report, sep="\t")
        euka = pd.concat([euka, _tmp], ignore_index=True)
    return (euka,)


@app.cell
def _(damage_dict, euka, re):
    euka["dataset"] = euka.RG.apply(
        lambda x: "1" if not "1" in x else x.replace(".fq", "").split("_")[-1]
    )
    euka["damage"] = euka.RG.apply(
        lambda x: re.search("(md|hd)", x).group()
        if bool(re.search("(md|hd)", x))
        else "nd"
    )
    euka["damage_hr"] = euka["damage"].map(damage_dict)
    euka["damage_hr"] = euka["damage_hr"].fillna("No Damage")
    return


@app.cell
def _(Path, pd):
    # read in each txt-file to get the number of correct reads per taxon
    def compare_readwise_euka(df):
        """
        In the df should be the columns RG, EukaTaxa. The function returns another
        df with the summary of the deduped family bam file

        looks per read of a family, if the sequence comes from that family or not, then
        calculates a false-positive rate.
        """

        readwise = []
        for i, dat in df.groupby(["dataset"]):
            factor = i[0]
            for rg in set(dat["RG"]):  # the damage profiles
                reads = {}
                p = Path(
                    f"assets/BENCH/euka_nf/{factor}/euka/{rg}/euka_output_FragNames.tsv"
                )
                with open(p) as infile:
                    for line in infile:
                        data = line.replace("\n", "").split("\t")
                        tmp = []
                        taxon = data[0]
                        for read in data[1:]:
                            tmp.append(
                                {
                                    "RG": rg,
                                    "EukaTaxa": taxon,
                                    "Origin": read.split(":")[0],
                                }
                            )
                        readwise.extend(tmp)
        readwise = pd.DataFrame(readwise)

        return readwise
    return (compare_readwise_euka,)


@app.cell
def _(compare_readwise_euka, euka):
    euka_reads = compare_readwise_euka(euka)
    euka_v2 = euka_reads.merge(euka, on=["RG", "EukaTaxa"], how="left", validate="m:1")
    return (euka_v2,)


@app.cell
def _(euka_v2):
    _family_dict = {
        "Bovidae": "DQ124371.1",
        "Suina": "KC469586.1",
        "Carnivora": "JF894377.1",
        "Proboscidea": "EU153448.1",
        "Hominidae": "KC879692.1",
    }
    _acc_dict = {v: k for k, v in _family_dict.items()}

    euka_v2["ExpectedTaxa"] = euka_v2["Origin"].map(_acc_dict)
    euka_v2["True_Positive"] = euka_v2["ExpectedTaxa"] == euka_v2["EukaTaxa"]
    euka_v2["All"] = 1

    # use only the ones, where euka was able to find a family!
    euka_v2_grouped = (
        euka_v2[euka_v2.EukaTaxa == euka_v2.EukaTaxa]
        .groupby(["dataset", "damage_hr"], as_index=False)
        .sum(numeric_only=True)
    )
    euka_v2_grouped["Accuracy"] = (
        100 * euka_v2_grouped["True_Positive"] / euka_v2_grouped["All"]
    )
    return (euka_v2_grouped,)


@app.cell
def _(euka_v2):
    euka_v2_sens = euka_v2.groupby(["RG", "Family"], as_index=False).sum(
        numeric_only=True
    )
    return


@app.cell
def _(mo):
    mo.md(r"""### parse quicksand runs""")
    return


@app.cell
def _(Path, pd):
    quicksand = pd.DataFrame()

    for _report in Path("assets/BENCH/quicksand_v2.3/").glob(
        "*/filtered_report_0.5p_0.5b.tsv"
    ):
        _tmp = pd.read_csv(_report, sep="\t")
        quicksand = pd.concat([quicksand, _tmp], ignore_index=True)
    return (quicksand,)


@app.cell
def _(damage_dict, quicksand, re):
    quicksand["dataset"] = quicksand.RG.apply(
        lambda x: "1" if not "1" in x else x.replace(".fq", "").split("_")[-1]
    )
    quicksand["damage"] = quicksand.RG.apply(
        lambda x: re.search("(md|hd)", x).group()
        if bool(re.search("(md|hd)", x))
        else "nd"
    )
    quicksand["damage_hr"] = quicksand["damage"].map(damage_dict)
    quicksand["damage_hr"] = quicksand["damage_hr"].fillna("No Damage")
    return


@app.cell
def _(pd, quicksand):
    quicksand_v2 = pd.read_csv(
        "assets/BENCH/quicksand_v2.3/quicksand_readwise.tsv", sep="\t"
    ).merge(quicksand, on=["RG", "Family"], how="left", validate="m:1")
    return (quicksand_v2,)


@app.cell
def _(quicksand_v2):
    quicksand_v2_grouped = quicksand_v2.groupby(
        ["dataset", "damage_hr"], as_index=False
    ).sum(numeric_only=True)
    quicksand_v2_grouped["Accuracy"] = (
        100 * quicksand_v2_grouped["True_Positive"] / quicksand_v2_grouped["All"]
    )
    return (quicksand_v2_grouped,)


@app.cell
def _(quicksand):
    quicksand_v2_sens = quicksand.groupby(["RG", "Family"], as_index=False).sum(
        numeric_only=True
    )
    return


@app.cell
def _(mo):
    mo.md(r"""### Parse HOLI runs""")
    return


@app.cell
def _(Path, damage_dict, pd, re):
    holi = pd.DataFrame()

    for _report in Path("assets/BENCH/HOLI_nf/").glob("*/final_report.tsv"):
        _tmp = pd.read_csv(_report, sep="\t")
        holi = pd.concat([holi, _tmp], ignore_index=True)

    holi["RG"] = holi["sample"]

    holi["dataset"] = holi.RG.apply(
        lambda x: "1" if not "1" in x else x.replace(".fq", "").split("_")[-1]
    )
    holi["damage"] = holi.RG.apply(
        lambda x: re.search("(md|hd)", x).group()
        if bool(re.search("(md|hd)", x))
        else "nd"
    )
    holi["damage_hr"] = holi["damage"].map(damage_dict)
    holi["damage_hr"] = holi["damage_hr"].fillna("No Damage")
    return (holi,)


@app.cell
def _():
    import statistics
    return (statistics,)


@app.cell
def _(holi, pd, statistics):
    holi_minreads = 5
    holi_minpercent = 1

    holi_readsfiltered = holi[
        (holi.nreads > holi_minreads)
        &(holi.taxa_path.str.contains('family'))
    ].copy()

    holi_final = pd.DataFrame()

    for _rg in set(holi["RG"]):
        _subset = holi_readsfiltered[holi_readsfiltered.RG == _rg].copy()

        #reduce to family-level (they are summed up, so otherwise the percentage and median is screwed)
        _subset = _subset[_subset['rank']=='family'].copy()
    
        _median_filter = statistics.median(_subset["nreads"]) / 2
        _subset_filtered = _subset[_subset["nreads"] > _median_filter].copy()

        _subset_filtered["TotalSum"] = sum(_subset_filtered["nreads"])
        _subset_filtered["FamPercentage"] = (
            _subset_filtered["nreads"] / _subset_filtered["TotalSum"]
        ) * 100
        _subset_filtered_final = _subset_filtered[_subset_filtered['FamPercentage']>holi_minpercent]

        holi_final = pd.concat([holi_final, _subset_filtered_final], ignore_index=True)

    holi_final["assigned"] = holi_final.taxa_path.str.extract(
        r"([A-Za-z]+idae)"
    )

    #remove 'superfamily' level
    holi_final = holi_final[holi_final.assigned == holi_final.assigned]
    return (holi_final,)


@app.cell
def _(Path, pd, re):
    holi_readwise = pd.DataFrame()

    for _report in Path("assets/BENCH/HOLI_nf/").glob("*/*.gz"):
        _rg = _report.parts[4].split(".")[0]
        _tmp = pd.read_csv(_report, sep="\t", skiprows=1)
        _tmp["RG"] = _rg
        holi_readwise = pd.concat([holi_readwise, _tmp], ignore_index=True)

    holi_readwise["dataset"] = holi_readwise.RG.apply(
        lambda x: "1" if not "1" in x else x.replace(".fq", "").split("_")[-1]
    )
    holi_readwise["damage"] = holi_readwise.RG.apply(
        lambda x: re.search("(md|hd)", x).group()
        if bool(re.search("(md|hd)", x))
        else "nd"
    )
    return (holi_readwise,)


@app.cell
def _(family_dict, holi_readwise):
    _family_dict = {v: k for k, v in family_dict.items()}


    # filter for everything that has a family
    holi_readwise_v2 = holi_readwise[
        holi_readwise.taxa_path.str.contains("family")
    ].copy()
    holi_readwise_v2["origin"] = holi_readwise_v2["queryid"].apply(
        lambda x: x.split(":")[0]
    )
    holi_readwise_v2["origin_family"] = holi_readwise_v2["origin"].map(_family_dict)
    holi_readwise_v2["assigned"] = holi_readwise_v2.taxa_path.str.extract(
        r"([A-Za-z]+idae)"
    )

    holi_readwise_v2['test'] = holi_readwise_v2['origin_family'] == holi_readwise_v2['assigned']
    holi_readwise_v2['sum'] = 1
    return (holi_readwise_v2,)


@app.cell
def _(damage_dict, holi_final, holi_readwise_v2):
    #(now apply the filters)
    _data = holi_final.merge(holi_readwise_v2, on=['dataset','damage','assigned'], how='left', validate='1:m')

    holi_readwise_v2_grouped = _data.groupby(['dataset','damage'],as_index=False).sum(numeric_only=True)
    holi_readwise_v2_grouped['Accuracy'] = ( holi_readwise_v2_grouped['test'] / holi_readwise_v2_grouped['sum'] ) * 100
    holi_readwise_v2_grouped["damage_hr"] = holi_readwise_v2_grouped["damage"].map(damage_dict)
    holi_readwise_v2_grouped["damage_hr"] = holi_readwise_v2_grouped["damage_hr"].fillna("No Damage")
    return (holi_readwise_v2_grouped,)


@app.cell
def _(mo):
    mo.md(r"""# Make Summaries""")
    return


@app.cell
def _(
    euka_v2_grouped,
    holi_readwise_v2_grouped,
    megan_v3_grouped_filtered,
    pd,
    quicksand_v2_grouped,
):
    megan_v3_grouped_filtered["Method"] = "BLAST/MEGAN"
    euka_v2_grouped["Method"] = "euka"
    quicksand_v2_grouped["Method"] = "quicksand"
    holi_readwise_v2_grouped['Method'] = 'HOLI'

    final_summary_accuracy = pd.concat(
        [
            megan_v3_grouped_filtered[
                ["Method", "dataset", "damage_hr", "Accuracy"]
            ].copy(),
            euka_v2_grouped[["Method", "dataset", "damage_hr", "Accuracy"]].copy(),
            quicksand_v2_grouped[["Method", "dataset", "damage_hr", "Accuracy"]].copy(),
            holi_readwise_v2_grouped[['Method','dataset','damage_hr','Accuracy']].copy()
        ],
        ignore_index=True,
    )
    return (final_summary_accuracy,)


@app.cell
def _(pd):
    runtime = pd.read_csv("assets/BENCH/runtime.tsv", sep="\t")
    runtime["seconds"] = runtime.runtime.apply(
        lambda x: float(x)
        if "m" not in x
        else 60 * float(x.split("m")[0]) + float(x.split("m")[1])
    )
    runtime["minutes"] = runtime.seconds / 60
    runtime["size"] = runtime.dataset.apply(lambda x: x.split("_")[-1])
    return (runtime,)


@app.cell
def _(final_summary_accuracy, get_palette, plt, runtime, sns, ticker):
    _fig = plt.figure(figsize=(10, 5))

    _ax1 = _fig.add_subplot(121)
    _ax3 = _fig.add_subplot(122)

    _palette = {
        m: c
        for m, c in zip(
            ["BLAST/MEGAN", "euka", "HOLI", "quicksand"], get_palette(4, r=True)
        )
    }

    final_summary_accuracy["Damage"] = final_summary_accuracy["damage_hr"].apply(
        lambda x: x.split()[0]
    )

    sns.boxplot(
        data=final_summary_accuracy,
        x="dataset",
        y="Accuracy",
        hue="Method",
        # order=["No", "Medium", "High"],
        hue_order=["BLAST/MEGAN", "euka", "HOLI", "quicksand"],
        ax=_ax1,
        palette=_palette,
        legend=False,
        width=0.75,
    )

    _ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    _ax1.set_ylabel("Accuracy (%)", size=15)
    _ax1.set_xlabel("Dataset Size", size=15)
    _ax1.set_title("A", size=21)

    sns.lineplot(
        data=runtime,
        x="size",
        y="minutes",
        hue="Method",
        ax=_ax3,
        palette=_palette,
        linewidth=4,
        hue_order=["BLAST/MEGAN", "euka", "HOLI", "quicksand"],
    )
    sns.despine()

    sns.move_legend(_ax3, "upper left")

    _ax3.set_ylabel("Runtime (Minutes)", size=15)
    _ax3.set_xlabel("Dataset Size", size=15)
    _ax3.set_title("B", size=21)

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""### Presence/Absence""")
    return


@app.cell
def _(euka_v2):
    euka_v2.groupby(["RG", "EukaTaxa"]).count()[["Origin"]].reset_index().pivot_table(
        index=["RG"], columns=["EukaTaxa"], values="Origin"
    )
    return


@app.cell
def _(megan_v3):
    megan_v3[(megan_v3.FamPercentage > 0.01)].groupby(["RG", "Family"]).sum()[
        ["ReadsMapped"]
    ].reset_index().pivot_table(index=["RG"], columns=["Family"], values="ReadsMapped")
    return


@app.cell
def _(quicksand):
    quicksand[
        (quicksand.FamPercentage > 0.5) & (quicksand.ProportionExpectedBreadth > 0.5)
    ].groupby(["RG", "Family"]).sum()[["ReadsMapped"]].reset_index().pivot_table(
        index=["RG"], columns=["Family"], values="ReadsMapped"
    )
    return


@app.cell
def _(holi_final):
    holi_final[holi_final['sample'] == 'simadna_hd_1000.fq']
    return


if __name__ == "__main__":
    app.run()
