

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Verification of quicksand using real data / SI 7: Denisova Cave Dataset

        We verified the performance of quicksand v2.3 and the kmer-size of 22 using a dataset published in Zavala et al. (2021).
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import seaborn as sns
    import re
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.patches as mpatches
    from matplotlib_venn import venn2, venn2_circles
    import requests
    from io import StringIO
    from pathlib import Path
    import json
    import numpy as np
    return Path, mo, mpatches, pd, plt, re, sns, ticker, venn2, venn2_circles


@app.cell
def _():
    from styles import set_style, get_palette

    set_style()
    return (get_palette,)


@app.cell
def _(mo):
    mo.md(r"""## Import the Zavala et al. (2021) dataset""")
    return


@app.cell
def _(pd):
    # import the data
    data= pd.read_csv("assets/SI7/20250429_Pleistocene sediment DNA reveals hominin and faunal turnovers at Denisova Cave_library_m_1.csv")
    return (data,)


@app.cell
def _(pd):
    # import family annotations
    familydata = pd.read_csv(   
        "assets/SI7/DENI_family_annotation.tsv",
        sep="\t"
    )

    familydata = { k:v for k,v in zip(familydata["Family"], familydata["Verdict"]) }
    return (familydata,)


@app.cell
def _(Path, data, familydata, pd):
    # fetch the quicksand results
    df = pd.DataFrame()

    for _i,_dat in data[data['Capture Probe']=='AA75'].groupby(['Sequencing Run','Sequencing Lane']):
        _run, _lane = _i
        _p = Path('quicksand_reports/')
        for _run in _p.glob(f"{_run}.{_lane}*.tsv"):
            df = pd.concat(
                [df, pd.read_csv(_run, sep='\t')],
                ignore_index=True
            )

    merged = df.merge(data, left_on='RG', right_on='Capture', how='left', validate='m:1' )
    aa75 = merged[merged['Capture Probe']=='AA75'].copy()

    aa75["Verdict"] = aa75.Family.map(familydata)
    aa75["Fauna"] = aa75["Verdict"].apply(lambda x: x if x != 'Likely' else 'Evidence')


    aa75_ancient = aa75[aa75.Ancientness == '++'].copy()
    return aa75, aa75_ancient


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Remove False Positives

        check the effect of the 0.5 percent and the 0.5 breadth filter
        """
    )
    return


@app.cell
def _(aa75, aa75_ancient, get_palette, plt, sns, ticker):
    _fig = plt.figure(figsize=(10,5))

    _ax1 = _fig.add_subplot(121)
    _ax2 = _fig.add_subplot(122)

    sns.scatterplot(
        data=aa75,
        y="ProportionExpectedBreadth",
        x="FamPercentage",
        ax=_ax1,
        hue='Fauna',
        hue_order=['Evidence','Possible','Unlikely'],
        palette=get_palette(3, r=True)
    )

    _ax1.set_xscale('log')
    _ax1.set_ylabel("PEB")
    _ax1.set_xlabel("PSF")
    _ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))
    _ax1.set_title('A')

    _ax1.axhline(y=0.5, color='red', lw=.5, ls='-.')
    _ax1.axvline(x=0.5, color='red', lw=.5, ls='-.')

    sns.scatterplot(
        data=aa75_ancient,
        y="ProportionExpectedBreadth",
        x="FamPercentage",
        ax=_ax2,
        hue='Fauna',
        hue_order=['Evidence','Possible','Unlikely'],
        palette=get_palette(3, r=True)
    )
    _ax2.set_xscale('log')
    _ax2.set_ylabel("PEB")
    _ax2.set_xlabel("PSF")
    _ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))
    _ax2.set_title('B')

    _ax2.axhline(y=0.5, color='red', lw=.5, ls='-.')
    _ax2.axvline(x=0.5, color='red', lw=.5, ls='-.')

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(aa75_ancient, get_palette, plt, sns):
    _fig = plt.figure(figsize=(13,10))

    _ax1 = _fig.add_subplot(211)
    _ax2 = _fig.add_subplot(212, sharex=_ax1)

    sns.boxplot(
        data=aa75_ancient.sort_values('FamPercentage', ascending=False),
        y="FamPercentage",
        x="Family",
        ax=_ax1,
        hue='Fauna',
        hue_order=['Evidence','Possible','Unlikely'],
        palette=get_palette(3, r=True)
    )

    _ax1.set_ylabel("PSF")
    _ax1.set_xlabel("")
    _ax1.set_yscale('log')
    _ax1.axhline(0.5, ls='-.', color=get_palette(1)[0])

    sns.boxplot(
        data=aa75_ancient.sort_values('FamPercentage', ascending=False)[aa75_ancient.FamPercentage > 0.5],
        y="ProportionExpectedBreadth",
        x="Family",
        ax=_ax2,
        hue='Fauna',
        hue_order=['Evidence','Possible','Unlikely'],
        palette=get_palette(3, r=True)
    )

    _ax2.set_ylabel("PEB")
    _ax2.axhline(0.5, ls='-.', color=get_palette(1)[0])
    _ax2.set_xlabel("")


    plt.setp(_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.setp(_ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")


    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(aa75_ancient):
    aa75_filtered = aa75_ancient[
        (aa75_ancient.FamPercentage > 0.5)
        &(aa75_ancient.ProportionExpectedBreadth > 0.5)
    ].copy()
    return (aa75_filtered,)


@app.cell
def _(mo):
    mo.md(r"""## Plot the Fauna on the Profile""")
    return


@app.cell
def _(aa75_filtered, pd):
    #The sheet is a table representing the coordinates in x and y **cells**
    # like this:
    #    x1  x2  x3  x4
    # Y1 NA  NA  NA M2
    # Y2 NA  M24 NA M3
    # Y3 NA  M23 NA M4
    # Y4 NA  M22 M1 M5

    coord = pd.read_csv("assets/SI7/MAIN_SE points_markers.csv", sep=";")
    coord['y']=coord.index
    coord.columns = [x.replace('0.','') for x in coord.columns]

    # What I want is simpy
    # Marker -> X -> Y
    # So I need to melt/depivot the dataframe

    coord = coord.melt(id_vars=['y'], value_vars=coord.columns)
    coord = coord[(coord.value != None)&(coord.value != 0)]
    coord = coord[coord.value.str.startswith("M")]
    coord.columns = ['y','x','marker']
    coord.marker = coord.marker.apply(lambda x: x.strip())

    #convert dtypes
    coord['x'] = coord.x.astype(int)
    coord['y'] = coord.y.astype(int)

    #find the relative points of coordx and coordy
    #this is required to plot them later on top of the image
    coord['relx'] = coord.x.apply(lambda x: x/max(coord['x']))
    coord['rely'] = coord.y.apply(lambda x: x/max(coord['y']))


    aa75_filtered['Marker'] = aa75_filtered['Sample Synonyms'].apply(lambda x: x.split(';')[1].split(":")[1] if x==x else '')

    samples = aa75_filtered[['Sample Name','Marker']].drop_duplicates(subset=['Sample Name','Marker'])
    samples = {k:v for k,v in zip(samples['Marker'], samples['Sample Name'])}
    coord['SampleID'] = coord['marker'].map(samples)
    return (coord,)


@app.cell
def _(aa75_filtered, coord, mpatches, plt):
    #use the colors from the Zavala2021 extended figure 8
    _colors_dict = {
        "Bovidae":"#f5deb3",
        "Canidae":"#cdcd00",
        "Cervidae":"#87cefa",
        "Elephantidae":"#ffc1c1",
        "Equidae":"#cd1076",
        "Felidae":"#ff4500",
        "Hominidae":"#68228b",
        "Hyaenidae":"#8b4513",
        "Rhinocerotidae":"#66cdaa",
        "Ursidae":"#006400",
        "Other":"grey"
    }

    #and create handles for the legend
    _handles = []
    for _label,_color in _colors_dict.items():
        _patch = mpatches.Patch(color=_color, label=_label)
        _handles.append(_patch)


    _fig = plt.figure(figsize=(9, 16))
    #main frame
    _back = _fig.add_axes([0,0,1,1])

    #print the background
    _img = plt.imread("assets/SI7/MAIN_SE_without column_white background.png")
    _back.imshow(_img)

    plt.legend(
        handles=_handles, 
        prop={'size':15}, 
        facecolor='white', 
        edgecolor='black',
        bbox_to_anchor=(0,0,0.95,0.9))

    #add the points
    for _n,(_i,_grp) in enumerate(coord.groupby('marker')):
        _x = set(_grp['relx']).pop() *0.85 +0.05 #this is to fit it to the padding of the image
        _y = 1-set(_grp['rely']).pop() *0.88 - 0.08 #same here

        _data=aa75_filtered[aa75_filtered.Marker==_i].copy()
        _colors = [_colors_dict[_x] if _x in _colors_dict else "grey" for _x in _data["Family"]]

        if _x==_x and _y==_y:
            _ax = _fig.add_axes([_x,_y,0.04,0.02])
            _sizes = _data['ReadsDeduped_x'].apply(lambda x: int(x))
            _ax.pie([1], colors=["white"], wedgeprops={"edgecolor":"k"}) #plot an empty white pie 
            _ax.pie(
                _sizes, 
                colors=_colors, 
                startangle=-270, 
                counterclock=False, 
                wedgeprops={"edgecolor":"k", "linewidth":0.5}
            )

    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""## Compare the number of reads to BLAST/MEGAN""")
    return


@app.cell
def _(pd):
    megan_aa75 = pd.read_excel("assets/SI7/zavala2021_supplements.xlsx", sheet_name="Mammalian MtDNA Summary", skiprows=2)
    megan_aa75.drop(["Sequences ≥ 35 bp",
                     "Duplication rate",
                     "Number of sequences assigned to a taxa",
                     "Average fragment size of sequences assigned to ancient taxa"
                    ], inplace=True, axis=1)

    megan_aa75 = megan_aa75.melt(
        id_vars=[
            "Field sample ID",
            "MPI-EVA Sample ID",
        ],
        value_vars=[x for x in megan_aa75.columns if x.endswith("idae")]
    ).copy()

    #remove NA
    megan_aa75 = megan_aa75[megan_aa75.value == megan_aa75.value ].copy()

    #filter for main chamber
    megan_aa75 = megan_aa75[megan_aa75["Field sample ID"].str.startswith("M")]

    #manage the columns for merging
    megan_aa75['origin'] = 'Zavala2021'
    megan_aa75.columns = ["marker", "SampleID", "Family", "ReadsDeduped","origin"]
    return (megan_aa75,)


@app.cell
def _(aa75_filtered, megan_aa75):
    qsaa75 = aa75_filtered[["Sample Name","Family","ReadsDeduped_x"]].copy()
    qsaa75.columns = ["SampleID","Family","ReadsDeduped"]
    qsaa75 = qsaa75[(qsaa75.Family.isin(megan_aa75['Family']))]

    #manage columns for merging
    qsaa75['origin'] = 'quicksand v2.3'
    return (qsaa75,)


@app.cell
def _(get_palette, megan_aa75, pd, plt, qsaa75, sns, ticker):
    # merge the datasets
    _fig = plt.figure(figsize=(6,5))
    _ax1 = _fig.add_subplot(111)

    _comp = pd.concat(
        [megan_aa75, qsaa75], 
        ignore_index=True
    )
    _comp.sort_values('ReadsDeduped', ascending=False, inplace=True)

    sns.boxplot(
        data=_comp, 
        x="Family", 
        y="ReadsDeduped", 
        hue="origin",
        hue_order=['Zavala2021','quicksand v2.3'],
        palette=get_palette(2),
        ax=_ax1,
        width=.7
    )

    _ax1.set_ylabel("Number of Sequences")
    _ax1.set_xlabel("")
    _ax1.set_yscale('log')
    _ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))

    plt.setp(_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor");

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(megan_aa75, pd, qsaa75, re):
    _trans = pd.read_csv("assets/SI7/20250429_Pleistocene sediment DNA reveals hominin and faunal turnovers at Denisova Cave_library_m_1.csv")

    _dict = {k:(re.search('SP[0-9]+',v).group() if v==v else v) for k,v in zip(_trans['Sample Name'], _trans['Sample Synonyms'])}
    # get the sample sample-ids
    _qsaa75 = qsaa75.copy()
    _qsaa75['SampleID'] = _qsaa75['SampleID'].map(_dict)

    _comp = _qsaa75.merge(megan_aa75, on=['SampleID','Family'], how='outer')

    #additional families in quicksand
    # The counter shows how often a family was found additionally by quicksand
    from collections import Counter
    print(Counter(_comp[_comp.origin_y != _comp.origin_y]['Family']))
    print(len(set(_comp[_comp.origin_y != _comp.origin_y]['SampleID'])))

    # and the other way around
    print(Counter(_comp[_comp.origin_x != _comp.origin_x]['Family']))
    print(len(set(_comp[_comp.origin_x != _comp.origin_x]['SampleID'])))

    #and positive in both analyses
    print(Counter(_comp[(_comp.origin_x==_comp.origin_x) & (_comp.origin_y==_comp.origin_y)]['Family']))
    print(len(set(_comp[(_comp.origin_x==_comp.origin_x) & (_comp.origin_y==_comp.origin_y)]['SampleID'])))

    _comp[_comp.origin_y == _comp.origin_y]
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Onwards to the human results

        1. Parse the table from Zavalas et al. (2021) supplement to get the "initial screening" data
        2. Get the "screening" tag from the hellforge data
        3. Get the correct capture library by IndexLibID and raw read count (The Capture ID is not reported)
        """
    )
    return


@app.cell
def _(data):
    aa163 = data[data['Capture Probe']=='AA163']
    return (aa163,)


@app.cell
def _(aa163, pd):
    zav2021 = pd.read_csv('assets/SI7/zavala2021_all_huMT_libraries.csv')
    zav2021['Marker'] = zav2021['Marker'].apply(lambda x: x.strip())

    #make sure that we look at the same libraries
    zav2021 = zav2021[zav2021.Library.isin(set(aa163['Library']))]

    def compare_lib_counts(library, raw):
        tmp = aa163[(aa163.Library ==library)&(aa163.ReadsRaw==raw)]
        #print(library, raw, len(tmp)==1, set(aa163[(aa163.Library ==library)]['ReadsRaw']))
        return len(tmp)==1

    zav2021 = zav2021[zav2021.apply(lambda x: compare_lib_counts(x['Library'],x['ReadsRaw']), axis=1)].copy()
    return (zav2021,)


@app.cell
def _(aa163, zav2021):
    # now filter the hellforge set again
    aa163_filtered = aa163[aa163.Library.isin(set(zav2021.Library))].copy()
    return (aa163_filtered,)


@app.cell
def _():
    MIN_SUPPORT = 10 #minimum % support per lineage
    TYPE = 'support' #deam for deaminated sequences or support for all fragments
    return MIN_SUPPORT, TYPE


@app.function
# reduce the datasheets and de-pivot

def convert_dataframe(df, style='quicksand', type='deam'):

    cols = ['Library','Sample Name','Sample Synonyms','MM_Ancient','H_support_deam', f'Deaminated(term3)']

    aa163 = df[[x for x in df.columns if x in cols or x.endswith(type) ]]

    aa163 = aa163.melt(id_vars=cols, value_vars=[x for x in aa163.columns if x.endswith(type)])

    # filter for ancient and significant
    if style=='quicksand':
        aa163 = aa163[aa163.MM_Ancient == '++'].copy()
        if type in ['deam', 'support']:
            aa163 = aa163[aa163.value.str.startswith('^^^')].copy()
    else:
        aa163 = aa163[aa163.MM_Ancient == 'Yes'].copy()
        if type in ['deam', 'support']:
            aa163 = aa163[aa163.value.str.startswith('**')].copy()

    # keep ancient human only if _deam is significant too
    aa163 = aa163[
        (aa163.apply(
            lambda x: x['variable'].startswith('H_') and (
                x['H_support_deam'].startswith("^^^") or x['H_support_deam'].startswith("**")
            ),axis=1)
        ) 
        | (aa163['variable'].str.startswith('H_') == False)
    ]

    return aa163


@app.cell
def _(MIN_SUPPORT, coord, mpatches, plt, re):
    def print_map(df):
    #this is for the legend
        _names_dict = {
             'D-S_support_deam':'Denisova-Sima',
             'D_support_deam':'Denisova',
             'H_support_deam':'Modern Human',
             'N-HST_support_deam':'Neanderthal',
             'N_support_deam':'Neanderthal',
             'D-S_support':'Denisova-Sima',
             'D_support':'Denisova',
             'H_support':'Modern Human',
             'N-HST_support':'Neanderthal',
             'N_support':'Neanderthal',
        }

        _colors_dict = {
             'Denisova-Sima':'pink',
             'Denisova':'red',
             'Modern Human':'yellow',
             'Neanderthal':'blue',
        }

        df['variable'] = df['variable'].map(_names_dict)

        _handles = []
        for _label,_color in _colors_dict.items():
            _patch = mpatches.Patch(color=_color, label=_label)
            _handles.append(_patch)

        # and now the figure
        _fig = plt.figure(figsize=(9, 16))
        #main frame
        _back = _fig.add_axes([0,0,1,1])

        #print the background
        _img = plt.imread("assets/SI7/MAIN_SE_without column_blue and orange background.png")
        _back.imshow(_img)

        plt.legend(
            handles=_handles, 
            prop={'size':17}, 
            facecolor='white', 
            edgecolor='black',
            bbox_to_anchor=(0,0,0.95,0.9))

        #add the points
        for _n,(_i,_grp) in enumerate(coord.groupby('marker')):
            _x = set(_grp['relx']).pop() *0.85 +0.05 #this is to fit it to the padding of the image
            _y = 1-set(_grp['rely']).pop() *0.88 - 0.07 #
            _spid = set(_grp['SampleID']).pop()

            _data=df[df.Marker==_i].copy()

            if _x==_x and _y==_y:        
                _wedgeprops = {"edgecolor":"k", "linewidth":1}

                _ax = _fig.add_axes([_x,_y,0.025,0.02])
                _ax.pie([1],colors=["white"], wedgeprops=_wedgeprops) #plot an empty white pie for negative

                #now check if 
                _data['support'] = _data['value'].apply(lambda x: int(re.search("[0-9]+",x).group()))
                _data['sizes'] = _data['value'].apply(lambda x: int(re.search("(?<=\()[0-9]+",x).group()))

                _data = _data[_data.support > MIN_SUPPORT].copy()        
                _data = _data.groupby("variable",as_index=False).mean(numeric_only=True) # check how elena did it. Take any positive?

                if 'sizes' in _data.columns:
                    _colors = [_colors_dict[x] if x in _colors_dict else "grey" for x in _data["variable"]]
                    _ax.pie(_data['sizes'], colors=_colors, wedgeprops=_wedgeprops)

        plt.show()
    return (print_map,)


@app.cell
def _(TYPE, aa163, coord, print_map):
    aa163_ancient = convert_dataframe(aa163, type=TYPE)
    aa163_ancient['Marker'] = aa163_ancient['Sample Synonyms'].apply(lambda x: x.split(';')[1].split(":")[1] if x==x else '')

    aa163_ancient = aa163_ancient.merge(coord, left_on='Marker', right_on='marker', how='left', validate='m:1')

    print_map(aa163_ancient)
    return


@app.cell
def _(mo):
    mo.md(r"""Here we see an overrepresentation of the Denisova-Sima shared lineage""")
    return


@app.cell
def _(TYPE, coord, print_map, zav2021):
    zav2021['Sample Synonyms'] = zav2021['Marker']
    zav2021['Deaminated(term3)'] = zav2021['ReadsDeam']
    zav2021['MM_Ancient'] = zav2021['Ancient']

    zav2021_ancient = convert_dataframe(zav2021, style='MEGAN', type=TYPE)

    zav2021_ancient['Marker'] = zav2021_ancient['Sample Synonyms']
    zav2021_ancient = zav2021_ancient.merge(coord, left_on='Marker', right_on='marker', how='left', validate='m:1')

    print_map(zav2021_ancient)
    return


@app.cell
def _(mo):
    mo.md(r"""## Compare the number of reads""")
    return


@app.cell
def _(aa163_filtered, pd, zav2021):
    zav2021['Origin'] = 'Zavala2021'
    zav2021['ReadsDeduped'] = zav2021['UniqueHominin']
    aa163_filtered['Origin']='quicksand v2.3'

    hum_combined = pd.concat(
        [zav2021,
         aa163_filtered
        ],
        ignore_index=True
    )

    hum_combined['Family'] = 'Hominidae'
    return (hum_combined,)


@app.cell
def _(get_palette, hum_combined, megan_aa75, pd, plt, qsaa75, sns, ticker):
    # This is the AA75 data
    _fig = plt.figure(figsize=(10,5))
    _ax1 = _fig.add_subplot(1,4,(1,3))
    _ax2 = _fig.add_subplot(144, sharey=_ax1)

    _comp = pd.concat(
        [megan_aa75, qsaa75], 
        ignore_index=True
    )

    _comp['Origin'] = _comp.origin

    _comp.sort_values('ReadsDeduped', ascending=False, inplace=True)

    sns.boxplot(
        data=_comp, 
        x="Family", 
        y="ReadsDeduped", 
        hue="Origin",
        hue_order=['Zavala2021','quicksand v2.3'],
        palette=get_palette(2),
        ax=_ax1,
        width=.7,
    )

    #_ax1.legend(bbox_to_anchor=(1.05, 1.05))
    _ax1.set_ylabel("Number of Sequences")
    _ax1.set_xlabel("")
    _ax1.set_yscale('log')
    _ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))
    _ax1.set_title('A')

    plt.setp(_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor");

    ## and the human reads

    sns.boxplot(
        data=hum_combined[hum_combined['ReadsDeduped']>0], 
        x='Family', 
        y='ReadsDeduped', 
        hue='Origin', 
        ax=_ax2,
        palette=get_palette(2),
        legend=False,
        width=.3
    )

    _ax2.set_title('B')
    _ax2.set_ylabel("")
    _ax2.set_xlabel("")
    _ax2.set_yscale('log')
    #_ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:g}'))
    plt.setp(_ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor");


    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(aa163_filtered, get_palette, plt, venn2, venn2_circles, zav2021):
    _count_ancient = aa163_filtered[['Library','MM_Ancient']].merge(zav2021[['Library','Ancient']], on='Library', how='left', validate='1:1').copy()
    _count_ancient['MM_Ancient'] = _count_ancient.MM_Ancient.apply(lambda x: 'Yes' if x=='++' else 'No')

    _ancient_in_zav_AND_quicksand = len(_count_ancient[_count_ancient.apply(lambda x: x['Ancient'] == x['MM_Ancient'], axis=1)]['Library'])
    _ancient_in_zav_ONLY = len(_count_ancient[_count_ancient.apply(lambda x: x['Ancient'] == 'Yes' and x['MM_Ancient']=='No', axis=1)]['Library'])
    _ancient_in_quicksand_ONLY = len(_count_ancient[_count_ancient.apply(lambda x: x['Ancient'] == 'No' and x['MM_Ancient']=='Yes', axis=1)]['Library'])

    _fig=plt.figure(figsize=(5,5))
    _ax1 = _fig.add_subplot(111)

    # Use the venn2 function
    _v = venn2(
        subsets = (_ancient_in_zav_ONLY, _ancient_in_quicksand_ONLY, _ancient_in_zav_AND_quicksand), 
        set_labels = ('Zavala2021', 'quicksand v2.3'),
    )
    _c = venn2_circles(
        subsets = (_ancient_in_zav_ONLY, _ancient_in_quicksand_ONLY, _ancient_in_zav_AND_quicksand), 
        linestyle='solid',
        lw=1
    )

    _v.get_patch_by_id('10').set_color(get_palette(2)[0])
    _v.get_patch_by_id('11').set_color(get_palette(3)[1])
    _v.get_patch_by_id('01').set_color(get_palette(2)[1])

    _v.get_patch_by_id('10').set_alpha(1.0)
    _v.get_patch_by_id('11').set_alpha(1.0)
    _v.get_patch_by_id('01').set_alpha(1.0)

    _v.get_label_by_id('10').set_text('Zavala2021'+"\n"+str(_ancient_in_zav_ONLY))
    _v.get_label_by_id('11').set_text('Same'+"\n"+str(_ancient_in_zav_AND_quicksand))
    _v.get_label_by_id('01').set_text('quicksand v2.3'+"\n"+str(_ancient_in_quicksand_ONLY))

    _v.get_label_by_id('A').set_text('')
    _v.get_label_by_id('B').set_text('')

    plt.tight_layout()
    plt.show()

    print(_count_ancient[_count_ancient.apply(lambda x: x['Ancient'] == 'Yes' and x['MM_Ancient']=='No', axis=1)])
    return


@app.cell
def _(aa163_filtered):
    # whats with the one 'negative' sample?
    # it turned from ++ to +
    aa163_filtered[aa163_filtered.Library == 'Lib.G.783'][['MM_Ancient', "5'CT(95%CI)","3'CT(95%CI)", 'Ancient','AncientTaxa', ]]
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Adjust the diagnostic positions

        there is this one position (2955) that is really overrepresented and kicks up the D-S lineage.
        """
    )
    return


@app.cell
def _(aa163):
    aa163_all = convert_dataframe(aa163, type='support')
    aa163_all_der = convert_dataframe(aa163, type='der')
    return aa163_all, aa163_all_der


@app.cell
def _(aa163_all_der):
    aa163_all_der['positions'] = aa163_all_der['value'].str.split(')')
    return


@app.cell
def _(aa163_all_der):
    aa163_positions = aa163_all_der.explode('positions')
    aa163_positions = aa163_positions[aa163_positions.positions.apply(lambda x: x==x and x!='')].copy()
    aa163_positions['pos'] = aa163_positions['positions'].apply(lambda x: int(x.replace('(','').split(':')[0]))
    aa163_positions['count'] = aa163_positions['positions'].apply(lambda x: int(x.split(':')[1]))
    return (aa163_positions,)


@app.cell
def _(aa163_positions, plt, sns):
    _fig = plt.figure(figsize=(10,5))
    _ax = _fig.add_subplot(111)

    sns.histplot(
        data = aa163_positions.groupby('pos', as_index=False).sum(numeric_only=True), 
        x='count',
        ax=_ax,
    )

    _ax.set_ylabel('Number of Positions')
    _ax.set_xlabel('Coverage')

    plt.show()
    return


@app.cell
def _(mo):
    mo.md(r"""Here we see that one position that is covered more than 2000 times. This is what we need to remove""")
    return


@app.cell
def _(aa163_positions):
    _test = aa163_positions.groupby('pos', as_index=False).sum(numeric_only=True)
    filterpos = list(_test[_test['count'] > 1000]['pos'])
    return (filterpos,)


@app.cell
def _(aa163_positions, filterpos):
    aa163_filteredpos = aa163_positions[aa163_positions.pos.isin(filterpos) == False].copy()
    aa163_filteredpos.drop(['value','positions','Sample Name'], axis=1, inplace=True)

    ## now group back together to find the significant calls
    aa163_filteredpos['Marker'] = aa163_filteredpos['Sample Synonyms'].apply(lambda x: x.split(';')[1].split(":")[1] if x==x else '')
    aa163_filteredpos['position'] = aa163_filteredpos.apply(lambda x: f"({x['pos']}:{x['count']})" ,axis=1)
    aa163_filteredpos = aa163_filteredpos.groupby(['Library','Marker','variable'], as_index=False).sum()
    aa163_filteredpos['n_pos'] = aa163_filteredpos.position.apply(lambda x: sum(1 for y in x if y=='('))
    aa163_filteredpos['variable'] = aa163_filteredpos['variable'].apply(
        lambda x: x.replace('_der','_support').replace('_deam_support','_support_deam')
    )
    aa163_filteredpos = aa163_filteredpos[['Library', 'Marker','variable','position','count', 'n_pos']].copy()
    return (aa163_filteredpos,)


@app.cell
def _(aa163_all, aa163_filteredpos, re):
    def get_sig_marker(x):
        if x['support'] >= 10 and x['n_pos']>=3:
            return "^^^"
        if x['support'] >= 10:
            return "^^"
        if x['support'] >= 5:
            return "^"
        return ""


    aa163_all['Marker'] = aa163_all['Sample Synonyms'].apply(lambda x: x.split(';')[1].split(":")[1] if x==x else '')
    aa163_newposmerged = aa163_all.merge(aa163_filteredpos, on=['Library','Marker','variable'], how='left', validate='1:1')

    aa163_newposmerged['overlapping'] = aa163_newposmerged['value'].apply(lambda x: int(re.search(r'(?<=\/.)[0-9]+', x).group()))
    aa163_newposmerged['initial'] = aa163_newposmerged['value'].apply(lambda x: int(re.search(r'(?<=\()[0-9]+', x).group()))

    # if we kicked out positions, the counts there need to be removed from the number of overlapping positions as well
    # we see that we kicked out positions, if the count (new) is not equal to the 'initial'

    aa163_newposmerged['adjusted_overlapping'] = aa163_newposmerged.apply(
        lambda x: x['overlapping'] if x['initial'] == x['count'] else x['overlapping']-(x['initial']-x['count']), axis=1
    )
    aa163_newposmerged['support'] = aa163_newposmerged.apply(lambda x: 100 * round(x['count']/x['adjusted_overlapping'], 2), axis=1)
    aa163_newposmerged['sig'] = aa163_newposmerged.apply(lambda x: get_sig_marker(x) , axis=1)

    ## now assemble new 'values'
    aa163_newposmerged['value'] = aa163_newposmerged.apply(
        lambda x:
        f"{x['sig']}{x['support']}% ({x['count']} / {x['adjusted_overlapping']})"       
        ,axis=1
    )
    return (aa163_newposmerged,)


@app.cell
def _(aa163_newposmerged, print_map):
    print_map(aa163_newposmerged[aa163_newposmerged['sig']=='^^^'])
    return


if __name__ == "__main__":
    app.run()
