from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import re

<<<<<<< HEAD
from staticinteract import StaticInteract, DropDownWidget, RangeWidget
=======
from staticinteract import StaticInteract
from stacticinteract import DropDownWidget, RangeWidget

>>>>>>> d7d6d2fc67508c3f91490fafcd2d5f03681bf40d

class PathwayPlotter:
    """
    Methods  to plot pathway representation results
    """
    def __init(self):
        pass


def extractKoPathwayName(Ko_str):
    return re.sub('\[.*?\]', '', Ko_str).strip()


# Plotting functions
def plotClusterData(pdata, cluster, ax=None, cluster_id=None):
    pdata[pdata.index.isin(cluster)].transpose().plot(
        legend=False, title=f'{cluster_id}, size={len(cluster)}',
        ax=ax, color='#9a9a9a', linewidth=0.8,
        marker='.', markerfacecolor='#ee9929', markersize=12)


def plotKEGGFrequencies(data: dict, color=None, axis=None):
    """
    Bar plot of sorted KEGG systems or subsystems
    data: dictionary wth keys being system or subsystem names and
    values, their frequencies/representation.
    """
    if color is None:
        color = 'C0'
    pvalues = [v[1] for v in data.values()]
    clean_name_data = 100 * pd.Series(
        {extractKoPathwayName(k): data[k][0] for k in data.keys()}
    )
    ax = clean_name_data.plot.bar(figsize=(12, 8), color=color, ax=axis)
    for i, p in enumerate(ax.patches):
        ax.annotate(f'({pvalues[i]:.4f})', (p.get_x() * 1.005, p.get_height() * 1.006))


def plotSystemsAndSubsystemsWebPage(clusters, cluster_data, permutation_results,
                                    plot_first_N=10, color=None,
                                    img_folder_name=None):
    """
    NOTE: Increase tick labels font size!
    """

    plt.rcParams.update({'figure.max_open_warning': 0})
    if color is None:
        color = 'C0'
    if img_folder_name is None:
        img_folder_name = 'iplot'

    system_types = ['system', 'pathway']
    cluster_ids = list(clusters.keys())

    def plot_fun(system_type, cluster_id):

        fig, ax = plt.subplot_mosaic(
            """
            A
            B
            """,
            gridspec_kw = {"height_ratios": [0.7, 1]}
        )

        ax['B'].set_ylabel('Pathway representation (%)')
        ax['B'].set_title(f'{system_type}s (sample p-value)')

        if cluster_id in clusters.keys():
            p_Data = permutation_results
            kdata = {k: v
                     for k,v in p_Data[system_type][cluster_id].items()}
            if len(kdata) > plot_first_N:
                kdata = {k: kdata[k] for k in list(kdata.keys())[:10]}
            plotClusterData(cluster_data, clusters[cluster_id],
                            ax['A'], cluster_id)
            plotKEGGFrequencies(kdata,
                                color=color, axis=ax['B'])
        fig.set_figwidth(20)
        fig.set_figheight(20)
        return fig

    i_fig = StaticInteract(plot_fun,
                           system_type=DropDownWidget(system_types,
                                        description='Level'),
                           cluster_id=DropDownWidget(cluster_ids,
                                        description='Cluster ID'),
                           interact_name=img_folder_name)
    return i_fig


def getRepresentationStackedPlotData(p_paths, pvalue_cutoff=1):
    """
    Make pandas dataframe of pathway representation across
    clusters
    """
    unique_systems = np.unique([
        system for c in p_paths['system'].values()
        for system in c.keys()
    ])
    unique_subsystems = np.unique([
        subsystem
        for c in p_paths['pathway'].values()
        for subsystem in c.keys()
    ])

    rep_systems, rep_subsystems = {}, {}
    for system in unique_systems:
        rep_systems[system] = []
        for cluster_id, values in p_paths['system'].items():
            rep, pvalue = values[system]
            if pvalue < pvalue_cutoff:
                rep_systems[system].append(rep)
            else:
                rep_systems[system].append(0.0)
        # Add rep outside clusters
        rep_systems[system].append(1 - sum(rep_systems[system]))

    for subsystem in unique_subsystems:
        rep_subsystems[subsystem] = []
        for cluster_id, values in p_paths['pathway'].items():
            rep, pvalue = values[subsystem]
            if pvalue < pvalue_cutoff:
                rep_subsystems[subsystem].append(rep)
            else:
                rep_subsystems[subsystem].append(0.0)
        # Add rep outside clusters
        rep_subsystems[subsystem].append(1 - sum(rep_subsystems[subsystem]))

        rep_subsystems = {extractKoPathwayName(k): v
                          for k,v in rep_subsystems.items()}

    return {'system': rep_systems, 'pathway': rep_subsystems}


def plotSystemsAndSubsystemsStacked(permutation_results, cluster_colors, img_folder_name):
    """
    """
    plt.rcParams.update({'figure.max_open_warning': 0})
    if img_folder_name is None:
        img_folder_name = 'iplot'

    # Remove No cluster assigned data
    p_Data = {system_type: {system: value for system, value in data.items() if system != 'No_cluster_assigned'}
              for system_type, data in permutation_results.items()}

    system_types = ['system', 'pathway']

    def plot_fun(system_type, tpvalue_cutoff):
        fig = plt.figure(figsize=(10, 8))
        ax = plt.gca()
        rep_data = getRepresentationStackedPlotData(p_Data, tpvalue_cutoff)

        df = pd.DataFrame(rep_data[system_type],
                          index=list(p_Data['system'].keys()) + ['No cluster'])
        df = df.loc[:, df.loc['No cluster'] < 1.0].sort_values('No cluster', axis=1)
        ax = (100 * df).transpose().plot.bar(stacked=True, legend=True,
                                             figsize=(10,8),
                                             color=list(cluster_colors.values()),
                                             ax=ax)
        return fig

    i_fig = StaticInteract(plot_fun,
                           system_type=DropDownWidget(system_types,
                                        description='Level'),
                           tpvalue_cutoff=RangeWidget(min=0.05, max=1,
                                                     step=0.05,
                                                     default=0.1,
                                                     description='p-value cutoff'),
                           interact_name=img_folder_name)
    return i_fig
