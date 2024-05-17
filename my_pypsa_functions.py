
# Functions modified from those at PYPSA: add_electricity.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Annuity formula
def annuity(n,r):
    """Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20,0.05)*20 = 1.6"""

    if r > 0:
        return r/(1. - 1./(1.+r)**n)
    else:
        return 1/n

def calculate_annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and.
    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6
    """

    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def load_tech_parameters(name_file_tech_parameters, config, Nyears=1.0):
    ''' set all asset costs and other parameters
     this a modified version of the function load_costs from PyPSA add_electricity.py
    Modifications:
    -  the gas fuel cost are not inserted in the marginal cost.  This has been done in order to model, for example, a gas turbine as  a link from a multi-source gas bus, where the CH4 gas could be form biogenic,  synthesis, or fossil origin. The respective fuel costs and emissions are accounted before the CH4 bus
    - there is not an aggregation of differet solar plants (utility, commercial, residential). This is avoided in order to make an analysis with different capacity factors for each type pf solar plants. The previous aggregation code is commented
    - N.B.: TO BE DONE land use data added in the csv and in the config files 
    '''


    tech_parameters = pd.read_csv(name_file_tech_parameters, index_col=[0, 1]).sort_index()

    # correct units to MW
    tech_parameters.loc[tech_parameters.unit.str.contains("/kW"), "value"] *= 1e3
    tech_parameters.unit = tech_parameters.unit.str.replace("/kW", "/MW")

    fill_values = config["fill_values"]
    tech_parameters = tech_parameters.value.unstack().fillna(fill_values)
    
    tech_parameters["capital_cost"] = (
        (
            calculate_annuity(tech_parameters["lifetime"], tech_parameters["discount rate"])
            + tech_parameters["FOM"] / 100.0
        )
        * tech_parameters["investment"]
        * Nyears
    )

    tech_parameters["marginal_cost_no_fuel"] = tech_parameters["VOM"] 

    tech_parameters = tech_parameters.rename(columns={"CO2 intensity": "co2_emissions"})

 #   tech_parameters.at["solar", "capital_cost"] = (config["rooftop_residential_share"] * tech_parameters.at["solar-rooftop residential", "capital_cost"]  + config["rooftop_commercial_share"] * tech_parameters.at["solar-rooftop commercial", "capital_cost"] + (1 - config["rooftop_residential_share"]  - config["rooftop_commercial_share"]) * tech_parameters.at["solar-utility", "capital_cost"])   
    # tech_parameters.at["solar", "marginal_cost_no_fuel"] = 0.01
    
    return tech_parameters


def colors_distinguishable_by_color_blind_person(num_of_categories):
    category_colors = plt.colormaps['inferno'](np.linspace(0., 1., num_of_categories))
    return category_colors

def my_stacked_horizontal_bar_plot(results, category_labels):
    """
    Parameters
    ----------
    results : dict
        A key entry is a label that describe the results 
        The list of values associated to a key are the results specifics to a category.
        The category labels are in the input category_labels
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_labels*.
    category_labels : list of str
        The category labels.
    """
    result_labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)    
    category_colors = colors_distinguishable_by_color_blind_person(data.shape[1])
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())
    cat_displayed = 0
    for i, (colname, color) in enumerate(zip(category_labels, category_colors)):
        widths = data[:, i]
        # Avoid displaying zero categories, approx to a very small values
        if widths.sum() > 0.001:
            cat_displayed +=1
            starts = data_cum[:, i] - widths
            rects = ax.barh(result_labels, widths, left=starts, height=0.2,
                        label=colname, color=color)

            r, g, b, _ = color
            text_color = 'white' if r * g * b < 0.2 else 'black'
            ax.bar_label(rects, label_type='center', color=text_color)
        
    ax.legend(ncol=cat_displayed, bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')

    return fig, ax
