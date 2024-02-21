import os
import pandas as pd
import argparse
from functools import reduce
from itertools import combinations

import utils.plotting as utpp

def get_combo_set(combo_dir, group):
    combo2_file = os.path.join(combo_dir, group, "combo2.csv")
    combo3_file = os.path.join(combo_dir, group, "combo3.csv")
    try:
        combo2_df = pd.read_csv(combo2_file, usecols=["uniq_items"])
        combo2_set = set(combo2_df.uniq_items.to_list())
    except FileNotFoundError:
        combo2_set = set()
    try:
        combo3_df = pd.read_csv(combo3_file, usecols=["uniq_items"])
        combo3_set = set(combo3_df.uniq_items.to_list())
    except FileNotFoundError:
        combo3_set =set()
    combo_set = combo2_set.union(combo3_set)
    return combo_set

def set_add(a, b):
    return a.intersection(b)

def set_subtract(a,b):
    return a.difference(b)

def get_intersects(include_sets, exclude_sets):
    intersect_sets = reduce(set_add, include_sets)
    if len(exclude_sets)>0:
        exclude_sets = [intersect_sets] + exclude_sets
        intersect_sets = reduce(set_subtract, exclude_sets)
    return intersect_sets

def get_upset_df(combo_dir, groups):
    # get the set of gene combinations for each group
    combo_dict = {g: get_combo_set(combo_dir, g) for g in groups}
    combo_boolean_dict = {g:[] for g in groups}
    counts = []
    unique_intersect_sets = dict()
    for i in range(1, len(groups) + 1):
        # get combinations of length i
        all_combos = list(combinations(groups, i))
        for combos in all_combos:
            # get counts of elements that are unique to the combinations
            include_combos = combos
            exclude_combos = tuple(g for g in groups if g not in include_combos)
            include_sets = [combo_dict[c] for c in include_combos]
            exclude_sets = [combo_dict[c] for c in exclude_combos]
            unique_intersects = get_intersects(include_sets, exclude_sets)
            counts.append(len(unique_intersects))
            unique_intersect_sets[(include_combos, exclude_combos)] = unique_intersects
            for c in include_combos:
                combo_boolean_dict[c].append(True)
            for c in exclude_combos:
                combo_boolean_dict[c].append(False)
    combo_boolean_dict["counts"] = counts
    df = pd.DataFrame(combo_boolean_dict)
    return df.set_index(groups), unique_intersect_sets

def get_unique_intersect_df(unique_intersects):
    uidfs = []
    for k in unique_intersects.keys():
        if (len(k[0])>1)&(len(unique_intersects[k])>0):
            uidf = pd.DataFrame(unique_intersects[k], columns=["uniq_items"])
            uidf["groups"] = "_w_".join(k[0])
            uidfs.append(uidf)

    unique_intersects_filtered_df = pd.concat(uidfs)
    return unique_intersects_filtered_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combo info")
    parser.add_argument("--combo_dir", type=str, help="Filepath of the combo files dir given by Rarecomb")
    parser.add_argument("--groups", type=str, help="Groups ran for combo files", nargs="+")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    upset_df, unique_intersects = get_upset_df(cli_args.combo_dir, cli_args.groups)
    parsed_upset_df = upset_df.loc[upset_df.counts>0]
    fig = utpp.get_upset_plot(parsed_upset_df)
    save_file = os.path.join(cli_args.save_dir, "overlap_combos.pdf")
    utpp.save_pdf(save_file, fig)

    uidf = get_unique_intersect_df(unique_intersects)
    save_file = os.path.join(cli_args.save_dir, "overlap_combos.csv")
    uidf.to_csv(save_file, index=False)
