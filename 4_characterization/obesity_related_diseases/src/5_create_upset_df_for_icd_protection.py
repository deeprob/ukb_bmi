from functools import reduce
from itertools import combinations
import pandas as pd
import json


def get_sorted_combos(ser, ncombo=3):
    return ";".join(sorted([i for i in ser.loc[[f"Item_{i}_symbol" for i in range(1, ncombo + 1)]] if not pd.isnull(i)]))

def get_combo_set(meta_df, group):
    group_df = meta_df.loc[meta_df.Phenotype==group]
    combo_set = set()
    if len(group_df)>0:
        group_df["combo_name"] = group_df.apply(get_sorted_combos, args=(3, ), axis=1)
        combo_set = set(group_df.combo_name.to_list())
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

def get_upset_df(meta_df, groups):
    # get the set of gene combinations for each group
    combo_dict = {g: get_combo_set(meta_df, g) for g in groups}
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


# protective
protective_combo_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/merged_combos/protective_combos.csv"
protective_df = pd.read_csv(protective_combo_file)
protective_df.Phenotype = protective_df.Phenotype.astype(str)
# # obesity risk 
# combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
# combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
# combo2_df = pd.read_csv(combo2_file)
# combo3_df = pd.read_csv(combo3_file)
# combo_df = pd.concat([combo2_df, combo3_df])
# meta_df = pd.concat([combo_df, protective_df])

icd_names = ['I10', 'E780', 'R074', 'I251', 'I259', 'E039', 'E11', 'BlockM15-M19', 'K80', 'K81', 'K82', 'F32', 'F33', 'G30']
group_names = icd_names #+ ["high_bmi"]
upset_df, unique_intersects = get_upset_df(protective_df, group_names)
upset_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/upset_files/protection.csv"
unique_intersects_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/upset_files/protection_intersects.json"
upset_df.to_csv(upset_file)

unique_intersects = {f"{','.join(k[0])}||{','.join(k[1])}":"||".join(list(v)) for k,v in unique_intersects.items() if ((len(v)>0) and (len(k[0])>1))}
# Serializing json
json_object = json.dumps(unique_intersects, indent=4)
 
# Writing to sample.json
with open(unique_intersects_file, "w") as outfile:
    outfile.write(json_object)
