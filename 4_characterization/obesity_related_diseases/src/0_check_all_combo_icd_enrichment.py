import pandas as pd
from statsmodels.stats.multitest import multipletests
import subprocess
import os
from tqdm import tqdm


def run_fishers_enrichment(contingency, table_out, fishers_bash_path):
    cmd = [
        "bash", fishers_bash_path, 
        contingency, table_out
        ]
    subprocess.run(cmd)
    return

def run_fishers_exact_in_r(contingency_table, fishers_bash_path, tmp_dir="/data5/deepro/tmp"):
    # create temporary contingency file
    contingency_file = os.path.join(tmp_dir, "tmp_contingency.csv")
    with open(contingency_file, "w") as f:
        f.write("data_type,condition,control\n")
        f.write(f"overlap,{contingency_table[0][0]},{contingency_table[0][1]}\n")
        f.write(f"nonoverlap,{contingency_table[1][0]},{contingency_table[1][1]}\n")
    results_file = os.path.join(tmp_dir, "tmp_res.csv")
    # run enrichment
    run_fishers_enrichment(contingency_file, results_file, fishers_bash_path)
    df = pd.read_csv(results_file)
    os.remove(contingency_file)
    os.remove(results_file)
    return df.iloc[0].values

def get_samples_with_combo(combo_dfs, ncombos, wes_file, selected_samples_file):
    all_enriched_combos = []
    all_cols = []
    for combo_df, nc in zip(combo_dfs, ncombos):
        # get unique items that came out of the enriched combinations
        enriched_combos = combo_df.loc[:, [f"Item_{c}" for c in range(1, nc+1)]].values
        columns_to_use = sorted(set(enriched_combos.flatten()))
        all_enriched_combos.extend(enriched_combos)
        all_cols.extend(columns_to_use)
    all_cols = ["Input_" + c for c in all_cols]
    # Filter samples who have these combinations identified by rarecomb
    cond = " | ".join(["(" + " & ".join([f"({c} == 1)" for c in ec_i]) + ")" for ec_i in all_enriched_combos])
    wes_df = pd.read_csv(wes_file, sep="\t", low_memory=False, index_col=0, usecols=all_cols+["Sample_Name"])
    wes_df.columns = [c[len('Input_'):] for c in wes_df.columns]
    with open(selected_samples_file, "r") as f:
        all_samples = set([l.strip() for l in f.readlines()])
    wes_df = wes_df.loc[wes_df.eval(cond)]
    combo_samples = set([str(i) for i in list(wes_df.index)])
    combo_samples_in_selected_samples =  combo_samples.intersection(all_samples)
    return combo_samples_in_selected_samples, all_samples

class Node:
    """
    Each ICD10 diagnosis is stored as a Node object
    """
    def __init__(self, node_id, code, meaning, parent=None, child=None):
        self.node = node_id
        self.parent = parent
        self.child = child
        self.code, self.meaning = code, meaning

    def add_child(self, child_node):
        if self.child:
            self.child.append(child_node)
        else:
            self.child = [child_node]
        return

    def add_parent(self, parent_node):
        if not self.parent:
            self.parent = parent_node
        else:
            assert self.parent == parent_node
        return

    def get_parent(self):
        return self.parent

    def get_child(self):
        return self.child

    def get_info(self):
        return self.code, self.meaning


class Tree:
    def __init__(self, root_node):
        self.root = root_node
        self.node_dict = {self.root.node : self.root}
        self.icd_df = pd.DataFrame()

    def update_node_dict(self, node_id, node):
        if node_id not in self.node_dict:
            self.node_dict[node_id] = node
        return

    def create_node_from_df_helper(self, node_id):
        assert self.icd_df.empty == False
        c, m, ni, pi =  self.icd_df.loc[self.icd_df.node_id==node_id].values[0]
        n = Node(ni, c, m)
        return n, pi

    def create_node_from_df(self, node_id):
        if node_id in self.node_dict:
            return self.node_dict[node_id]

        # creating a node and providing parent information
        mn, mnpi = self.create_node_from_df_helper(node_id)
        # if parent is not present in the tree
        if mnpi not in self.node_dict:
            # create the parent node and get its parent
            mnp = self.create_node_from_df(mnpi)
            # add that parent info to the created node
            mn.add_parent(mnp)
        else:
            mnp = self.node_dict[mnpi]
            # add that parent info to the created node
            mn.add_parent(mnp)

        # update the node dict with the created node
        self.update_node_dict(node_id, mn)
        # add the created node as a child of the parent node
        mnp.add_child(mn)
        return mn

    def print_node(self, curr_node, node_level, tree_file):
        curr_node_info = curr_node.get_info()
        tree_file.write(f"{'-' * node_level}{curr_node.node}\t{curr_node_info[1]}\n")
        return

    def print_tree(self, curr_node, tree_file, node_level=0, max_node_level=2):
        if node_level>max_node_level:
            return
        
        if curr_node:
            self.print_node(curr_node, node_level, tree_file)

            if curr_node.child:
                for c in curr_node.child:
                    self.print_tree(c, tree_file, node_level+1, max_node_level)
        return

def get_icd_assigned_samples(code, sample_dir):
    samples = []
    file = os.path.join(sample_dir, f"{code[:1]}", f"{code}.txt")
    if os.path.exists(file):
        with open(file, "r") as f:
            samples = [l.strip() for l in f.readlines()]
    return set(samples)


def get_samples_with_icd_helper(node, sample_dir):
    all_samples = get_icd_assigned_samples(node.code, sample_dir)
    if node.child:
        for child in node.child:
            child_samples = get_samples_with_icd_helper(child, sample_dir)
            all_samples.update(child_samples)
    assert type(all_samples) == set
    return all_samples


def get_samples_with_icd(icd_tree, icd_code, icd2sample_dir):
    # Create code to node id reference
    code2nodeid_dict = {n.code:nid for nid,n in icd_tree.node_dict.items()}
    interested_node = icd_tree.node_dict[code2nodeid_dict[icd_code]]
    samples_with_code = get_samples_with_icd_helper(interested_node, icd2sample_dir)
    return samples_with_code


def create_icd_enrichment_table(combo_dfs, ncombos_list, wes_file, selected_samples_file, icd_df, icd2sample_dir, fishers_bash_path, save_file):
    # get samples with combinations
    combo_samples, all_samples = get_samples_with_combo(combo_dfs, ncombos_list, wes_file, selected_samples_file)
    noncombo_samples = all_samples.difference(combo_samples)
    # create the icd tree
    root_pheno = Node(0, "0", "Root Phenotype")
    icd_tree = Tree(root_pheno)
    # fill the tree with icd hierarchical info
    icd_tree.icd_df = icd_df
    for ni in icd_df.node_id:
        icd_tree.create_node_from_df(ni)
    
    stats=[]
    for icd_code in tqdm(icd_df.coding):
        # print(icd_code)
        icd_samples = get_samples_with_icd(icd_tree, icd_code, icd2sample_dir)
        icd_samples_selected = icd_samples.intersection(all_samples)
        # only do this if there are greater than 2000 individuals for a diagnosis
        if len(icd_samples_selected)>2000:
            combo_with_icd_samples = combo_samples.intersection(icd_samples_selected)
            combo_without_icd_samples = combo_samples.difference(icd_samples_selected)
            noncombo_with_icd_samples = noncombo_samples.intersection(icd_samples_selected)
            noncombo_without_icd_samples = noncombo_samples.difference(icd_samples_selected)

            contingency_table = [[len(combo_with_icd_samples), len(noncombo_with_icd_samples)],
                                [len(combo_without_icd_samples), len(noncombo_without_icd_samples)]]
            
            result = run_fishers_exact_in_r(contingency_table, fishers_bash_path)
            oddsratio = result[0]
            pvalue = result[1]
            ci_low = result[2]
            ci_high = result[3]
            stats.append([icd_code, oddsratio, pvalue, ci_low, ci_high, len(combo_with_icd_samples), len(combo_without_icd_samples), len(noncombo_with_icd_samples), len(noncombo_without_icd_samples)])
    # multiple testing for stats df
    stats = pd.DataFrame(stats, columns=['icd_phenotype', 'oddsratio', 'pvalue', 'conf_int_lower', 'conf_int_upper', 'Num_samples_with_combo_and_phenotype', 'Num_samples_with_combo_and_without_phenotype', 'Num_samples_without_combo_and_with_phenotype', 'Num_samples_without_combo_and_without_phenotype'])
    stats['FDR'] = multipletests(stats.pvalue, method='fdr_bh')[1]
    # save file
    stats.sort_values("FDR").to_csv(save_file, index=False)
    return stats

if __name__ == "__main__":
    icd_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_icd10/data/coding19.tsv"
    combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
    combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
    fishers_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/src/scripts/fishers_exact.sh"
    wes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/tables/wes.tsv"
    icd10codes2samples_dir = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/icd2sample"
    selected_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british/data/cases_controls/cases.txt"
    save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/icd_enrichment/icd_enrichment_all_risk_combos.csv"
    

    icd_df = pd.read_csv(icd_file, sep="\t", usecols=["coding", "meaning", "node_id", "parent_id"])
    combinations2_df = pd.read_csv(combinations2_file, low_memory=False)
    combinations3_df = pd.read_csv(combinations3_file, low_memory=False)


    stats = create_icd_enrichment_table([combinations2_df, combinations3_df], [2, 3], wes_file, selected_samples_file, icd_df, icd10codes2samples_dir, fishers_bash_path, save_file)
