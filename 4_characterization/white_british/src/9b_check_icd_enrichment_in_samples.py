import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import subprocess
import os
from tqdm import tqdm


def get_samples(samples_file):
    with open(samples_file, "r") as f:
        samples = set([l.strip() for l in f.readlines()])
    return samples

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


def create_icd_enrichment_table(combo_samples_file, all_samples_file, icd_df, icd2sample_dir, save_file):
    # get samples with combinations
    combo_samples = get_samples(combo_samples_file)
    all_samples = get_samples(all_samples_file)
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
            
            result = fisher_exact(contingency_table)
            oddsratio = result[0]
            pvalue = result[1]
            # ci_low = result[2]
            # ci_high = result[3]
            stats.append([icd_code, oddsratio, pvalue, len(combo_with_icd_samples), len(combo_without_icd_samples), len(noncombo_with_icd_samples), len(noncombo_without_icd_samples)])
    # multiple testing for stats df
    stats = pd.DataFrame(stats, columns=['icd_phenotype', 'oddsratio', 'pvalue', 'Num_samples_with_combo_and_phenotype', 'Num_samples_with_combo_and_without_phenotype', 'Num_samples_without_combo_and_with_phenotype', 'Num_samples_without_combo_and_without_phenotype'])
    stats['FDR'] = multipletests(stats.pvalue, method='fdr_bh')[1]
    # save file
    stats.sort_values("FDR").to_csv(save_file, index=False)
    return stats

if __name__ == "__main__":
    icd_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_icd10/data/coding19.tsv"
    combo_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/icd/combosamples_bmidecile1to3.txt"
    all_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/icd/allsamples_bmidecile1to3.txt"
    icd10codes2samples_dir = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/icd2sample"
    save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/icd/bmidecile1to3.csv"
    

    icd_df = pd.read_csv(icd_file, sep="\t", usecols=["coding", "meaning", "node_id", "parent_id"])

    stats = create_icd_enrichment_table(combo_samples_file, all_samples_file, icd_df, icd10codes2samples_dir, save_file)
