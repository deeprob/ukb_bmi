import os
import pandas as pd
from functools import reduce


def get_gene_set_from_file(gene_file):
    with open(gene_file, "r") as f:
        genes = set([g.strip() for g in f.readlines()])
    return genes


def get_samples_with_gene_mutation(genotype_df, genes):
    all_gene_samples = set(",".join(genotype_df.loc[genotype_df.gene.isin(genes)].samples.values).split(","))
    return all_gene_samples


def get_case_cont_samples(case_cont_df):
    case_samples = set(case_cont_df.loc[case_cont_df.Output_BMI==1, "Sample_Name"].astype("str").values)
    control_samples = set(case_cont_df.loc[case_cont_df.Output_BMI==0, "Sample_Name"].astype("str").values)
    return case_samples, control_samples


def get_combo_samples_from_file(combo_file):
    combo_df = pd.read_csv(combo_file)
    all_combo_samples = set("|".join(combo_df.combo_samples.dropna().values).split("|"))
    return all_combo_samples


def get_combo_genes_from_file(combo_file):
    combo_df = pd.read_csv(combo_file)
    combo_df["uniq_items"] = combo_df.uniq_items.apply(lambda x: x.replace("Input_", ""))
    combo_genes = set("|".join(combo_df.uniq_items.values).split("|"))
    return combo_genes


def get_combo_info_from_files(combo_files):
    combo_genes = [get_combo_genes_from_file(cf) for cf in combo_files]
    combo_genes = reduce(lambda x,y: x.union(y), combo_genes)
    combo_samples = [get_combo_samples_from_file(cf) for cf in combo_files]
    combo_samples = reduce(lambda x,y: x.union(y), combo_samples)
    return combo_genes, combo_samples


def concat_multiple_combo_files(combo_files):
    combo_df = pd.concat([pd.read_csv(cf) for cf in combo_files])
    return combo_df

############
# ICD tree #
############
class Node:
    """
    Each ICD10 diagnosis is stored as a Node object
    """
    def __init__(self, node_id, code, meaning, parent=None, child=None):
        self.node = node_id
        self.parent = parent
        self.child = child
        self.code, self.meaning = code, meaning
        self.samples = set()

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
    
    def get_samples(self):
        return self.samples
    
    def get_samples_number(self):
        return len(self.samples)

class Tree:
    def __init__(self, root_node, coding_df):
        self.root = root_node
        self.node_dict = {self.root.node : self.root}
        self.coding_df = coding_df

    def update_node_dict(self, node_id, node):
        if node_id not in self.node_dict:
            self.node_dict[node_id] = node
        return

    def create_node_from_df_helper(self, node_id):
        c, m, ni, pi =  self.coding_df.loc[self.coding_df.node_id==node_id].values[0]
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
        if tree_file:
            tree_file.write(f"{'-' * node_level}{curr_node.node}\t{curr_node_info[1]}\n")
        else:
            print(f"{'-' * node_level}{curr_node.node}\t{curr_node_info[1]}\n")
        return

    def print_tree(self, curr_node, tree_file="", node_level=0, max_node_level=2):
        if node_level>max_node_level:
            return
        
        if curr_node:
            self.print_node(curr_node, node_level, tree_file)

            if curr_node.child:
                for c in curr_node.child:
                    self.print_tree(c, tree_file, node_level+1, max_node_level)
        return
    
    def add_sample_info(self, node_id, samples):
        curr_node = self.node_dict[node_id]
        curr_node.samples = samples.union(curr_node.samples)
        if curr_node.parent:
            self.add_sample_info(curr_node.parent.node, samples)
        return
    
def create_tree(icd_codes_df, icd_samples_df):
    # create tree
    # plant the tree
    root_pheno = Node(0, "0", "Root Phenotype")
    pheno_tree = Tree(root_pheno, icd_codes_df)
    # fill the tree with leaves and branches - takes 6 secs
    for ni in icd_codes_df.node_id:
        pheno_tree.create_node_from_df(ni)
    c2nodeid_dict = dict(zip(icd_codes_df.coding, icd_codes_df.node_id))
    # add sample info
    for icd_code, samples in zip(icd_samples_df.index, icd_samples_df.sample_names):
        pheno_tree.add_sample_info(c2nodeid_dict[icd_code], set(samples.split("|")))
    return pheno_tree, root_pheno, c2nodeid_dict

def create_icd_samples_file(icd_raw_dir):
    dfs = []
    for file in os.scandir(icd_raw_dir):
        filepath = os.path.join(icd_raw_dir, file)
        df = pd.read_csv(filepath)
        dfs.append(df)
    icd_samples_df = pd.concat(dfs)
    icd_samples_df["icd"] = icd_samples_df.icd.str.split("|")
    icd_samples_df = icd_samples_df.explode("icd").groupby("icd").agg(lambda x: "|".join(map(str,x)))
    return icd_samples_df
