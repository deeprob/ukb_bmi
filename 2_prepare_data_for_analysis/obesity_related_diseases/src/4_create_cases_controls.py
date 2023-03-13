#!/bin/python3
import pandas as pd


########################
# select obese samples #
########################

corrected_file_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/samples_with_residuals.csv"

df = pd.read_csv(corrected_file_samples_file)
df = df.sort_values('bmi_residuals')
df['bmi_decile'] = pd.qcut(df['bmi_residuals'], q=10)

lowest_group = df['bmi_decile'].unique()[[0,1,2]]
highest_group = df['bmi_decile'].unique()[[-2, -1]]

samples_in_lowest_group = set(df[df['bmi_decile'].isin(lowest_group)]['eid'].astype(str).to_list())
samples_in_highest_group = set(df[df['bmi_decile'].isin(highest_group)]['eid'].astype(str).to_list())


#########################################################
# get samples from the UKBiobank with related disorders #
#########################################################

import os
import pandas as pd

# input files
icd_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_icd10/data/coding19.tsv"
icd10codes2samples_dir = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/icd2sample"
annot_dir = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/annotated_by_sample"


df_pheno = pd.read_csv(icd_file, usecols=["coding", "meaning", "node_id", "parent_id"], sep="\t")


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
        self.code2samples = dict()

    def update_node_dict(self, node_id, node):
        if node_id not in self.node_dict:
            self.node_dict[node_id] = node
        return

    def create_node_from_df_helper(self, node_id):
        c, m, ni, pi =  df_pheno.loc[df_pheno.node_id==node_id].values[0]
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


# plant the tree
root_pheno = Node(0, "0", "Root Phenotype")
pheno_tree = Tree(root_pheno)

# fill the tree with leaves and branches - takes 6 secs
for ni in df_pheno.node_id:
    pheno_tree.create_node_from_df(ni)

# print the tree
tree_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/icd_tree.txt"
pt = open(tree_file, "w")
pt.close()
tf = open(tree_file, "a")
pheno_tree.print_tree(root_pheno, tf, max_node_level=4)
tf.close()

# Create code to node id reference
code2nodeid_dict = {n.code:nid for nid,n in pheno_tree.node_dict.items()}
# get all samples of a node

def get_icd_assigned_samples(code, sample_dir):
    samples = []
    file = os.path.join(sample_dir, f"{code[:1]}", f"{code}.txt")
    if os.path.exists(file):
        with open(file, "r") as f:
            samples = [l.strip() for l in f.readlines()]
    return set(samples)

def get_samples(node, sample_dir):
    all_samples = get_icd_assigned_samples(node.code, sample_dir)
    if node.child:
        for child in node.child:
            child_samples = get_samples(child, sample_dir)
            all_samples.update(child_samples)
    assert type(all_samples) == set
    return all_samples

icd_codes_of_interest = ['I10', 'E780', 'R074', 'I251', 'I259', 'E039', 'E11', 'Block M15-M19', 'K80', 'K81', 'K82', 'F32', 'F33', 'G30']
# save icd codes to look at table
icd_interest_save_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/tables/obesity_related_diseases_icd.csv"
df_pheno_interest = df_pheno.loc[df_pheno.coding.isin(icd_codes_of_interest)]

icd2nsamples_dict = {}

for icd_code in icd_codes_of_interest:
    print(icd_code)
    interested_node = pheno_tree.node_dict[code2nodeid_dict[icd_code]]
    samples_with_code = get_samples(interested_node, icd10codes2samples_dir)
    print(len(samples_with_code))
    # from the obese group select samples with the code
    cases_obese = samples_in_highest_group.intersection(samples_with_code)
    controls_obese = samples_in_highest_group.difference(samples_with_code)
    assert len(cases_obese)+len(controls_obese) == len(samples_in_highest_group)
    icd2nsamples_dict[icd_code] = len(cases_obese)
    icd_code = icd_code.replace(" ", "")
    cases_file = f"/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/cases_controls/cases_{icd_code}.txt"
    controls_file = f"/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/cases_controls/controls_{icd_code}.txt"
    with open(cases_file, 'w') as f:
        for sample in list(cases_obese):
            f.write(sample+'\n')
    with open(controls_file, 'w') as f:
        for sample in list(controls_obese):
            f.write(sample+'\n')

df_pheno_interest["ncases"] = df_pheno_interest.coding.map(icd2nsamples_dict)
df_pheno_interest.to_csv(icd_interest_save_file, index=False)
