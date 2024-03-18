import pandas as pd
import os
import argparse
import tqdm

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


def get_samples_in_top_parent(node, tree, root):
    if node.parent==root:
        return node.get_samples()
    samples = get_samples_in_top_parent(node.parent, tree, root)
    return samples

def get_samples_in_parent(node, tree, root):
    if node.parent==root:
        return node.get_samples()
    return node.parent.get_samples()

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
    for icd_code, samples in tqdm.tqdm(zip(icd_samples_df.index, icd_samples_df.sample_names)):
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

def create_case_controls_file(icd_codes_of_interest, cohort_samples, pheno_tree, root_pheno, c2nodeid_dict, save_dir, case_cont_mode):
    for icdc in icd_codes_of_interest:
        all_samples = root_pheno.get_samples()
        icdc_node = pheno_tree.node_dict[c2nodeid_dict[icdc]]
        comorbid_samples = icdc_node.get_samples()
        comorbid_samples = cohort_samples.intersection(comorbid_samples)
        comorbid_samples_parent = get_samples_in_parent(icdc_node, pheno_tree, root_pheno)
        non_comorbid_samples = cohort_samples.difference(comorbid_samples)
        non_comorbid_parent_samples = cohort_samples.difference(comorbid_samples_parent)
        
        if case_cont_mode=="risk":
            # people in the cohort with the comorbidity
            case_samples = comorbid_samples
            # people in the cohort without the comorbidity or related parent phenotypes
            control_samples = non_comorbid_samples
        else:
            # this is not ready at all
            # people in the cohort without the comorbidity
            case_samples = non_comorbid_parent_samples
            # people in the cohort with the comorbidity
            control_samples = comorbid_samples
        
        assert len(case_samples.intersection(control_samples)) == 0

        print(icdc, len(case_samples), len(control_samples))
        icd_meaning = icdc_node.meaning.replace(" ", "")
        data_dict = {
            "Sample_Name": list(case_samples) + list(control_samples),
            f"Output_{icd_meaning}": [1 for i in range(len(case_samples))] + [0 for i in range(len(control_samples))]
        }
        df = pd.DataFrame(data_dict)
        save_file = os.path.join(save_dir, icd_meaning, "case_controls.csv")
        os.makedirs(os.path.dirname(save_file), exist_ok=True)
        df.to_csv(save_file, index=False)
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Rarecomb pipeline.')
    parser.add_argument("--icd_raw_dir", type=str, help="Filepath of the icd codes with sample info dir")
    parser.add_argument("--icd_codes_file", type=str, help="Filepath of the icd codes with parent info")
    parser.add_argument("--hes_info_file", type=str, help="Filepath of the in patient info file")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort samples")
    parser.add_argument("--save_dir", type=str, help="Filepath where case control files will be stored")
    parser.add_argument("--case_cont_mode", type=str, help="whether mining for risk or protective combos: one of risk/protective", default="risk")

    cli_args = parser.parse_args()

    icd_samples_df = create_icd_samples_file(cli_args.icd_raw_dir)
    cohort_samples_df = pd.read_csv(cli_args.cohort_file)
    icd_codes_df = pd.read_csv(cli_args.icd_codes_file, usecols=["coding", "meaning", "node_id", "parent_id"], sep="\t")
    hes_info_df = pd.read_csv(cli_args.hes_info_file, dtype={"sample_names": str, "hes_info": float})

    icd_codes_of_interest = ['Block M15-M19', 'Block K80-K87', 'E780',  'Block I20-I25', 'I10', 'E11', 'E039'] # Arthrosis, Gall bladder biliary tract pancreas, Pure hypercholesterolaemia, Ischaemic heart disease, Essential (primary) hypertension, Non-insulin-dependent diabetes mellitus, Hypothyroidism
    cohort_samples = set(cohort_samples_df.loc[:, "sample_names"].astype(str).values)
    all_icd_samples = set(hes_info_df.loc[hes_info_df.hes_info>0, "sample_names"].values)
    cohort_samples = cohort_samples.intersection(all_icd_samples)

    pheno_tree, root_pheno, c2nodeid_dict = create_tree(icd_codes_df, icd_samples_df)
    create_case_controls_file(icd_codes_of_interest, cohort_samples, pheno_tree, root_pheno, c2nodeid_dict, cli_args.save_dir, cli_args.case_cont_mode)
    