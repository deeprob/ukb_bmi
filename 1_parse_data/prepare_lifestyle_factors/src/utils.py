import os
import json
import pandas as pd
import numpy as np
import subprocess
from sklearn.decomposition import PCA


#######################
# reading input files #
#######################

def read_phenos_of_interest_data(file, min_samples=2000):
    """
    This function reads interesting phenotype file of type xlsx prepared manually, 
    selects the shortlisted phenotypes with at least 2000 exome samples and
    returns a filtered dataframe 
    """
    df = pd.read_excel(file)
    df = df.loc[df.shortlist=="X"]
    df = df.loc[df.Num_exome_samples_with_phenotype>=min_samples]
    return df

def read_exome_data(sample_to_exome_file):
    """
    This function reads the sample mapped to their exome vcf files,
    selects only those sample ids, which have both vcf and vcf index,
    return a filtered dataframe
    """
    df = pd.read_csv(sample_to_exome_file, index_col=0)
    df = df.dropna(how="all").reset_index()
    return df

def get_pheno_table_filepath(root_dir, pheno_type, pheno_cat, pheno_id):
    """
    This function accepts
    1) root_dir: where the sample to phenotype value table is stored under type -> category -> id hierarchy
    2) pheno_type: the type which the field belongs to
    3) pheno_cat: the category which the field belongs to 
    4) pheno_id: the UKBiobank id of the field
    It returns
    the concatenated filepath of the sample to field value table
    """
    pheno_table_path = os.path.join(
        root_dir, pheno_type, pheno_cat, "tables", f"{pheno_id}.csv"
        )
    assert os.path.exists(pheno_table_path)
    return pheno_table_path


def read_pheno_table(pheno_table_path):
    df =  pd.read_csv(pheno_table_path, index_col=0)
    return df


##############################
# filtering phenotype values #
##############################

def filter_pheno_table_no_negs(pheno_df_all, pheno_type):
    """
    Filter all samples which have negative field values
    Note: negative field values are associated with irrelavant information
    except in case of categorical multiple type of fields where -7 
    which denotes None of the above has a meaning 
    """
    if pheno_type in {"categorical_multiple"}:
        df = pheno_df_all.loc[(pheno_df_all==-7).any(axis=1)|(pheno_df_all>=0).any(axis=1)]
    else:
        df = pheno_df_all.loc[(pheno_df_all>=0).any(axis=1)]
    return df


########################
# merging field values #
########################

def merge_values_categorical(pheno_df_no_negative_vals, pheno_type):
    """
    For categorical_single:
    Merges all the field values for UKBiobank into a single field value 
    by taking the consensus of all field values. If there is not consensus, 
    it is marked as NaN.

    For categorical_multiple:
    Since these types of fields can take multiple values, there is not need to merge
    into a consensus value. Therefore it is returned as is for downstream processing
    """
    if pheno_type == "categorical_single":
        pheno_df_no_negative_vals_nona = pheno_df_no_negative_vals.fillna(method='bfill', axis=1).iloc[:,0]
        pheno_df_no_negative_vals_consistent = pheno_df_no_negative_vals.nunique(axis=1)==1
        pheno_df_no_negative_vals_consensus = pd.Series([i if j==True else np.nan for i,j in zip(pheno_df_no_negative_vals_nona, pheno_df_no_negative_vals_consistent)])
        assert len(pheno_df_no_negative_vals) ==  len(pheno_df_no_negative_vals_consensus)
        pheno_df_no_negative_vals["merged"] = pheno_df_no_negative_vals_consensus.values
        # drop rows that do not reach a consensus
        pheno_df_no_negative_vals = pheno_df_no_negative_vals.dropna(subset=["merged"])
    return pheno_df_no_negative_vals

def merge_values_numerical(pheno_df_no_negative_vals):
    """
    Merges all the field values for UKBiobank into a single field value 
    by taking the mean of all field values.
    """
    pheno_df_no_negative_vals["merged"] = pheno_df_no_negative_vals.mean(axis=1, numeric_only=True)
    return pheno_df_no_negative_vals


#####################
# binarizing values #
#####################

def binarize_numericals(df, strategy="median", quantile_low=0.25, quantile_high=0.75):

    ser = df["merged"]

    if strategy == "median":
        thresh = ser.median()
        binarized_low = (ser<=thresh).astype(int)
        binarized_high = (ser>thresh).astype(int)
        df[f"binarized_{thresh}_low"] = binarized_low
        df[f"binarized_{thresh}_high"] = binarized_high

    elif strategy == "quantile":
        qlow = ser.quantile(quantile_low)
        qhigh = ser.quantile(quantile_high)

        binarized_low = (ser <= qlow).astype(int)
        binarized_high = (ser >= qhigh).astype(int)
        df[f"binarized_{qlow}_low"] = binarized_low
        df[f"binarized_{qhigh}_high"] = binarized_high

    return df


def get_ordinal_categorical_bins(df, field_encodings):
    
    ser = df["merged"]

    field_encodings_relevant = sorted([int(fe) for fe in field_encodings.keys() if int(fe)>=0])
    # case 1, if there are two types of relevant categories for this field,
    # just return low and high
    if len(field_encodings_relevant) == 2:
        # low category
        fe_low = field_encodings_relevant[0]
        fe_low_val = field_encodings[str(fe_low)]
        fe_low_val = "-".join(fe_low_val.split())
        # high category
        fe_high = field_encodings_relevant[-1]
        fe_high_val = field_encodings[str(fe_high)]
        fe_high_val = "-".join(fe_high_val.split())
        # binarize
        binarized_low = (ser == fe_low).astype(int)
        binarized_high = (ser == fe_high).astype(int)

    # case 2, if there are more than two types of relevant categories for this field
    if len(field_encodings_relevant) > 2:
        # check how many low bins will give us at least 10% samples 
        # excluding the highest bin
        for lower_bin_range in range(1, len(field_encodings_relevant)):
            fe_lows = field_encodings_relevant[:lower_bin_range]
            fe_low_val = "|".join(["-".join(field_encodings[str(fe_low)].replace(",", "").split()) for fe_low in fe_lows])
            binarized_low = ser.isin(fe_lows).astype(int)
            if (sum(binarized_low)/len(binarized_low)) > 0.1:
                break
        
        # check how many high bins will give us at least 10% samples 
        # excluding all the bins selected by the lower bin range
        for higher_bin_range in range(len(field_encodings_relevant) - 1, lower_bin_range - 1, -1):
            fe_highs = field_encodings_relevant[higher_bin_range:]
            fe_high_val = "|".join(["-".join(field_encodings[str(fe_high)].replace(",", "").split()) for fe_high in fe_highs])            
            binarized_high = ser.isin(fe_highs).astype(int)
            if (sum(binarized_high)/len(binarized_high)) > 0.1:
                break
    
    # check if the number of samples in the lowest and highest bins 
    # are less than 10% of the original number of samples
    if (sum(binarized_low)/len(binarized_low)) < 0.1:
        print(f"Warning:: lowest category has less than 10% samples for field id {df.columns[0]}")
        print(f"Warning:: lowest category has {sum(binarized_low)} samples")

    if (sum(binarized_high)/len(binarized_high)) < 0.1:
        print(f"Warning:: highest category has less than 10% samples for field id {df.columns[0]}")
        print(f"Warning:: highest category has {sum(binarized_high)} samples w/ or w/o exomes")

    df[f"binarized_{fe_low_val}_low"] = binarized_low
    df[f"binarized_{fe_high_val}_high"] = binarized_high

    return df


def get_ohe_type_encoding(df, field_encodings):
    # only get the relevant field encoding values
    field_encodings_relevant = sorted([int(fe) for fe in field_encodings.keys() if int(fe)>=0])
    
    # add "None of the above", encoding value "-7", here if it is present
    if "-7" in field_encodings.keys():
        field_encodings_relevant.append(-7)
    df_copy = df.copy()
    for fe in field_encodings_relevant:
        fe_ser = df.isin([fe]).any(axis=1).astype(int)
        fe_val = field_encodings[str(fe)]
        fe_val = "-".join(fe_val.replace(",", "").split())
        df_copy[f"binarized_{fe_val}"] = fe_ser    
    return df_copy


def binarize_categoricals(df, field_type, field_encodings, ordinal_status, ohe_encodings, ordinal_encodings):

    if field_type == "categorical_single":

        if not pd.isnull(ordinal_status):
            # this categorical single phenotype is not ordinal
            encoding_type = ordinal_status
            if encoding_type == "O":
                # one hot encoded type
                df = get_ohe_type_encoding(df, field_encodings)

            elif encoding_type == "B":
                print("Warning:: Still working on it!!")
                print("Warning:: Encoding might be erroneous if ohe type field encodings have 1s!!")
                # read new field encodings for high low 
                df = get_ordinal_categorical_bins(df, ordinal_encodings)
                # read new field encodings for ohe type
                df = get_ohe_type_encoding(df, ohe_encodings)
        
        else:
            df = get_ordinal_categorical_bins(df, field_encodings)

    elif field_type == "categorical_multiple":
        if type(field_encodings) == dict:
            df = get_ohe_type_encoding(df, field_encodings)
        else:
            print(f"Warning :: Field id {df.columns[0]} has incorrect field encoding type: {field_encodings}.")
            print(f"Warning :: It will not be binarized.")    

    return df


########################
# filtering for exomes #
########################

def filter_pheno_with_exomes(binarized_df, exome_df):
    binarized_df = binarized_df.loc[binarized_df.index.isin(exome_df.eid)]
    return binarized_df


##########################
# saving binarized table #
##########################

def save_pheno_table(
    binarized_pheno_df, storage_root_dir, 
    pheno_type, pheno_cat, pheno_id, 
    method):

    table_basename = f"{pheno_id}_{method}" if method else str(pheno_id)
    binarized_pheno_path = os.path.join(
        storage_root_dir, pheno_type, pheno_cat, "tables", f"{table_basename}.csv"
        )
    os.makedirs(os.path.dirname(binarized_pheno_path), exist_ok=True)
    binarized_pheno_df.to_csv(binarized_pheno_path)
    return


##########################
# field encodings parser #
##########################

def get_pheno_encoding_filepath(root_dir, pheno_type, pheno_cat):
    """
    This function accepts
    1) root_dir: where the sample to phenotype value table is stored under type -> category -> id hierarchy
    2) pheno_type: the type which the field belongs to
    3) pheno_cat: the category which the field belongs to 
    It returns
    the filepath of the json file that contains fields encodings of all fields 
    that fall under the type and category specified
    """
    pheno_json_path = os.path.join(
        root_dir, pheno_type, pheno_cat, f"fields_data_coding.json"
        )
    assert os.path.exists(pheno_json_path)
    return pheno_json_path


def get_modified_pheno_encoding_filepath(root_dir, enc_type):
    """
    This function accepts
    1) root_dir: where the sample to phenotype value table is stored under type -> category -> id hierarchy
    2) pheno_type: the type which the field belongs to
    3) pheno_cat: the category which the field belongs to 
    It returns
    the filepath of the json file that contains fields encodings of all fields 
    that fall under the type and category specified
    """
    pheno_json_path = os.path.join(
        root_dir, "modified_field_encodings", enc_type, f"field_encodings.json"
        )
    assert os.path.exists(pheno_json_path)
    return pheno_json_path


def read_pheno_encodings(pheno_json_path, pheno_field_id):
    pheno_field_id = str(pheno_field_id)
    with open(pheno_json_path, "r") as f:
        field_encoding_dict = json.load(f)
    return field_encoding_dict[pheno_field_id]


#########################
# meta table formatting #
#########################

def get_exome_index(exome_file):
    return pd.read_csv(exome_file, index_col=0).dropna(how="all").index


def get_binarized_table_path(root_dir, pheno_type, pheno_cat, pheno_id, strategy):
    pheno_basename = f"{pheno_id}_{strategy}.csv" if strategy else f"{pheno_id}.csv"
    pheno_binarized_table_path = os.path.join(
        root_dir, 
        pheno_type, 
        pheno_cat, 
        "tables", 
        pheno_basename
        )
    return pheno_binarized_table_path


def read_binarized_table(table_path):
    df = pd.read_csv(table_path, index_col=0)
    df_columns = [c for c in list(df.columns) if c.startswith("binarized_")]
    # checking to see that there are at least 2 binarized columns
    if len(df_columns)<2:
        print(f"Warning:: {table_path} was not binarized")
    df = df.loc[:, df_columns]
    return df


def reindex_binarized_table1(pheno_df, pheno_id, exome_index):
    pheno_df_old_columns = pheno_df.columns
    pheno_df_new_columns = [f"Input_{pheno_id}_{c.split('_')[-1]}" for c in list(pheno_df_old_columns)]
    # keeping track of old column names 
    col_data = {
        "old": pheno_df_old_columns,
        "new": pheno_df_new_columns
    }
    col_df = pd.DataFrame(col_data)
    pheno_df.columns = pheno_df_new_columns
    return pheno_df.reindex(exome_index), col_df


def get_field_encodings(root_dir, pheno_type, pheno_cat, pheno_id):
    pheno_json_path = get_pheno_encoding_filepath(root_dir, pheno_type, pheno_cat)
    pheno_encodings = read_pheno_encodings(pheno_json_path, pheno_id)
    if type(pheno_encodings) == dict:
        pheno_encodings = {"-".join(v.replace(",", "").split()):k for k,v in pheno_encodings.items()}
    return pheno_encodings


def reindex_binarized_table2(pheno_df, pheno_id, pheno_encodings, exome_index):
    pheno_df_old_columns = pheno_df.columns
    pheno_df_new_columns = [f"Input_{pheno_id}_{pheno_encodings[c.split('_')[-1]]}" for c in list(pheno_df_old_columns)]
    # keeping track of old column names 
    col_data = {
        "old": pheno_df_old_columns,
        "new": pheno_df_new_columns
    }
    col_df = pd.DataFrame(col_data)
    pheno_df.columns = pheno_df_new_columns    
    return pheno_df.reindex(exome_index), col_df


def get_modified_field_encodings(root_dir, enc_type, pheno_id):
    pheno_json_path = get_modified_pheno_encoding_filepath(root_dir, enc_type)
    pheno_encodings = read_pheno_encodings(pheno_json_path, pheno_id)
    if type(pheno_encodings) == dict:
        pheno_encodings = {"-".join(v.replace(",", "").split()):k for k,v in pheno_encodings.items()}
    return pheno_encodings


def reindex_binarized_table3(pheno_df, pheno_id, ohe_encodings, exome_index):
    pheno_df_old_ohe_columns = [c for c in pheno_df.columns if all(x!=c.split("_")[-1] for x in ["low", "high"])]
    pheno_df_new_ohe_columns = [f"Input_{pheno_id}_{ohe_encodings[c.split('_')[-1]]}" for c in list(pheno_df_old_ohe_columns)]
    pheno_df_old_ordinal_columns = [c for c in pheno_df.columns if any(x==c.split("_")[-1] for x in ["low", "high"])]
    pheno_df_new_ordinal_columns = [f"Input_{pheno_id}_{c.split('_')[-1]}" for c in list(pheno_df_old_ordinal_columns)]    
    # keeping track of old column names 
    pheno_df_old_columns = pheno_df_old_ohe_columns + pheno_df_old_ordinal_columns
    pheno_df_new_columns = pheno_df_new_ohe_columns + pheno_df_new_ordinal_columns
    col_data = {
        "old": pheno_df_old_columns,
        "new": pheno_df_new_columns
    }
    col_df = pd.DataFrame(col_data)
    pheno_df = pheno_df.rename(columns = dict(zip(pheno_df_old_columns, pheno_df_new_columns)))    
    return pheno_df.reindex(exome_index), col_df

##############################
# read the lifestyle factors #
##############################

def get_continuous_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="continuous")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}_quantile.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict


def get_integer_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="integer")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}_quantile.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict

def get_catsingle_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="categorical_single")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict

def get_catmultiple_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="categorical_multiple")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}.csv"), usecols=lambda x: x.startswith("eid") or x.startswith("binarized"), index_col=["eid"]).rename(columns=lambda col: col.replace("binarized", f"{i}")) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict


#############
# normalize #
#############

def standardize_df(df):
    # df = (df - df.mean())/df.std()
    df = (df-df.min())/(df.max()-df.min())
    return df

def replace_negs(df):
    return df.mask(df<0)

def feature_pca(df):
    pca = PCA(n_components=2)
    Xt = pca.fit_transform(df.values)
    plot_df = pd.DataFrame()
    plot_df.index = df.index
    plot_df["component_1"] = Xt[:,0]
    plot_df["component_2"] = Xt[:,1]
    plot_df = plot_df.merge(df, left_index=True, right_index=True)
    return plot_df

def normalize_single_field_function(standardized_df):
    pca_df = pd.DataFrame()
    pca_df["component_1"] = standardized_df.iloc[:, 0]
    pca_df = pca_df.merge(standardized_df, left_index=True, right_index=True)
    return pca_df

def normalize_integers(int_df, sel_fields):
    # select relevant fields
    selected_df = int_df.loc[:, sel_fields]
    # get rid of negative values
    selected_noneg_df = replace_negs(selected_df).dropna(axis=0)
    # standardize df
    selected_noneg_standardized_df = standardize_df(selected_noneg_df)
    # get pca components
    if len(sel_fields) == 1:
        pca_df = normalize_single_field_function(selected_noneg_standardized_df)
    else:
        pca_df = feature_pca(selected_noneg_standardized_df)    
    return pca_df

def normalize_catsingle_helper(df, field_encodings, onesided_fields):
    # TODO: high and low encodings 2 columns per field if flag else no
    df_encoded = pd.DataFrame()
    for coln, cold in df.iteritems():
        fe_coln = field_encodings[coln]
        field_encodings_relevant = sorted([int(fe) for fe in fe_coln.keys() if int(fe)>=0])
        # case 1, if there are two types of relevant categories for this field,
        # just return low and high
        if len(field_encodings_relevant) == 2:
            df_encoded[coln] = cold
        # case 2, if there are more than two types of relevant categories for this field
        elif len(field_encodings_relevant) > 2:
            # check if the field is in onesided fields
            if coln in onesided_fields:
                # ignore the first encoding assuming it's answer is no and take all the other encodings into account
                # example bipolar disorder: highs include all types of bipolar disorder
                fe_highs = field_encodings_relevant[1:]
                fe_high_val = "|".join(["-".join(fe_coln[str(fe_high)].replace(",", "").split()) for fe_high in fe_highs])          
                binarized_high = cold.isin(fe_highs).astype(int) 
                df_encoded[f"{coln}_{fe_high_val}_high"] = binarized_high              
            else:
                # check how many low bins will give us at least 25% samples 
                # excluding the highest bin
                for lower_bin_range in range(1, len(field_encodings_relevant)):
                    fe_lows = field_encodings_relevant[:lower_bin_range]
                    fe_low_val = "|".join(["-".join(fe_coln[str(fe_low)].replace(",", "").split()) for fe_low in fe_lows])
                    binarized_low = cold.isin(fe_lows).astype(int)
                    if (sum(binarized_low)/len(binarized_low)) >= 0.25:
                        break        
                # check how many high bins will give us at least 25% samples 
                # exclude the lowest bin
                for higher_bin_range in range(len(field_encodings_relevant) - 1, lower_bin_range - 1, -1):
                    fe_highs = field_encodings_relevant[higher_bin_range:]
                    fe_high_val = "|".join(["-".join(fe_coln[str(fe_high)].replace(",", "").split()) for fe_high in fe_highs])          
                    binarized_high = cold.isin(fe_highs).astype(int)
                    if (sum(binarized_high)/len(binarized_high)) >= 0.25:
                        break
                # check to see if binarized low is the exact opposite of binarized high, then only keep binarized high
                if (binarized_low!=1).astype(int).equals(binarized_high):
                    df_encoded[f"{coln}_{fe_high_val}_high"] = binarized_high
                else:
                    df_encoded[f"{coln}_{fe_low_val}_low"] = binarized_low
                    df_encoded[f"{coln}_{fe_high_val}_high"] = binarized_high
    return df_encoded

def normalize_catsingle(catsingle_df, sel_fields, field_encodings, onesided_fields):
    # select relevant fields
    selected_df = catsingle_df.loc[:, sel_fields]
    # get rid of negative values
    selected_noneg_df = replace_negs(selected_df).dropna(axis=0)
    # get dummies
    selected_noneg_dummies_df = normalize_catsingle_helper(selected_noneg_df, field_encodings, onesided_fields) #pd.get_dummies(selected_noneg_df, columns=sel_fields, drop_first=True)
    # get pca components
    if selected_noneg_dummies_df.shape[1] == 1:
        pca_df = normalize_single_field_function(selected_noneg_dummies_df)
    else:
        pca_df = feature_pca(selected_noneg_dummies_df)    
    return pca_df

def normalize(field_name, sel_fields, field_type_dict, onesided_fields):
    if "integer" in field_name:
        return normalize_integers(field_type_dict["integer"][0], sel_fields)
    
    elif "catsingle" in field_name:
        return normalize_catsingle(field_type_dict["catsingle"][0], sel_fields, field_type_dict["catsingle"][1], onesided_fields)
    
    elif "catmultiple" in field_name:
        return normalize_catsingle(field_type_dict["catmultiple"], sel_fields, onesided_fields)

    return

##########
# encode #
##########

def create_labels(df, lowq=0.25, highq=0.75):
    binary_columns = list(df.columns[df.isin([0, 1, np.nan]).all()])
    non_binary_columns = [c for c in df.columns if c not in binary_columns]
    non_binary_df = df.loc[:, non_binary_columns]
    qdf = non_binary_df.quantile([lowq, highq])
    label_df = non_binary_df.copy()
    label_df = label_df.mask(non_binary_df<qdf.iloc[0, :], other="low") # stringent cutoff for the low field because often these are normal
    label_df = label_df.mask(non_binary_df>=qdf.iloc[1, :], other="high")
    label_df = label_df.mask((non_binary_df>=qdf.iloc[0, :])&(non_binary_df<qdf.iloc[1, :]), other="medium")
    return label_df, df.loc[:, binary_columns]
