import os
import pandas as pd


def create_gmt_file(go_terms, gmt_file, out_file):
    save_out = open(out_file, "w")
    with open(gmt_file, "r") as gmtf:
        for line in gmtf:
            go_term_label = line.split('\t')[1]
            term_id = line.split('\t')[0]
            if any([x in term_id for x in ["%GOBP%", "%GOCC%", "%GOMF%"]]):
                if go_term_label in go_terms:
                    save_out.write(line)
    save_out.close()
    return


if __name__ == "__main__":

    gmt_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/misc/enrichment_map_files/Human_GO_AllPathways_no_GO_iea_February_08_2023_symbol.gmt"
    input_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/enrichment/go/go_enrichment.csv"
    out_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/enrichment/go/go_enrichment.gmt"

    df = pd.read_csv(input_file)
    go_terms_to_include = df["Description"].to_list()
    create_gmt_file(go_terms_to_include, gmt_file, out_file)
