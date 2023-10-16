import pandas as pd


def parse_gaf_file(gaf_file, save_file):
    df = pd.read_csv(gaf_file, sep="\t", comment="!", header=None, usecols=[2, 4], names=["gene", "go_term"])
    df.groupby("gene").aggregate(lambda x: ";".join(x)).to_csv(save_file, index=True, sep="\t", header=False)
    return


if __name__=="__main__":
    gaf_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/go/goa_human.gaf"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/go/goa_human.annot.tab"
    parse_gaf_file(gaf_file, save_file)

    