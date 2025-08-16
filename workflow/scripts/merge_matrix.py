__author__ = "Feel Liao"
__copyright__ = "Copyright 2025, Feel Liao"
__email__ = "feel2027@outlook.com"
__license__ = "MIT"

from pathlib import Path
import pandas as pd

input_dir = snakemake.params.get("input_dir")
merge_matrix = snakemake.output.get("merged")
tpm_matrix = snakemake.output.get("tpm")
# fpkm_matrix = snakemake.output.get("fpkm")

expression_files = Path(input_dir).glob("*.txt")
merge_matrix = Path(merge_matrix)
tpm_matrix = Path(tpm_matrix)
expressionL = [str(i) for i in expression_files]
assert len(expressionL) > 0, "No files found in the {}".format(input_dir)

print(expressionL)
# Merge the expression files
read1 = pd.read_csv(expressionL[0], sep="\t", header=0, comment="#")
read1 = read1.iloc[:, [0, 5, 6]]


def rename_columns_temp(df, num):
    df_columns = df.columns.tolist()
    columns_rename = str(expressionL[num]).split(
        "/")[-1].split(".")[0]
    return df_columns[len(df_columns)-1], columns_rename


def rename_columns(df):
    df_columns = df.columns.tolist()
    columns_rename = df_columns[len(df_columns)-1].split(".")[0]
    return df_columns[len(df_columns)-1], columns_rename


def simple_df(df):
    df = df.iloc[:, [0, 6]]
    return df


read1_names = rename_columns_temp(read1, 0)
read1.rename(columns={read1_names[0]: read1_names[1]}, inplace=True)

for i in range(1, len(expressionL)):
    reads = pd.read_csv(
        expressionL[i], sep="\t", header=0, comment="#")
    reads = simple_df(reads)
    reads_names = rename_columns_temp(reads, i)
    reads.rename(columns={reads_names[0]: reads_names[1]}, inplace=True)
    read1 = pd.merge(read1, reads, on="Geneid")

read1.rename(columns={"Geneid": "gene_id"}, inplace=True)

counts_matrix = read1.copy()
counts_matrix.drop(columns=["Length"], inplace=True)
counts_matrix.to_csv(merge_matrix, index=False)


# TPM calculation
def tpm_calc(counts: pd.DataFrame):

    counts.set_index(counts.columns[0], inplace=True)

    cols = counts.columns[1:len(counts.columns)]
    rpk = counts[cols].div(counts.iloc[:, 0], axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1)
    return tpm.map(lambda x: round(x*1e6, 5))


tpm_out = tpm_calc(read1)
tpm_out.to_csv(tpm_matrix)
