from pathlib import Path

import pandas as pd
from pathlib import Path
from snakemake.utils import validate

if not config["SRA"]["activate"]:
    validate(config, schema="../schemas/config.schema.yaml")
    
samples = (
    pd.read_csv(config["samples"], dtype="string")
    .set_index("sample", drop=False)
    .sort_index()
)

# validate(samples, schema="../schemas/samples.schema.yaml")

def get_accession():
    acc_list_file = config["SRA"]["acc_list"]
    assert Path(acc_list_file).exists(), f"{acc_list_file} does not exist"
    with open(acc_list_file) as f:
        acc_list = f.read().splitlines()
    return acc_list

def get_samples():
    return samples.index.tolist()

SAMPLES=get_samples()

def get_upstream_output():
    final_output = []
    final_output.append("out/quantification/samples_merged_counts.csv")
    final_output.append("out/quantification/samples_merged_tpm.csv")
    # final_output.extend(expand("out/featurecounts/{sample}.txt",sample=SAMPLES))

    if config["newGene"]["activate"]:
      final_output.append("out/newGene/merged_assembly.gtf")
      final_output.append("out/newGene/gffcompare.annotated.gtf")

    if config["reports"]:
      final_output.append("out/reports/report.html")

    return final_output

def get_post_input():
    post_output=[]
    post_output.append("out/quantification/samples_merged_counts.csv")
    post_output.append("out/quantification/samples_merged_tpm.csv")

    if config["newGene"]["activate"]:
        post_output.append("out/newGene/merged_assembly.gtf")
        post_output.append("out/newGene/gffcompare.annotated.gtf")
    return post_output



def get_fq(wildcards):
  file1 = samples["read1"].iloc[0]
  file2 = samples["read2"].iloc[0]
  file_path = Path(file1).parent
  end1 = file1.split("_")[-1]
  end2 = file2.split("_")[-1]
  fq1 = "{}/{}_{}".format(file_path,wildcards.sample,end1)
  fq2 = "{}/{}_{}".format(file_path,wildcards.sample,end2)
  return fq1, fq2