__author__ = "Feel Liao"
__copyright__ = "Copyright 2025, Feel Liao"
__email__ = "feel2027@outlook.com"
__license__ = "MIT"

from pathlib import Path
import pandas as pd
import tempfile
from snakemake.shell import shell


log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra = snakemake.params.get("extra", "")
acc_list_file = snakemake.input[0]
fastq_path = snakemake.output[0]
acc_list_file = Path(acc_list_file)
fastq_path = Path(fastq_path)

assert acc_list_file.exists(), f"File not found: {acc_list_file}"

# Download SRA data from ncbi
# TODO: Add parallel download option: When one SRA file is downloaded,
# it is converted to fastq format immediately.
with tempfile.TemporaryDirectory(prefix="sra_download_") as tmpdir:
    tmp_path = Path(tmpdir)

    try:
        shell("prefetch -O {tmp_path} "
              "--option-file {acc_list_file} "
              "{log} ")

        if not fastq_path.exists():
            fastq_path.mkdir(parents=True)

        downloaded_sras = tmp_path.rglob('*.sra')
        for sra in downloaded_sras:
            shell("parallel-fastq-dump -O {fastq_path} "
                  "-t {snakemake.threads} "
                  "--split-files -s {sra} "
                  "{log}")

    except Exception as e:
        print(f"Error processing SRR files : {str(e)}")
        raise


# generate sample sheet
fqFiles = fastq_path.glob("*.fastq")
samples = list(set([str(i).split("/")[-1].split("_")[0] for i in fqFiles]))

assert len(samples) > 0, "no files found in {} \
  please check download step".format(fastq_path)

# create sample sheet with pandas
sample_sheet = pd.DataFrame(samples, columns=["sample"])
sample_sheet["group"] = ""
sample_sheet["read1"] = sample_sheet["sample"].apply(
    lambda x: "{}_1.fastq".format(str(fastq_path)+"/"+x))
sample_sheet["read2"] = sample_sheet["sample"].apply(
    lambda x: "{}_2.fastq".format(str(fastq_path)+"/"+x))
sample_sheet["extra"] = ""

acc_parent = Path(acc_list_file).parent
outputFile = Path(acc_parent, "sample_sheet_sra.csv")

# save the sample sheet to the output path
sample_sheet.to_csv(outputFile, index=False)
