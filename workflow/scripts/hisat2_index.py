"""Snakemake wrapper for HISAT2 index"""

__author__ = "Feel Liao"
__copyright__ = "Copyright 2025, Feel Liao"
__email__ = "feel2027@outlook.com"
__license__ = "MIT"

from snakemake.shell import shell
from pathlib import Path
import tempfile

# Creating log
# snakemake.script
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Placeholder for optional parameters
extra = snakemake.params.get("extra", "")

# Allowing for multiple FASTA files
fasta = snakemake.input.get("fasta")
assert fasta is not None, "input-> a FASTA-file or a sequence is required"


compressed = False
if Path(fasta).suffix == ".gz":
    # If the input is a gzipped file, we need to decompress it
    compressed = True

input_seq = fasta

# get common prefix
prefix = Path(snakemake.output[0])

if not prefix.exists():
    prefix.mkdir(parents=True)

if compressed:
    with tempfile.NamedTemporaryFile(
            mode="w+b", delete=True, suffix=".fa") as tmp_fa:
        tmp_seq = tmp_fa.name
        shell("pigz -p {snakemake.threads} -d {input_seq} "
              "-c > {tmp_seq} && "
              "hisat2-build -p {snakemake.threads} "
              "{extra} {tmp_seq} {prefix}/genome {log}")
else:
    shell("hisat2-build -p {snakemake.threads} "
          "{extra} {input_seq} {prefix}/genome {log}")
