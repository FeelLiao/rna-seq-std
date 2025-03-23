__author__ = "Sebastian Kurscheid"
__copyright__ = "Copyright 2019, Sebastian Kurscheid"
__email__ = "sebastian.kurscheid@anu.edu.au"
__license__ = "MIT"

from snakemake.shell import shell

# dict.get() returns "" if key not found
extra = snakemake.params.get("extra", "")
adapters = snakemake.params.get("adapters", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


# Assert input
n = len(snakemake.input)
assert (
    n == 1 or n == 2
), "input->sample must have 1 (single-end) or 2 (paired-end) elements."


# Input files
if n == 1:
    reads = "--in1 {}".format(snakemake.input)
else:
    reads = "--in1 {} --in2 {}".format(*snakemake.input)


# Output files
trimmed_paths = snakemake.output.get("trimmed", None)
if trimmed_paths:
    if n == 1:
        trimmed = "--out1 {}".format(snakemake.output.trimmed)
    else:
        trimmed = "--out1 {} --out2 {}".format(*snakemake.output.trimmed)
else:
    trimmed = ""

# Output failed reads
failed = snakemake.output.get("failed", None)
if failed:
    trimmed += f" --failed_out {failed}"


# Stats
html = "--html {}".format(snakemake.output.html)
json = "--json {}".format(snakemake.output.json)


shell(
    "(fastp --thread {snakemake.threads} "
    "{extra} "
    "{adapters} "
    "{reads} "
    "{trimmed} "
    "{json} "
    "{html} ) {log}"
)