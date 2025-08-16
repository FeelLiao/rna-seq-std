__author__ = "Feel Liao"
__copyright__ = "Copyright 2025, Feel Liao"
__email__ = "feel2027@outlook.com"
__license__ = "MIT"

from pathlib import Path
import pandas as pd
import tempfile
import os
import threading
import queue
from snakemake.shell import shell


log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra = snakemake.params.get("extra", "")
acc_list_file = snakemake.input[0]
fastq_path = snakemake.output[0]
acc_list_file = Path(acc_list_file)
fastq_path = Path(fastq_path)
threads = snakemake.threads

assert acc_list_file.exists(), f"File not found: {acc_list_file}"

fastq_path.mkdir(parents=True, exist_ok=True)

# Download SRA data from ncbi
# TODO: Add parallel download option: When one SRA file is downloaded,
# it is converted to fastq format immediately.
with tempfile.TemporaryDirectory(prefix="sra_download_") as tmpdir:
    tmp_path = Path(tmpdir)

    with open(acc_list_file, 'r') as f:
        accessions = [line.strip() for line in f if line.strip()]

    dump_threads = max(threads, os.cpu_count()//4)
    consumers = max(1, dump_threads // 2)

    q = queue.Queue(maxsize=consumers * 2)
    errors = []

    def find_sra(acc: str) -> Path | None:
        matchs = list(tmp_path.rglob(f"{acc}*.sra"))
        return matchs[0] if matchs else None

    def download_sra():
        for acc in accessions:
            try:
                shell("prefetch -O {tmp_path} "
                      "{acc} "
                      "{log} ")
                sra = find_sra(acc)
                if not sra or not sra.exists():
                    msg = f"Failed to download SRA file for accession: {acc}"
                    errors.append(msg)
                    continue
                q.put(sra)
            except Exception as e:
                errors.append(
                    f"Error downloading SRA file for accession {acc}: {str(e)}")

        for _ in range(consumers):
            q.put(None)

    def convert_sra_to_fastq(worker_id: int):
        while True:
            sra = q.get()
            if sra is None:
                break
            try:
                shell("parallel-fastq-dump -O {fastq_path} "
                      "-t {dump_threads} "
                      "--split-files -s {sra} "
                      "{log}")
                os.remove(sra)
            except Exception as e:
                errors.append(f"Error fastq-dump failed for {sra}: {e}")

    try:
        prod_t = threading.Thread(target=download_sra,name="Producer",daemon=True)
        cons_ts = [
            threading.Thread(target=convert_sra_to_fastq, args=(i,), name=f"Consumer-{i}", daemon=True)
            for i in range(consumers)
        ]
        prod_t.start()
        for t in cons_ts:
            t.start()
        prod_t.join()
        for t in cons_ts:
            t.join()

    except Exception as e:
        errors.append(f"Error occurred in processing the sra: {e}")


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
