import pandas as pd
from pathlib import Path
import json
import re
import shutil

# get the parameters
trimdir = snakemake.params.get("trimdrp")
aligndir = snakemake.params.get("hisat2log")
trimout = snakemake.output[0]
alignout = snakemake.output[1]
clean = snakemake.params.get("clean")
newGene = snakemake.params.get("newGene")

# check the parameters
aligndir = Path(aligndir)
alignout = Path(alignout)
trimdir = Path(trimdir)
trimout = Path(trimout)

# define the temp directories 
hisat2_index = Path("out/hisat2_index")
hisat2_mapped = Path("out/mapped")
fastp_report_json = Path("out/reports/json")
fastp_trimmed = Path("out/trimmed")
featureCounts = Path("out/featurecounts")
stringtie = Path("out/stringtie")

assert trimdir.exists(), f"Trimmed reports directory {trimdir} does not exist"
assert aligndir.exists(), f"HISAT2 log directory {aligndir} does not exist"


# process the trimmed report
reportJson = trimdir.rglob("*.json")

reportAll = {}
for report in reportJson:
		# Read the JSON file
		with open(report, 'r') as f:
				data = json.load(f)
				data_filtered = [data["summary"]["before_filtering"],
												 data["summary"]["after_filtering"]]
				data_out = {"total_reads": data_filtered[0]["total_reads"],
										"after_filtering": data_filtered[1]["total_reads"],
										"GC_content": data_filtered[1]["gc_content"]}

				sample = report.stem
				reportAll[sample] = data_out

trim_report = pd.DataFrame.from_dict(reportAll, orient='index')
trim_report.index.name = "Sample"
trim_report.reset_index(inplace=True)

trim_report["PassRate"] = trim_report["after_filtering"] / \
		trim_report["total_reads"] * 100
trim_report["PassRate"] = trim_report["PassRate"].round(2)
trim_report["GC_content"] = trim_report["GC_content"].round(2)
trim_report["PassRate"] = trim_report["PassRate"].astype(str) + "%"
trim_report["GC_content"] = trim_report["GC_content"].astype(str) + "%"
trim_report = trim_report.rename(
		columns={"total_reads": "Total Reads", "after_filtering": "After Filtering",
						 "GC_content": "GC Content", "PassRate": "Pass Rate"})

trim_report = trim_report.sort_values(
		by=["Sample"], ascending=True)


# process the HISAT2 log report
def extract_hisat2_metrics(hisat2_log):

		metrics = {
				"total_reads": 0,
				"mapped_reads": 0,
				"unique_mapping": 0.0
		}

		hisat2_output = Path(hisat2_log).read_text(encoding="utf-8")

		total_reads = re.search(r"(\d+) reads;", hisat2_output).group(1)
		al_con_e_1 = re.search(r"(\d+)\s*\([^)]*\) aligned concordantly exactly 1 time",
													 hisat2_output).group(1)
		al_con_g_1 = re.search(r"(\d+)\s*\([^)]*\) aligned concordantly >1 times",
													 hisat2_output).group(1)
		al_dis_1 = re.search(r"(\d+)\s*\([^)]*\) aligned discordantly 1 time",
												 hisat2_output).group(1)
		al_ex_1 = re.search(r"(\d+)\s*\([^)]*\) aligned exactly 1 time",
												hisat2_output).group(1)
		al_g_1 = re.search(r"(\d+)\s*\([^)]*\) aligned >1 times",
											 hisat2_output).group(1)

		total = int(total_reads) * 2
		mapped = (int(al_con_e_1)*2 + int(al_con_g_1)*2 +
							int(al_dis_1)*2 +
							int(al_ex_1) + int(al_g_1))

		unq_map = (int(al_con_e_1)*2 + int(al_dis_1)*2 +
							 int(al_ex_1))

		metrics["total_reads"] = total
		metrics["mapped_reads"] = mapped
		metrics["unique_mapping"] = unq_map

		return metrics


hisat2_logs = aligndir.glob("*.log")

alignAll = {}
for report in hisat2_logs:
		# Read the HISAT2 log file
		sample = report.stem
		metrics = extract_hisat2_metrics(report)
		alignAll[sample] = metrics

align_report = pd.DataFrame.from_dict(alignAll, orient='index')
align_report.index.name = "Sample"
align_report.reset_index(inplace=True)
align_report["Mapped Rate"] = align_report["mapped_reads"] / \
		align_report["total_reads"] * 100
align_report["Mapped Rate"] = align_report["Mapped Rate"].round(2)
align_report["Mapped Rate"] = align_report["Mapped Rate"].astype(str) + "%"
align_report["Unique Mapped Rate"] = align_report["unique_mapping"] / \
		align_report["total_reads"] * 100
align_report["Unique Mapped Rate"] = align_report["Unique Mapped Rate"].round(
		2)
align_report["Unique Mapped Rate"] = align_report["Unique Mapped Rate"].astype(
		str) + "%"
align_report = align_report.rename(
		columns={"total_reads": "Total Reads", "mapped_reads": "Reads Mapped",
						 "unique_mapping": "Unique Mapped"})

align_report = align_report[["Sample", "Total Reads",
														 "Reads Mapped", "Mapped Rate",
														 "Unique Mapped", "Unique Mapped Rate"]]

align_report = align_report.sort_values(
		by=["Sample"], ascending=True)

# Save two reports to CSV files
trim_report.to_csv(trimout, index=False)
align_report.to_csv(alignout, index=False)

# Clean up the temporary files
if clean:
		shutil.rmtree(hisat2_index)
		shutil.rmtree(fastp_report_json)
		shutil.rmtree(fastp_trimmed)
		shutil.rmtree(featureCounts)
		shutil.rmtree(hisat2_mapped)
		if newGene:
				shutil.rmtree(stringtie)
