rule post_process:
	input:
		get_post_input(),
	output:
		"out/reports/trim_report.csv",
		"out/reports/hisat2_align_report.csv",
	conda:
		"../envs/reports.yaml",
	log:
		"out/logs/post_process.log",
	params:
		trimdrp =lambda w, output: Path(output[0]).parent,
		hisat2log = "out/logs/hisat2_align",
		clean = config["clean"],
		newGene = config["newGene"]["activate"],
	script:
		"../report/post_process.py"

rule reports:
	input:
		**get_report_input(),
	output:
		"out/reports/report.html",
	params:
		newGene = config["newGene"]["activate"],
		pca = config["PCA"]["activate"],
		deg = config["edgeR"]["activate"],
		# wgcna = config["WGCNA"]["activate"],
	log:
		"out/logs/report.log",
	conda:
		"../envs/reports.yaml",
	script:
		"../report/reports.Rmd"