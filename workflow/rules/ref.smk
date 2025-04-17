rule hisat2_index:
		input:
				fasta=config["ref"]["genome"],
		output:
				directory("out/hisat2_index"),
		params:
				extra="--seed 42",
		priority: 50,
		cache: True,
		conda:
				"../envs/hisat2.yaml",
		log:
				"out/logs/hisat2_index.log",
		threads: 10
		script:
				"../scripts/hisat2_index.py"

rule download_sra:
		input:
				config["SRA"]["acc_list"]
		output:
				directory(config["SRA"]["output_dir"]),
		conda:
				"../envs/sratools.yaml"
		log:
				"out/logs/sra_download.log",
		params:
				extra="--skip-technical",
		threads: 8	# defaults to 6
		script:
				"../scripts/sra_download.py"