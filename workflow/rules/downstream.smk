rule pca:
	input:
		"out/quantification/samples_merged_tpm.csv",
	output:
		pdf="out/downstream/pca_plot.pdf",
		png="out/downstream/pca_plot.png",
		eig="out/downstream/pca_eigenvalues.csv",
		result="out/downstream/pca_metadata.csv",
	threads:1
	log:
		"out/logs/pca.log"
	conda:
		"../envs/rdata.yaml"
	params:
		ncp = config["PCA"]["ncp"],
	script:
		"../scripts/pca.R"

rule deg:
	input:
		"out/quantification/samples_merged_counts.csv"
	output:
		deg_table="out/downstream/deg_table.csv",
	threads:1
	log:
		"out/logs/deg.log"
	conda:
		"../envs/rdata.yaml"
	params:
		normalization = config["edgeR"]["normalization"],
		# contrast = config["edgeR"]["contrast"],
		fdr_threshold = config["edgeR"]["fdr_threshold"],
		log2fc_threshold = config["edgeR"]["logFC_threshold"],
	script:
		"../scripts/deg.R"

# rule wgcna:
# 	input:
# 		"out/quantification/samples_merged_tpm.csv",
# 		"out/quantification/samples_merged_counts.csv"
# 	output:
# 		""
# 	threads:1
# 	log:
# 		"out/logs/wgcna.log"
# 	conda:
# 		"../envs/rdata.yaml"
# 	script:
# 		"../scripts/wgcna.R"