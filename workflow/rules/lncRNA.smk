rule lnc_init:
		input:
				"out/stringtie/merged_assembly.gtf"
		output:
				"out/lncRNA/candidate_transcript.txt"
		conda:
				"../envs/lnctools.yaml"
		shell:
				"""
				grep 'transcript_id "MSTRG' {input} > candidate_transcript.gtf && \
				gffread -w candidate_transcript.fasta -g genome.fasta candidate_transcript.gtf && \
				grep '>' candidate_transcript.fasta | awk '{print \$1}' | sed 's/>//g' | sort -u > candidate_transcript.txt
				"""

rule lnc_identification:
		input:
				"out/lncRNA/candidate_transcript.txt"
		output:
				"out/lncRNA/lncRNA_candidates.txt"
		shell:
				"""
				python3 lnc_identification.py {input} {output}
				"""

rule lnc_classifying:
		input:
				"out/lncRNA/lncRNA_candidates.txt"
		output:
				"out/lncRNA/lncRNA_classified.txt"
		shell:
				"""
				python3 lnc_classifying.py {input} {output}
				"""

rule lnc_TE_derived:
		input:
				"out/lncRNA/lncRNA_classified.txt"
		output:
				"out/lncRNA/lncRNA_TE_derived.txt"
		shell:
				"""
				python3 lnc_TE_derived.py {input} {output}
				"""