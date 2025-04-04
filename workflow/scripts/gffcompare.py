from snakemake.shell import shell
from pathlib import Path
import pandas as pd

st_gtf = snakemake.input.get("st_transcripts")
anno = snakemake.input.get("anno")
out = snakemake.output[0]
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output = out.split(".")[0]

shell("gffcompare -G -r {anno} "
      "-o {output} "
      "{st_gtf}")


def newGeneGtf(gffout: Path):
    columns = ['seqname', 'source', 'feature', 'start', 'end',
               'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gffout,  sep='\t', comment='#',
                     names=columns, header=None, low_memory=False)

    mask_feature = df['feature'] == 'transcript'
    mask_class = df['attribute'].str.contains(r'class_code\s+"u"', regex=True)
    filtered_df = df[mask_feature & mask_class]

    transcript_ids = filtered_df['attribute'].str.extract(
        r'transcript_id "([^"]+)"', expand=False).dropna()

    return transcript_ids


out_gtf = Path(out)
novel_ids = Path(out_gtf.parent, "novel_ids.txt")
novel_gtf = Path(out_gtf.parent, "novel_with_exons.gtf")
new_gene_ids = newGeneGtf(out_gtf)
new_gene_ids.to_csv(novel_ids,  index=False,
                    header=False, sep='\t')

if new_gene_ids.empty:
    shell("echo 'No novel genes found'")
else:
    shell("grep -wFf {novel_ids} {out_gtf} > {novel_gtf}")
