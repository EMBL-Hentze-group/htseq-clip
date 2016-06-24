import os
from snakemake.utils import report

# ========================================================================================
# Pipeline by Marko Fritz <marko.fritz@embl.de>, Thomas Schwarzl <schwarzl@embl.de>, Nadia Ashraf <nadia.ashraf@embl.de>
# __author__  = "Marko Fritz, Thomas Schwarzl, Nadia Ashraf"
# __license__ = "EMBL"
# ========================================================================================
# Options and paths

"""
htseq-CLIP workflow

This workflow does standard htSeq-CLIP processing for single-read. 
"""

# ----------------------------------------------------------------------------------------
# CONFIG FILE
# ----------------------------------------------------------------------------------------
configfile: "config.json"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# TOOL PATHS
# ----------------------------------------------------------------------------------------

OUTDIR = config["outdir"]
GTF = config["gtfdir"]
GTFN = config["annotation_file"]
BAMDIR = config["bamdir"]
FASTA = config["fasta"]
FASTAN = config["fasta_name"]
END_PATTERN = config["bam_end"]
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# SAMPLES
# ----------------------------------------------------------------------------------------

# Automatically read in all samples

SAMPLES, = glob_wildcards(expand("{sample_dir}/{{samples}}.{pattern}",pattern = END_PATTERN, sample_dir = config["bamdir"])[0])
# ----------------------------------------------------------------------------------------

SITES = "MS SS ES DEL INS".split()
# ----------------------------------------------------------------------------------------	
# :::::::: ALL :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ----------------------------------------------------------------------------------------	

print("PROCESSING samples %s" % (', '.join(SAMPLES)))
print("PROCESSING sites %s" % (', '.join(SITES)))

rule all:
	input:
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# ::::::: PROCESS FASTA :::::::::
		#expand("../{outdir}/fastq/{fastan}.fastq.gz", fastan=FASTAN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# ::::::: PROCESS GTF :::::::::
		expand("../{outdir}/gtf/{gtfn}.bed.gz", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# ::::::: SORT PROCESSED GTF :::::::::
		expand("../{outdir}/gtf/{gtfn}.sorted.bed.gz", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: PROCESS SLIDING WINDOW :::::::::
		expand("../{outdir}/gtf/{gtfn}.sw.bed.gz", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: SORT SLIDING WINDOW :::::::::
		expand("../{outdir}/gtf/{gtfn}.sw.sorted.bed.gz", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: EXTRACT SITES :::::::::
		expand("../{outdir}/extract/{samples}_{sites}.bed.gz", samples=SAMPLES, outdir=OUTDIR, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: COUNT :::::::::
		expand("../{outdir}/counts/{samples}_{sites}to{gtfn}.count.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: JUNCTION :::::::::
		expand("../{outdir}/junction/{samples}_{sites}to{gtfn}.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: DO SLIDING WINDOW :::::::::
		expand("../{outdir}/slidingWindow/{samples}_{sites}to{gtfn}.sw.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: DO CONVERTION TO DEXSEQ :::::::::
		expand("../{outdir}/dexseq/{samples}_{sites}to{gtfn}.dex.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: BOKEH READ PLOTS :::::::::
		expand("../{outdir}/plots/reads/{samples}_{sites}/{samples}_{sites}to{gtfn}_reads.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: BOKEH COUNT PLOTS :::::::::
		expand("../{outdir}/plots/counts/{samples}_{sites}/{samples}_{sites}to{gtfn}_counts.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: BOKEH JUNCTION PLOTS :::::::::
		expand("../{outdir}/plots/junction/{samples}_{sites}/{samples}_{sites}to{gtfn}_junction.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# ----------------------------------------------------------------------------------------	
# FASTA
# ----------------------------------------------------------------------------------------	

rule process_fasta:
	input:
		fasta=FASTA
	output:
		expand("../{outdir}/fastq/{fastan}.fastq.gz", fastan=FASTAN, outdir=OUTDIR)
	log:
		expand("../logs/{fastan}.log", fastan=FASTAN, outdir=OUTDIR)
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py genomeToReads -i {input.fasta} -o {output} -x 42 2> {log}"
		
				
# ----------------------------------------------------------------------------------------	
# GTF
# ----------------------------------------------------------------------------------------	
				
rule process:
	input:
		gtf=GTF
	output:
		expand("../{outdir}/gtf/{gtfn}.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.log", gtfn=GTFN, outdir=OUTDIR)
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py process -g {input.gtf} -t gene_type -o {output} 2> {log}"
		
rule sort_gtf:
	input:
		expand("../{outdir}/gtf/{gtfn}.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	output:
		expand("../{outdir}/gtf/{gtfn}.sorted.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.sorted.log", gtfn=GTFN , outdir=OUTDIR)
	shell:
		"zcat {input} | sort -k1,1 -k2,2n | gzip > {output} 2> {log}"
		
				
rule sliding_window:
	input:
		expand("../{outdir}/gtf/{gtfn}.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	output:
		expand("../{outdir}/gtf/{gtfn}.sw.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.sw.log", gtfn=GTFN, outdir=OUTDIR)
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py slidingWindow -i {input} -o {output} -w 50 -s 20 2> {log}"
		
rule sort_sliding_window:
	input:
		expand("../{outdir}/gtf/{gtfn}.sw.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	output:
		expand("../{outdir}/gtf/{gtfn}.sw.sorted.bed.gz", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.sw.sorted.log", gtfn=GTFN , outdir=OUTDIR)
	shell:
		"zcat {input} | sort -k1,1 -k2,2n | gzip > {output} 2> {log}"
		
				
# ----------------------------------------------------------------------------------------	
# EXTRACT
# ----------------------------------------------------------------------------------------	

#START SITES
rule extract_SS:
	input:
		expand("../{outdir}/extract/{samples}_SS.temporary.bed.gz", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_SS:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_SS.temporary.bed.gz")
	log:
		"../logs/{sample}_SS.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py extract -i {input} -o {output} -c s 2> {log}"

#MIDDLE SITES
rule extract_MS:
	input:
		expand("../{outdir}/extract/{samples}_MS.temporary.bed.gz", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_MS:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_MS.temporary.bed.gz")
	log:
		"../logs/{sample}_MS.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py extract -i {input} -o {output} -c m 2> {log}"

#END SITES
rule extract_ES:
	input:
		expand("../{outdir}/extract/{samples}_ES.temporary.bed.gz", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_ES:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_ES.temporary.bed.gz")
	log:
		"../logs/{sample}_ES.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py extract -i {input} -o {output} -c e 2> {log}"

#DELETION SITES
rule extract_DEL:
	input:
		expand("../{outdir}/extract/{samples}_DEL.temporary.bed.gz", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_DEL:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_DEL.temporary.bed.gz")
	log:
		"../logs/{sample}_DEL.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py extract -i {input} -o {output} -c d 2> {log}"

#INSERTION SITES
rule extract_INS:
	input:
		expand("../{outdir}/extract/{samples}_INS.temporary.bed.gz", samples=SAMPLES, outdir=OUTDIR)

		
rule do_extract_INS:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_INS.temporary.bed.gz")
	log:
		"../logs/{sample}_INS.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py extract -i {input} -o {output} -c i 2> {log}"
		
		
		
# ----------------------------------------------------------------------------------------	
# EXTRACT SITES
# ----------------------------------------------------------------------------------------	
rule sort:
	input:
		expand("../{outdir}/extract/{samples}_{sites}.bed.gz", samples=SAMPLES, outdir=OUTDIR, sites=SITES)

rule do_sort:		
	input:
		"../{OUTDIR}/extract/{sample}_{site}.temporary.bed.gz"
	output:
		"../{OUTDIR}/extract/{sample}_{site,\w+}.bed.gz"
	log:
		"{../logs/sample}_{site}.extract.sort.log"
	shell:
		"zcat {input} | sort -k1,1 -k2,2n | gzip > {output} 2> {log}"
		
# ----------------------------------------------------------------------------------------
# JUNCTION
# ----------------------------------------------------------------------------------------

rule junction:
	input:
		expand("../{outdir}/junction/{samples}_{sites}to{gtfn}.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_junction:
	input:
		bed="../{OUTDIR}/extract/{sample}_{site}.bed.gz",
		gtf="../{OUTDIR}/gtf/{gtfn}.sorted.bed.gz"
	output:
		"../{OUTDIR}/junction/{sample}_{site}to{gtfn}.txt.gz"
	log:
		"../logs/{sample}_{site}.junction.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py junction -i {input.bed} -f {input.gtf} -o {output} 2> {log}"
		
		
# ----------------------------------------------------------------------------------------
# COUNT
# ----------------------------------------------------------------------------------------

rule count:
	input:
		expand("../{outdir}/counts/{samples}_{sites}to{gtfn}.count.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_count:
	input:
		bed="../{OUTDIR}/extract/{sample}_{site}.bed.gz",
		gtf="../{OUTDIR}/gtf/{gtfn}.sorted.bed.gz"
	output:
		"../{OUTDIR}/counts/{sample}_{site}to{gtfn}.count.txt.gz"
	log:
		"../logs/{sample}_{site}.count.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py count -i {input.bed} -f {input.gtf} -o {output} -c o 2> {log}"
		
		
# ----------------------------------------------------------------------------------------
# SLIDING WINDOW
# ----------------------------------------------------------------------------------------

rule sw:
	input:
		expand("../{outdir}/slidingWindow/{samples}_{sites}to{gtfn}.sw.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_sw:
	input:
		bed="../{OUTDIR}/extract/{sample}_{site}.bed.gz",
		gtf="../{OUTDIR}/gtf/{gtfn}.sw.sorted.bed.gz"
	output:
		"../{OUTDIR}/slidingWindow/{sample}_{site}to{gtfn}.sw.txt.gz"
	log:
		"../logs/{sample}_{site}.sw.count.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py countSlidingWindow -i {input.bed} -f {input.gtf} -o {output} 2> {log}"
				
# ----------------------------------------------------------------------------------------
# DO CONVERTION TO DEXSEQ
# ----------------------------------------------------------------------------------------

rule dexseq:
	input:
		expand("../{outdir}/dexseq/{samples}_{sites}to{gtfn}.dex.txt.gz", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_dexseq:
	input:
		"../{OUTDIR}/slidingWindow/{sample}_{site}to{gtfn}.sw.txt.gz"
	output:
		"../{OUTDIR}/dexseq/{sample}_{site}to{gtfn}.dex.txt.gz"
	log:
		"../logs/{sample}_{site}.dex.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py toDexSeq -i {input} -o {output} 2> {log}"
		
		
# ----------------------------------------------------------------------------------------	
# BOKEH READ PLOTS
# ----------------------------------------------------------------------------------------	

rule rplot:
	input:
		expand("../{outdir}/plots/reads/{samples}_{sites}/{samples}_{sites}to{gtfn}_reads.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		
rule do_rplot:
	input:
		"../{OUTDIR}/extract/{sample}_{site}.bed.gz"
	output:
		"../{OUTDIR}/plots/reads/{sample}_{site}/{sample}_{site}to{gtfn}_reads.html"
	log:
		"../logs/{sample}_{site}.rplot.log"
	shell:
		"PYTHONPATH="" "
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py plot -i {input} -o {output} -c r 2> {log}"
		
# ----------------------------------------------------------------------------------------	
# BOKEH COUNT PLOTS
# ----------------------------------------------------------------------------------------	

rule cplot:
	input:
		expand("../{outdir}/plots/counts/{samples}_{sites}/{samples}_{sites}to{gtfn}_counts.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		
rule do_cplot:
	input:
		"../{OUTDIR}/counts/{sample}_{site}to{gtfn}.count.txt.gz"
	output:
		"../{OUTDIR}/plots/counts/{sample}_{site}/{sample}_{site}to{gtfn}_counts.html"
	log:
		"../logs/{sample}_{site}.cplot.log"
	shell:
		"PYTHONPATH="" "
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py plot -i {input} -o {output} -c c 2> {log}"
		
# ----------------------------------------------------------------------------------------	
# BOKEH JUNCTION PLOTS
# ----------------------------------------------------------------------------------------	

rule jplot:
	input:
		expand("../{outdir}/plots/junction/{samples}_{sites}/{samples}_{sites}to{gtfn}_junction.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		
rule do_jplot:
	input:
		"../{OUTDIR}/junction/{sample}_{site}to{gtfn}.txt.gz"
	output:
		"../{OUTDIR}/plots/junction/{sample}_{site}/{sample}_{site}to{gtfn}_junction.html"
	log:
		"../logs/{sample}_{site}.jplot.log"
	shell:
		"PYTHONPATH="" "
		"python /g/hentze/projects/Software/htseq-clip/clip/clip.py plot -i {input} -o {output} -c j 2> {log}"
		