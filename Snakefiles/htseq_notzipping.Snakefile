import os
from snakemake.utils import report

# ========================================================================================
# Pipeline by Marko Fritz <marko.fritz@embl.de>, Thomas Schwarzl <schwarzl@embl.de>
# __author__  = "Marko Fritz, Thomas Schwarzl"
# __license__ = "EMBL"
# ========================================================================================
# Options and paths

"""
htseq-CLIP workflow

This workflow does standard htSeq-CLIP processing for single-read or paired-end. 
"""

# ----------------------------------------------------------------------------------------
# TOOL PATHS
# ----------------------------------------------------------------------------------------

GTF = "/g/hentze/projects/Software/htseq-clip/Genomes/Homo_sapiens.GRCh37.82.gtf"
GTFN = "Homo_sapiens.GRCh37.82"

#dir for output
OUTDIR = "output"
BAMDIR = "/g/hentze/projects/Software/htseq-clip/GRCh37_bam"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# SAMPLES
# ----------------------------------------------------------------------------------------

# Automatically read in all samples

SAMPLES, = glob_wildcards("/g/hentze/projects/Software/htseq-clip/GRCh37_bam/{samples}.bam")
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
		# ::::::: PROCESS GTF :::::::::
		expand("../{outdir}/gtf/{gtfn}.bed", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# ::::::: SORT PROCESSED GTF :::::::::
		expand("../{outdir}/gtf/{gtfn}.sorted.bed", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: PROCESS SLIDING WINDOW :::::::::
		expand("../{outdir}/gtf/{gtfn}.sw.bed", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: SORT SLIDING WINDOW :::::::::
		expand("../{outdir}/gtf/{gtfn}.sw.sorted.bed", gtfn=GTFN, outdir=OUTDIR),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: EXTRACT SITES :::::::::
		expand("../{outdir}/extract/{samples}_{sites}.bed", samples=SAMPLES, outdir=OUTDIR, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: COUNT :::::::::
		expand("../{outdir}/counts/{samples}_{sites}to{gtfn}.count.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: JUNCTION :::::::::
		expand("../{outdir}/junction/{samples}_{sites}to{gtfn}.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: DO SLIDING WINDOW :::::::::
		expand("../{outdir}/slidingWindow/{samples}_{sites}to{gtfn}.sw.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: DO CONVERTION TO DEXSEQ :::::::::
		expand("../{outdir}/dexseq/{samples}_{sites}to{gtfn}.dex.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES),
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: BOKEH JUNCTION PLOTS :::::::::
		#expand("../{outdir}/plots/junction/{samples}_{sites}/{samples}_{sites}to{gtfn}_junction.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# ::::::: PLOT :::::::::
		#expand("../{outdir}/plots/{samples}/{samples}_{sites}to{gtfn}_exonStart.pdf", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		
				
				
# ----------------------------------------------------------------------------------------	
# GTF
# ----------------------------------------------------------------------------------------	
				
rule process:
	input:
		gtf=GTF
	output:
		expand("../{outdir}/gtf/{gtfn}.bed", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.log", gtfn=GTFN, outdir=OUTDIR)
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py process -g {input.gtf} -o {output} 2> {log}"
		
rule sort_gtf:
	input:
		expand("../{outdir}/gtf/{gtfn}.bed", gtfn=GTFN, outdir=OUTDIR)
	output:
		expand("../{outdir}/gtf/{gtfn}.sorted.bed", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.sorted.log", gtfn=GTFN , outdir=OUTDIR)
	shell:
		"sort -k1,1 -k2,2n {input} > {output} 2> {log}"
		
				
rule sliding_window:
	input:
		expand("../{outdir}/gtf/{gtfn}.bed", gtfn=GTFN, outdir=OUTDIR)
	output:
		expand("../{outdir}/gtf/{gtfn}.sw.bed", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.sw.log", gtfn=GTFN, outdir=OUTDIR)
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py slidingWindow -i {input} -o {output} -w 50 -s 20 2> {log}"
		
rule sort_sliding_window:
	input:
		expand("../{outdir}/gtf/{gtfn}.sw.bed", gtfn=GTFN, outdir=OUTDIR)
	output:
		expand("../{outdir}/gtf/{gtfn}.sw.sorted.bed", gtfn=GTFN, outdir=OUTDIR)
	log:
		expand("../logs/{gtfn}.sw.sorted.log", gtfn=GTFN , outdir=OUTDIR)
	shell:
		"sort -k1,1 -k2,2n {input} > {output} 2> {log}"
		
				
# ----------------------------------------------------------------------------------------	
# EXTRACT
# ----------------------------------------------------------------------------------------	

#START SITES
rule extract_SS:
	input:
		expand("../{outdir}/extract/{samples}_SS.temporary.bed", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_SS:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_SS.temporary.bed")
	log:
		"../logs/{sample}_SS.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py extract -i {input} -o {output} -c s 2> {log}"

#MIDDLE SITES
rule extract_MS:
	input:
		expand("../{outdir}/extract/{samples}_MS.temporary.bed", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_MS:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_MS.temporary.bed")
	log:
		"../logs/{sample}_MS.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py extract -i {input} -o {output} -c m 2> {log}"

#END SITES
rule extract_ES:
	input:
		expand("../{outdir}/extract/{samples}_ES.temporary.bed", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_ES:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_ES.temporary.bed")
	log:
		"../logs/{sample}_ES.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py extract -i {input} -o {output} -c e 2> {log}"

#DELETION SITES
rule extract_DEL:
	input:
		expand("../{outdir}/extract/{samples}_DEL.temporary.bed", samples=SAMPLES, outdir=OUTDIR)
		
rule do_extract_DEL:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_DEL.temporary.bed")
	log:
		"../logs/{sample}_DEL.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py extract -i {input} -o {output} -c d 2> {log}"

#INSERTION SITES
rule extract_INS:
	input:
		expand("../{outdir}/extract/{samples}_INS.temporary.bed", samples=SAMPLES, outdir=OUTDIR)

		
rule do_extract_INS:
	input:
		expand("{bamdir}/{{sample}}.bam", bamdir=BAMDIR)
	output:
		temp("../{OUTDIR}/extract/{sample}_INS.temporary.bed")
	log:
		"../logs/{sample}_INS.extract.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py extract -i {input} -o {output} -c i 2> {log}"
		
		
		
# ----------------------------------------------------------------------------------------	
# EXTRACT SITES
# ----------------------------------------------------------------------------------------	
rule sort:
	input:
		expand("../{outdir}/extract/{samples}_{sites}.bed", samples=SAMPLES, outdir=OUTDIR, sites=SITES)

rule do_sort:		
	input:
		"../{OUTDIR}/extract/{sample}_{site}.temporary.bed"
	output:
		"../{OUTDIR}/extract/{sample}_{site,\w+}.bed"
	log:
		"{../logs/sample}_{site}.extract.sort.log"
	shell:
		"sort -k1,1 -k2,2n {input} > {output} 2> {log}"
		
		
# ----------------------------------------------------------------------------------------
# JUNCTION
# ----------------------------------------------------------------------------------------

rule junction:
	input:
		expand("../{outdir}/junction/{samples}_{sites}to{gtfn}.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_junction:
	input:
		bed="../{OUTDIR}/extract/{sample}_{site}.bed",
		gtf="../{OUTDIR}/gtf/{gtfn}.sorted.bed"
	output:
		"../{OUTDIR}/junction/{sample}_{site}to{gtfn}.txt"
	log:
		"../logs/{sample}_{site}.junction.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py junction -i {input.bed} -f {input.gtf} -o {output} 2> {log}"
		
		
# ----------------------------------------------------------------------------------------
# COUNT
# ----------------------------------------------------------------------------------------

rule count:
	input:
		expand("../{outdir}/counts/{samples}_{sites}to{gtfn}.count.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_count:
	input:
		bed="../{OUTDIR}/extract/{sample}_{site}.bed",
		gtf="../{OUTDIR}/gtf/{gtfn}.sorted.bed"
	output:
		"../{OUTDIR}/counts/{sample}_{site}to{gtfn}.count.txt"
	log:
		"../logs/{sample}_{site}.count.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py count -i {input.bed} -f {input.gtf} -o {output} c- o 2> {log}"
		
		
# ----------------------------------------------------------------------------------------
# SLIDING WINDOW
# ----------------------------------------------------------------------------------------

rule sw:
	input:
		expand("../{outdir}/slidingWindow/{samples}_{sites}to{gtfn}.sw.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_sw:
	input:
		bed="../{OUTDIR}/extract/{sample}_{site}.bed",
		gtf="../{OUTDIR}/gtf/{gtfn}.sw.sorted.bed"
	output:
		"../{OUTDIR}/slidingWindow/{sample}_{site}to{gtfn}.sw.txt"
	log:
		"../logs/{sample}_{site}.sw.count.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py countSW -i {input.bed} -f {input.gtf} -o {output} 2> {log}"
		
		
		
# ----------------------------------------------------------------------------------------
# DO CONVERTION TO DEXSEQ
# ----------------------------------------------------------------------------------------

rule dexseq:
	input:
		expand("../{outdir}/dexseq/{samples}_{sites}to{gtfn}.dex.txt", samples=SAMPLES, outdir=OUTDIR, gtfn=GTFN, sites=SITES)
		
rule do_dexseq:
	input:
		"../{OUTDIR}/slidingWindow/{sample}_{site}to{gtfn}.sw.txt"
	output:
		"../{OUTDIR}/dexseq/{sample}_{site}to{gtfn}.dex.txt"
	log:
		"../logs/{sample}_{site}.dex.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py toDexSeq -i {input} -o {output} 2> {log}"
				
# ----------------------------------------------------------------------------------------	
# BOKEH JUNCTION PLOTS
# ----------------------------------------------------------------------------------------	

rule jplot:
	input:
		expand("../{outdir}/plots/junction/{samples}_{sites}/{samples}_{sites}to{gtfn}_junction.html", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)
		
rule do_jplot:
	input:
		"../{OUTDIR}/junction/{sample}_{site}to{gtfn}.txt"
	output:
		"../{OUTDIR}/plots/junction/{sample}_{site}/{sample}_{site}to{gtfn}_junction.html"
	log:
		"../logs/{sample}_{site}.jplot.log"
	shell:
		"python /g/hentze/projects/Software/htseq-clip/htseq-clip/htseqClip.py bokeh -i {input} -o {output} -c j 2> {log}"

# ----------------------------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------------------------
		
rule plot:
	input:
		expand("../{outdir}/plots/{samples}/{samples}_{sites}to{gtfn}_exonStart.pdf", samples=SAMPLES, gtfn=GTFN, sites=SITES, outdir=OUTDIR)

rule do_plot:
	input:
		expand("../{outdir}/junction/{{sample}}_{{site}}to{{gtfn}}.txt", outdir=OUTDIR)
	output:
		file=expand("../{outdir}/plots/{{sample}}/{{sample}}_{{site}}to{{gtfn}}_exonStart.pdf", outdir=OUTDIR)
		#dir=expand("../{outdir}/plots/{{sample}}/", outdir=OUTDIR)
	log:
		"../logs/{sample}_{site}.plot.log"
	shell:
		"Rscript /g/hentze/projects/Software/htseq-clip/htseq-clip/junction.r {input} $(dirname {output}) 2> {log}"
		
rule mkplotdir:
	output:
		"../plots/{sample}/"
	shell:
		"mkdir -p {output}"