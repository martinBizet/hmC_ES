The scripts run without any installation as long as the required dependencies and languages are properly installed.
(Tested on Ubuntu 16.04 linux system)

#1° Venn\_Peak (expected time on a normal computer: few minutes)
--------------------------------------------------------------
Tool to overlap several peak lists (require a bed format). Typically used to found reproducible peaks among replicates.

expected output:
Many files (tsv formats) corresponding to intermediate intersections will be generated. To select the file corresponding to the overlap of all samples please use the file referred at the last line of "_selectedFiles.txt"
a pdf file containing the venn diagramm will also be generated.

usage:
./Venn_peaks.R -i=<comma separated peaks files in bed format> -p=<directory where input files are located, if not provided absolute path should be used in -i, use '.' if the inputs are in the same directory as this code> -N=<comma separated names of the peaks> -o=<output directory>

example:
./Venn\_peaks.R -i=GSM3825620\_hMeRIP\_t1608\_ES\_WT\_hmC\_D\_star\_map.bam\_cInput\_peaks\_foldchange\_summits100.bed,GSM3825621\_hMeRIP\_t1703\_ES\_WT\_1\_hmC\_star\_map.bam\_cInput\_peaks\_foldchange\_summits100.bed,GSM3825622\_hMeRIP\_t1703\_ES\_WT\_2\_hmC\_star\_map.bam\_cInput\_peaks\_foldchange\_summits100.bed -p=. -N=R1,R2,R3 -o=Results/

dependencies:
- R (version 3.2.3)
- "VennDiagram" package of R (version 1.6.20)
- the provided python code "compare_bed.py" and its dependencies:
  - python 2 (version 2.7.12)
  - python modules: pybedtools (version 0.7.10), numpy (version 1.11.0), pandas (version 0.17.1)
  - bedtools (version 2.25.0)

#2° annotate-peaks-bed (expected time on a normal computer: few minutes)
-----------------------------------------------------------------------
Tool for RNA peak annotation (require a bed format)

expected output:
a text file (bed format) listing the peaks positions with an extra-column 'annotation' containing notably the transcript(s) ID, the gene(s) name(s) and the transcriptomic region(s) associated to the peak.

usage:
./annotate-peaks-bed <MACS2 peak in bed format> <ouput file> <transcriptome refFlat format> <comma-separated names for the columns available after chromosom, start, stop in the input file bedfile (typically 'name')>

example:
./annotate-peaks-bed GSM3825620\_hMeRIP\_t1608\_ES\_WT\_hmC\_D\_star\_map.bam\_cInput\_peaks\_foldchange\_summits100.bed Results/GSM3825620\_annoted.bed mm9\_RefSeq\_refFlat.txt name

dependencies:
- python 2 (version 2.7.12)
- the provided python code "commandify.py" (has to be in "lib" folder) and its dependencies:
	- python modules: argparse (version 1.2.1)
- the provided python code "config.py" (has to be in "lib" folder)
- the provided python code "base.py" (has to be in "lib" folder)
- the provided python code "ngs.py" (has to be in "lib" folder)

#3° REANNOT\_DNA (expected time on a normal computer: few minutes)
----------------------------------------------------------------
Tool for DNA peak annotation (require a bed format)

expected output:
Many text files (bed format) will be generated. The files with "\_promoNamed.bed" and "\_gbNamed.bed" extensions associate each peak to a transcript at promoter and gene body level respectively.

usage:
set-up the following variables with the appropriate path:
 inputDir= <directory containing MACS2 peaks in bed format>
 outputDir= <output directory>
 BEDTOOLS= <path to bedtools executable>
then, run using bash language.

example of settings:
inputDir="DIP-seq-tERR031632\_less.bam\_peaks\_foldchange\_nocomments\_summits100.bed"
outputDir="Results/"
BEDTOOLS=/usr/bin/bedtools

(Note: the "promo" and "gB" variables are the provided promoter and gene body annotations ("mm9\_promoter.bed", "mm9\_geneBody.bed") but can also be set-up)

dependencies:
- bash (version 4.3.48)
- awk (version 4.1.3)
- bedtools (version 2.25.0)
