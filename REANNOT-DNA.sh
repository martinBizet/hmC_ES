#need be set-up (see comment below for an example)
inputDir="<directory containing MACS2 peaks in bed format>"
outputDir="<output directory (including terminal /)>"
BEDTOOLS="<path to bedtools executable>"

#EXAMPLE:
#--------
#inputDir="DIP-seq-tERR031632_less.bam_peaks_foldchange_nocomments_summits100.bed"
#outputDir="Results/"
#BEDTOOLS=/usr/bin/bedtools

mkdir ${outputDir}

promo="Data/mm9_promoter.bed"
gB="Data/mm9_geneBody.bed"

header_sup="name"

for inputFile in $( ls ${inputDir} )
do
	outNam=$( basename ${inputFile} )

	##annotate
	function promo_bed(){
		awk 'FNR > 1 {if ($4 != 0) { print $1"\t"$2"\t"$3 }}' ${1} > ${2}
	}

	function gb_bed(){
		awk 'FNR > 1 {if ($5 != 0) { print $1"\t"$2"\t"$3 }}' ${1} > ${2}
	}

	# Annotate peaks
	outNam=${outputDir}${outNam}
	awk 'BEGIN{OFS=FS="\t"}{if($2 < 0) {print $1,int($3/2),int($3/2)} else {print $1,int(($2+$3)/2),int(($2+$3)/2)}}' ${inputFile}  > ${outNam}.bed
	cmd=${BEDTOOLS}" annotate -i "${outNam}.bed" -counts -files "${promo}" "${gB}" > "${outNam}_annotated.bed

	echo ${cmd}
	eval ${cmd}

	# Annotate promoter: gene name
	promo_bed ${outNam}_annotated.bed ${outNam}_promo.bed
	cmd=${BEDTOOLS}" intersect -a "${outNam}_promo.bed" -b "${promo}" -wo > "${outNam}_promoNamed.bed
	echo ${cmd}
	eval ${cmd}

	# Annotate gene body: gene name
	gb_bed ${outNam}_annotated.bed ${outNam}_gb.bed
	cmd=${BEDTOOLS}" intersect -a "${outNam}_gb.bed" -b "${gB}" -wo > "${outNam}_gbNamed.bed
	echo ${cmd}
	eval ${cmd}

done
