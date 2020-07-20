#!/usr/bin/env Rscript
argu <- commandArgs(TRUE)

###-------------
##	ARGUMENTS & HELP
#---------------

##Defaults
print.help <- FALSE
peaks <- NULL; nam <- NULL; path.peaks <- '';
sel <- 'percentage'
n <- 53962 #arbitrarily set to 1 peaks by transcripts in UCSC hg19

##Treat argu
argu <- strsplit(argu, split='=')
argu.1 <- argu[which(sapply(argu, length) == 1)]
if(length(argu.1) > 0) {
	for(el in argu.1) {
		if(is.element(el[1], c('-h', '--help'))) {
			print.help <- TRUE
		}
	}
}
if(print.help) {
	cat('Venn_peaks.R draw Venn for 2 to 5 peak files comparison by calling compare_bed.py for pairwise comparisons.\n
		\t-i=|--inputFiles=  comma separated peaks files (+path if not specified in -p). [mandatory]\n
		\t-p=|--pathToFiles=  directory where input files are located, if not provided absolute path should be used in -i. [default= root]\n
		\t-N=|--names=  comma separated names of the peaks (should be the same size as -i). [default= -i]\n
		\t-o=|--outputDir=  directory where to put results and figures. [mandatory]\n
		\t-s=|--selectionMethod=  way to select the best overlap. can be \'percentage\' (highest percentage of peaks overlapping) or \'absolute\' (highest absolute number of peaks overlapping) [default= \'percentage\']\n
		\t-n=|--setSize=  assumption about the total number of possible peaks for hypergeometric computation (Only for pairwise comparison). Either a number either \'hg\' or \'mm\' for number of UCSC transcripts from human or muse respectively [default= \'hg\']\n
		\t-f=|--multiply=  multiplying factor for -n. Should be set to 5 for m6A experiment if -n is not modified [default= 1]\n
		\t-h|--help  print this help\n
		Example\n
		-------\n
		Rscript Venn_peaks.R -i=peaksFile1,peaksFile2 -p=</path/to/peaks/files/> -N=p1,p2 -o=</path/to/output/directory/> -s=percentage -n=5362 -f=1')
} else {
	argu <- argu[which(sapply(argu, length) > 1)]
	argu.keys <- sapply(argu, '[[', 1)
	argu.vals <- sapply(argu, '[[', 2)
	if(length(intersect(c('-i', '-inputFiles'), argu.keys)) > 0) {
		peaks <- argu.vals[which(is.element(argu.keys, c('-i', '--inputFile')))]
		peaks <- strsplit(peaks, split=',')[[1]]
	} else {
		stop('Please provide an input file (-i or --inputFile argument)')
	}
	if(length(intersect(c('-p', '--pathToFiles'), argu.keys)) > 0) { path.peaks <- argu.vals[which(is.element(argu.keys, c('-p', '--pathToFiles')))] }
	if(length(intersect(c('-N', '--names'), argu.keys)) > 0) {
		nam <- argu.vals[which(is.element(argu.keys, c('-N', '--names')))]
		nam <- strsplit(nam, split=',')[[1]]
	} else {
		nam <- strsplit(peaks, split='/')
		for(i in 1:length(nam)) { nam[[i]] <- nam[[i]][length(nam)] }
		nam <- unlist(nam)
	}
	if(length(intersect(c('-o', '--outputDir'), argu.keys)) > 0) {
		outdir <- argu.vals[which(is.element(argu.keys, c('-o', '--outputDir')))]
	} else {
		stop('Please provide an output file (-o or --outputDir argument)')
	}
	if(length(intersect(c('-s', '--selectionMethod'), argu.keys)) > 0) { sel <- argu.vals[which(is.element(argu.keys, c('-s', '--selectionMethod')))] }
	if(length(intersect(c('-n', '--setSize'), argu.keys)) > 0) {
		n <- argu.vals[which(is.element(argu.keys, c('-n', '--setSize')))]
		if(n == 'hg') {
			n <- 53962
		} else if (n == 'mm') {
			n <- 31314
		} else {
			n <- as.numeric(n)
		}
	}
	if(length(intersect(c('-f', '--multiply'), argu.keys)) > 0) { n <- as.numeric( argu.vals[which(is.element(argu.keys, c('-f', '-multiply')))] )*n }
	
	###-------------
	##	FUNCTIONS
	#---------------

	require('VennDiagram')
	modif.draw <- function(draw, add) {
		f <- function(el){return( data.class(el[[1]]) == "character")}
		idx <- which(sapply(draw, f))
		for(i in idx) {
			if(is.element(draw[[i]][[1]], names(add))) {
				draw[[i]][[1]] <- paste(draw[[i]][[1]], ' (', add[ draw[[i]][[1]] ], '%)', sep='')
			}
		}
		print
		return(draw)
	}	
	extract.counts <- function(prefix, r1, r2) {
		r <- c(r1, r2)
		nam <- c(paste(prefix, '_specific_', r[1], '.tsv', sep=''), paste(prefix, '_specific_', r[2], '.tsv', sep=''), NA, NA)
		r1.spec <- length(readLines(con=nam[1]))
		r2.spec <- length(readLines(con=nam[2]))
		inter <- rep(NA, 2)
		for(i in 1:length(inter)) {
			nam[i+2] <- paste(prefix, '_intersect_', r[i], '.tsv', sep='')
			inter[i] <- length(readLines(con=nam[i+2]))
		}
		res <- c(r1.spec, r2.spec, inter)
		names(res) <- nam
		return(res)
	}
	get.best <- function(counts, method='percentage') {
		r1.spec <- counts[1]; r2.spec <- counts[2];
		if(method == 'percentage') { #best = highest percentage of overlap (overlap from the small set is favored)
			res <- c(-Inf, -Inf, counts[3]/sum(counts[c(1,3)]), counts[4]/sum(counts[c(1,4)]))
		} else if(method == 'absolute') { #best = highest absolute overlap (overlap from the large set is favored)
			res <- c(-Inf, -Inf, counts[3], counts[4])
		} else if(method == 'minimal') { #best = smallest overlap (specific sides are favored)
			res <- c(-Inf, -Inf, -1*counts[3], -1*counts[4])
		} else {
			stop('method not available')
		}
		return(which.max(res))
	}
	get.label <- function(nam, peaks) {
		nam <- unique(strsplit(nam, split='-vs-')[[1]])
		label <- which(is.element(peaks, nam))
		label <- paste('n', paste(label, collapse=''), sep='')
		return(label)
	}
	check.conflict <- function(inter, label) {
		f <- function(v1, v2) { return(sum(is.element(v1, v2)) == length(v1)) }
		label <- strsplit(label, split='')[[1]]
		nam <- strsplit(names(inter), split='')
		inter <- inter[which(sapply(nam, f, v2=label))]
		return(min(inter))
	}
	
	###-------------------------------------
	##	MAIN
	#---------------------------------------
	
	peaks <- file.path(path.peaks, peaks)
	names(peaks) <- nam
	print('peaks:')
	print(peaks)
	if(length(peaks) == 2) { # 2 is treated independently because it is different: the two intersect are drawn and the pvalue is computed
		prefix <-  file.path(outdir, paste(names(peaks)[1], names(peaks)[2], sep='-vs-'))
		#run compar_bed
		cmd=paste('./compare_bed.py', '-f', peaks[1], '-F', peaks[2], '-o', prefix, '-n', names(peaks)[1], '-N', names(peaks)[2])
		print(cmd)
		system(cmd)
		#counts
		r <- rep(names(peaks), 2)		
		counts <- extract.counts(prefix=prefix, r1=r[1], r2=r[2])
		#draw venn	
		draws <- list() #1 and 2 will remain NULL
		pvs <- rep(NA, length(counts)) #1 and 2 will remain NA
		pdf(paste(prefix, 'pdf', sep='.')) #Fake plot just to extract draws and directly overwrite
			for(i in 3:length(counts)) {
				r1.spec <- counts[1]; r2.spec <- counts[2]; inter <- counts[i]; 
				area1 <- r1.spec + inter
				area2 <- r2.spec + inter
				cross.area <- inter
				if(i == 3) {
					perc <- round(c(r1.spec, inter)*100/area1, 2); names(perc) <-  c(r1.spec, inter);
				} else if (i == 4) {
					perc <- round(c(r2.spec, inter)*100/area2, 2); names(perc) <-  c(r2.spec, inter); 
				}				
				#hypergeometrical pvalue
				pvs[i] <- phyper(q=cross.area-1, m=area1, n=n-area1, k=area2, lower.tail=FALSE)
				draws[[i]] <- modif.draw(draw.pairwise.venn(area1=area1, area2=area2, cross.area=cross.area, category=c(r[1], r[2]), col=c('royalblue4', 'grey'), lwd=2, ext.length=0.9), add=perc)
			}
		dev.off()
		pdf(paste(prefix, 'pdf', sep='.')) #re-draw venn with appropriate text (see modif.draw)
			for(i in 3:length(counts)) {
				plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), main=paste(r[1], 'versus', r[2], 'intersect from', r[i], '\n(pv=', signif(pvs[i], 3), ')'), xlab='', xaxt='n', ylab='', yaxt='n', bty='n')
				grid.draw(draws[[i]])
			}
		dev.off()
	} else {
		#generate combination
		if(length(peaks) == 3) {
			label <- c('n12', 'n23', 'n13', 'n123')
		} else if(length(peaks) == 4) {
			label <- paste('n', c(12:14, 23:24, 34, 123:124, 134, 234, 1234), sep='')
		} else if(length(peaks) == 5) {
			label <- paste('n', c(12:15, 23:25, 34:35, 45, 123:125, 134:135, 145, 234:235, 245, 345, 1234:1235, 1245, 1345, 2345, 12345), sep='')
		}
		inter <- rep(NA, length(label)); names(inter) <- label;
		#run all compar_bed
		compar <- peaks
		i <- 1
		while(sum(is.na(inter)) > 0) {
			j <- i+1
			while(j <= length(compar)) {
				current <- compar[c(i, j)]
				nam <- paste(names(current)[1], names(current)[2], sep='-vs-')
				label <- get.label(nam, names(peaks))
				if(is.na(inter[label])) {
					print(paste(label, nam, sep=': '))
					prefix <-  file.path(outdir, nam)
					#run compar_bed
					cmd=paste('./compare_bed.py', '-f', current[1], '-F', current[2], '-o', prefix, '-n', names(current)[1], '-N', names(current)[2])
					print(cmd)
					system(cmd)
					#counts
					counts <- extract.counts(prefix=prefix, r1=names(current)[1], r2=names(current)[2])
					#select best
					idx <- get.best(counts, method=sel)
					best <- names(counts)[idx]; names(best) <- nam;
					inter[label] <- counts[idx]				
					#Avoid confict
					inter[label] <- check.conflict(inter, label)
					#Add to queue
					compar <- c(compar, best)
				}
				j <- j + 1
			}
			i <- i + 1
		}
		prefix <- paste(names(peaks), collapse='-')		
		compar <- cbind(names(compar), compar)
		write.table(compar, file.path(outdir, paste(prefix, '_selectedFiles.txt', sep='')), quote=FALSE, row.names=FALSE, col.names=FALSE)
		#draw venn
		areas <- rep(NA, length(peaks))
		for(i in 1:length(peaks)) { areas[i] <- length(readLines(con=peaks[i])) }
		print(areas)
		print(inter)
		pdf(file.path(outdir, paste(prefix, 'pdf', sep='.')))
			plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), main=prefix, xlab='', xaxt='n', ylab='', yaxt='n', bty='n')
			if(length(peaks) == 3) {
				draw.triple.venn(area1=areas[1], area2=areas[2], area3=areas[3], n12=inter[1], n23=inter[2], n13=inter[3], n123=inter[4],
									category=names(peaks), col=c('royalblue4', 'springgreen4', 'red3'), lwd=2)
			} else if(length(peaks) == 4) {
				draw.quad.venn(area1=areas[1], area2=areas[2], area3=areas[3], area4=areas[4], n12=inter[1], n13=inter[2], n14=inter[3], n23=inter[4], n24=inter[5], n34=inter[6],
								n123=inter[7], n124=inter[8], n134=inter[9], n234=inter[10], n1234=inter[11], category=names(peaks), col=c('royalblue4', 'springgreen4', 'red3', 'orange'), lwd=2)
			} else if(length(peaks) == 5) {
				draw.quintuple.venn(area1=areas[1], area2=areas[2], area3=areas[3], area4=areas[4], area5=areas[5], n12=inter[1], n13=inter[2], n14=inter[3], n15=inter[4],
									n23=inter[5], n24=inter[6], n25=inter[7], n34=inter[8], n35=inter[9], n45=inter[10], n123=inter[11], n124=inter[12], n125=inter[13],
									n134=inter[14], n135=inter[15], n145=inter[16], n234=inter[17], n235=inter[18], n245=inter[19], n345=inter[20],
									n1234=inter[21], n1235=inter[22], n1245=inter[23], n1345=inter[24], n2345=inter[25], n12345=inter[26],
									category=names(peaks), col=c('royalblue4', 'springgreen4', 'red3', 'orange', 'purple'), lwd=2)
			}
		dev.off()
	}
}
