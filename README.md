DER Finder Beta
=========

This is the beta version of the R package for DER Finder, a method for analyzing differential expression with RNA-seq data.  DER Finder is described in detail in the paper [Differential expression analysis of RNA-seq data at single-base resolution](http://biostatistics.oxfordjournals.org/content/early/2014/01/06/biostatistics.kxt053.full). Code for the analysis used in the paper is in `analysis_code.R`. 

The package in this repository is an exact implementation of the methods described in the DER Finder manuscript. The final version of the software is significantly more efficient than this version and is [available on Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/derfinder.html). For research projects, **please use the Bioconductor package** instead of the package in this repository, which mostly exists for historical/reproducibility reasons. The official package is written and maintained by Leonardo Collado Torres, and issues can be submitted [to its GitHub repository](https://github.com/lcolladotor/derfinder). 

# installation
You can install this package directly from GitHub using the `install_github` function from the `devtools` package:
```S
library(devtools)
install_github('derfinder', 'alyssafrazee')
```
You may also need to install `Genominator`, `limma`, `HiddenMarkov`, `splines`, and `locfdr`.

# preprocessing data
Before using the `derfinder` R package, raw RNA-seq data should be processed as follows:

### read alignment
We use [TopHat](http://tophat.cbcb.umd.edu/) to align reads, but any junction-aware aligner that writes alignments in `.bam` format is appropriate. 

After aligning reads, `.bam` files need to be indexed with `samtools index`. (See the [SAMtools page](http://samtools.sourceforge.net/) for information about using and installing SAMtools.)

### coverage calculation
DER Finder operates at single-base resolution, so we need per-nucleotide read coverage for every genomic position. We provide a Python script, `countReads.py`, that does this on a per-chromosome basis. For each sample, and for each chromosome you would like per-base counts for, you need to run:
```
python countReads.py --file accepted_hits.bam --output <OUTFILE> --kmer <READ_LENGTH> --chrom <CHROMOSOME> --stranded <IS_STRANDED>
```
where `accepted_hits.bam` is the file of read alignments for the sample (this is the default output name from TopHat, but if you used a different aligner, this filename might be different), `<OUTFILE>` is the name of the fill that will contain the per-nt counts, `<READ_LENGTH>` is the length of the RNA-seq reads in your sample, and `<CHROMOSOME>` is the chromosome you want to count. The chromosome name should match the chromosome names in `accepted_hits.bam`. `<IS_STRANDED>` should be one of `TRUE`, `REVERSE`, or `FALSE`. If you used a stranded protocol, `<IS_STRANDED>` should be `TRUE`. If you used a reversely-stranded protocol (e.g., an Illumina protocol where strand labels are reversed), `<IS_STRANDED>` should be `REVERSE`. If your protocol wasn't stranded, you can eliminate this option (it's `FALSE` by default).

`countReads.py` depends on the [pysam](https://code.google.com/p/pysam/) module, which requires a Python version >2.4 and <3.0.

`<OUTFILE>` will be a 2-column tab-delimited text file, with genomic position in the first column and coverage in the second column. If `--stranded` is `TRUE`, two output files, called `OUTFILE_plus` and `OUTFILE_minus`, will be produced (one for each strand). If you used a stranded protocol and would like to find separate DERs for each strand, you will need to run the pipeline twice, starting here (once on each of the resulting coverage matrices).

### merging coverage files
The `<OUTFILES>` for each sample will need to be merged to create a nucleotide-by-sample matrix for analysis. If you have a giant computer (i.e., if you can hold this matrix in memory), you can do this in with some code like this: 
```S
# which chromosome?
chr = "chr22"

# sample IDs
samps = c("orbFrontalF1", "orbFrontalF2", "orbFrontalF3", "orbFrontalF11", "orbFrontalF23", "orbFrontalF32", "orbFrontalF33", "orbFrontalF40", "orbFrontalF42", "orbFrontalF43", "orbFrontalF47", "orbFrontalF53", "orbFrontalF55", "orbFrontalF56", "orbFrontalF58")

# read in each sample:
countlist = list()
for(s in 1:length(samps)){
    print(paste("reading: sample",samps[s]))
    y = read.table(paste0(samps[s], "-outfile.txt"), sep="\t", header=FALSE)
    print("done reading.")
    countlist[[s]] = y$V2
    if(s==1) pos = y$V1
    if(s>1){
        if(length(y$V1)>length(pos)) pos = y$V1
    }
    rm(y);gc();gc();gc();gc()
  }

# put samples together and zero-pad as needed:
thelen = length(pos)
for(i in 1:length(countlist)){
    countlist[[i]] = c(countlist[[i]],rep(0,thelen-length(countlist[[i]])))
}
names(countlist) = samps
chr.table = as.data.frame(countlist)
chr.table = data.frame(pos, chr.table)
write.table(chr.table, file=paste0(chr, '_allcounts.txt'), row.names=FALSE, quote=FALSE, sep="\t")
```
This isn't super efficient, but it gets the job done. You can also use something like the Unix `paste` utility, making sure to account for the fact that each sample may start/end at slightly different positions.

# loading data
After aligning and counting reads and merging each sample's read counts, we can begin using the `derfinder` package. We don't want to read the huge matrices into R, so we utilize SQLite databases - but everything is built into the R package, so no SQL commands are actually needed. The `makeDb` command creates the database, and you can explore the database with `Genominator`:

```S
# create database from merged file:
library(derfinder)
makeDb(dbfile = 'chr22_allcounts.db', textfile='chr22_allcounts.txt', tablename = 'chr22')

# explore data:
library(Genominator)
dat = ExpData(dbFilename = 'chr22_allcounts.db', tablename = 'chr22')
head(dat)
getColnames(dat)
length(dat[,1]$pos)
```

# nucleotide-level statistical analysis:
As described in the manuscript, at each nucleotide, we see if a covariate of interest is associated with expression (i.e., coverage) using a linear model with coverage as outcome. The code to do this is:
```S
# define covariate:
sex = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1)

# fit initial model (no shrinkage):
limma.input = getLimmaInput(dbfile = 'chr22_allcounts.db', tablename = 'chr22', group = sex)

# get positions where model was fit:
pos = limma.input$pos

# calculate moderated t-statistics and estimated fold changes:
tstats = getTstats(fit = limma.input$ebobject, trend=TRUE)
tt = tstats$tt
logfchange = tstats$logfchange
```

# region-level statistical analysis
We now fit a Hidden Markov Model to the test statistics from each base-pair to get regions of similar differential expression signal:
```S
params = getParams(tt)
regions = getRegions(method='HMM', chromosome='chr22', pos = pos, 
  tstats=tt, stateprobs=params$stateprobs, params=params$params,
  includet=TRUE, includefchange=TRUE, fchange=logfchange)
#### ^ may take some time
head(regions$states)
ders = subset(regions$states, state==3|state==4)
head(ders)
```
We can get region-level p-values using a permutation test (details in the manuscript). P-values can be adjusted for multiple testing (for false discovery rate control) by being converted to q-values.
```S
pvals = get.pvals(regions=regions$states, dbfile='chr22_allcounts.db',
  tablename='chr22', num.perms=1000, group=sex, est.params=params,
  chromosome='chr22') #takes a while
regions.withp = data.frame(regions$states, pvals=pvals, qvals=p.adjust(pvals, 'fdr'))
head(reg.withp)
```

# connecting results to annotation
The `getAnnotation` function allows users to download exon annotation information from UCSC, though changes to this database often break this function, so a more stable version will be available in the release version of the package. You can match regions to annotated exons using `getExons`:
```S
exons = getAnnotation('hg19', 'knownGene')
chr22exons = subset(exons, chr=='chr22')
getExons(region=c('chr22', 60234591, 60234791), annotation=chr22exons)
```

# visualizing results
You can also plot specific regions, exons, or genes using one of the available plotting functions.  For example, to plot the second region in your `regions` object:
```S
plotRegion(regions, ind=2, tstats=tt, pos=pos, annotation=chr22exons,
  counts='chr22_allcounts.db', group=ifelse(sex==1, 'male', 'female'),
  tabname='chr22', chromosome='chr22')
```

# support
This version of _derfinder_ is no longer actively maintained, since the official Bioconductor package is much more efficient. Please see [the derfinder GitHub repo](https://github.com/lcolladotor/derfinder) for support on that version.

