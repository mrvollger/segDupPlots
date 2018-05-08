TITLE=asm
ADJPLOTS=adjustingExtra
all: $(TITLE).mash \
   $(ASM).fai \
   $(TITLE).bed5 \
   $(TITLE).resolved \
   $(TITLE)-SegdupResolutionByPoint.pdf \
   $(TITLE).euchromatic.resolved \
   $(TITLE).euchromatic.unresolved \
   $(TITLE).euchromatic.all \
   $(TITLE)-Mean-SegdupResolutionByDensity.pdf \
   $(TITLE).mean.resolved \
   $(TITLE).bb \
   $(ADJPLOTS)Mean/zeroTo250k.png \
   $(ADJPLOTS)/zeroTo250k.png \
   $(TITLE).unresolved.gene.intersection \

# make -f ~/projects/segDupPlots/ResolvedSegdups.mak ASM=../Yoruban.fasta TITLE=asm
help:
	@echo "Usage: make -f CompareAssembly.mak ASM=assembly.fa "	
	@echo "  [REF=/path/to/ref]  The files REF.sa, and REF.contig_sizes (2nd col of .fai file) must exist."
	@echo "  [SEGDUPS=file.segdups]  The segdups file. Must have 'chrom start end identity' in the format"
	@echo "  [TITLE= name of plot title]  "

HOME=~mchaisso

#REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/indexes/GRCh38.fasta
REF=~mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta

#MAX_SEGDUPS=/net/eichler/vol5/home/mchaisso/projects/SegDups/GRCh38/grch38_superdups.max_identity.bed
MAX_SEGDUPS=~mvollger/assemblies/hg38/ucsc.merged.max.segdups.bed

#MEAN_SEGDUPS=/net/eichler/vol5/home/mchaisso/projects/SegDups/GRCh38/grch38_superdups.mean_identity.bed
MEAN_SEGDUPS=~mvollger/assemblies/hg38/ucsc.merged.mean.segdups.bed

# collapsed segdups
COL_SEGDUPS=/net/eichler/vol2/home/mvollger/assemblies/hg38/ucsc.collapsed.sorted.segdups

#EUCHROMATIC?=/net/eichler/vol5/home/mchaisso/projects/GenomeTracks/GRCh38/euchromatic.hg38.bed
EUCHROMATIC?=~mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta.euchromatic.bed

BLASR=/net/eichler/vol21/projects/bac_assembly/nobackups/scripts/testResolvedSegDups/blasr/alignment/bin/blasr
MASHMAP=/net/eichler/vol5/home/mchaisso/software/MashMap/mashmap
TITLE?=$(ASM)

$(ASM).fai: $(ASM)
	samtools faidx $(ASM)

$(TITLE).mash: $(ASM)
	#blasr $(ASM) $(REF) -sa $(REF).sa -alignContigs -nproc 8 -minMapQV 30 -out $@ -sam
	$(MASHMAP) -r $(REF) -q $(ASM) -t 8
	mv mashmap.out $(TITLE).mash	
	# output is nammed mashmap.out
	#--perc_identity 95 -t 8 

$(TITLE).bed5: $(TITLE).mash
	# the 250 number is jsut a fake mapq
	awk '{ if ($$10 > 85) print  $$6"\t"$$8"\t"$$9"\t"$$1"\t"250  }' ${TITLE}.mash > $@
	#cat $< | grep -v "^@" | awk '{ if ($$5 > 20) print $$3"\t"$$4"\t"$$4+$$9"\t"$$1"\t"$$5 }'> $@

$(TITLE).bb: $(TITLE).bed5
	bedtools sort -i $(TITLE).bed5 > $(TITLE).bed5.sorted
	module load ucsc && bedToBigBed $(TITLE).bed5.sorted $(REF).fai $(TITLE).bb -type=bed5

$(TITLE).resolved: $(TITLE).bed5
		#awk -F"\t" '{{print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9*100}}' 
	bedtools intersect -a $^ -b $(MAX_SEGDUPS) -wa -wb  \
		> $(TITLE).segdup_intersect
	~mvollger/projects/segDupPlots/PrintResolvedSegdups.py $(TITLE).segdup_intersect $(MAX_SEGDUPS) --resolved $@ --unresolved $(TITLE).unresolved --allseg $(TITLE).all --extra 50000

$(TITLE).mean.resolved: $(TITLE).bed5
		#awk -F"\t" '{{print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9*100}}' 
	bedtools intersect -a $^ -b $(MEAN_SEGDUPS) -wa -wb \
		> $(TITLE).segdup_intersect.mean
	~mvollger/projects/segDupPlots/PrintResolvedSegdups.py $(TITLE).segdup_intersect.mean $(MEAN_SEGDUPS) --resolved $@ --unresolved $(TITLE).mean.unresolved --allseg $(TITLE).mean.all --extra 50000


$(TITLE).col.resolved $(TITLE).col.all: $(TITLE).bed5
	bedtools intersect -a $^ -b $(COL_SEGDUPS) -wa -wb \
		> $(TITLE).segdup_intersect.col
	~mvollger/projects/segDupPlots/PrintResolvedSegdups.py $(TITLE).segdup_intersect.col \
	   	$(MEAN_SEGDUPS) --resolved $@ --unresolved $(TITLE).col.unresolved \
		--allseg $(TITLE).col.all --extra 50000


$(TITLE)-SegdupResolutionByPoint.pdf: $(TITLE).resolved
	Rscript $(HOME)/projects/AssemblyAnalysis/scripts/PlotResolved.R $(TITLE).resolved $(TITLE).unresolved  $(TITLE) $(TITLE)
	convert -density 150  $(TITLE)-SegdupResolutionByDensity.pdf $(TITLE)-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-SegdupResolutionByPoint.pdf $(TITLE)-SegdupResolutionByPoint.png


$(TITLE)-Mean-SegdupResolutionByDensity.pdf: $(TITLE).mean.resolved
	Rscript $(HOME)/projects/AssemblyAnalysis/scripts/PlotResolved.R $(TITLE).mean.resolved $(TITLE).mean.unresolved  $(TITLE)-Mean $(TITLE)-Mean
	convert -density 150  $(TITLE)-Mean-SegdupResolutionByDensity.pdf $(TITLE)-Mean-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-Mean-SegdupResolutionByPoint.pdf $(TITLE)-Mean-SegdupResolutionByPoint.png


$(TITLE)-Col-SegdupResolutionByDensity.pdf: $(TITLE).col.resolved
	Rscript $(HOME)/projects/AssemblyAnalysis/scripts/PlotResolved.R $(TITLE).col.resolved $(TITLE).col.unresolved  $(TITLE)-Col $(TITLE)-Col
	convert -density 150  $(TITLE)-Col-SegdupResolutionByDensity.pdf $(TITLE)-Col-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-Col-SegdupResolutionByPoint.pdf $(TITLE)-Col-SegdupResolutionByPoint.png


$(TITLE).euchromatic.resolved: $(TITLE).resolved
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

$(TITLE).euchromatic.unresolved: $(TITLE).unresolved
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

$(TITLE).euchromatic.all: $(TITLE).all
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

$(TITLE).euchromatic.mean.all: $(TITLE).mean.all
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

$(TITLE).euchromatic.col.all: $(TITLE).col.all
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

#
# some additional plots that examine how sensative everything is to the --extra parameter
#
$(ADJPLOTS)/zeroTo250k.png: $(TITLE).all $(TITLE).euchromatic.all
	mkdir -p $(ADJPLOTS) 
	Rscript ~mvollger/projects/segDupPlots/SegDupPlots.R --dest $(ADJPLOTS) --all $(TITLE).all --euchro $(TITLE).euchromatic.all

$(ADJPLOTS)Mean/zeroTo250k.png: $(TITLE).mean.all $(TITLE).euchromatic.mean.all
	mkdir -p $(ADJPLOTS)Mean 
	Rscript ~mvollger/projects/segDupPlots/SegDupPlots.R --dest $(ADJPLOTS)Mean --all $(TITLE).mean.all --euchro $(TITLE).euchromatic.mean.all


$(ADJPLOTS)Col/zeroTo250k.png: $(TITLE).col.all $(TITLE).euchromatic.col.all
	mkdir -p $(ADJPLOTS)Col
	Rscript ~mvollger/projects/segDupPlots/SegDupPlots.R --dest $(ADJPLOTS)Col --all $(TITLE).col.all --euchro $(TITLE).euchromatic.col.all




$(TITLE).unresolved.gene.intersection: $(TITLE).resolved:
	bedtools intersect -a ~mvollger/assemblies/hg38/hg38.gene.locations.bed -b $(TITLE).unresolved -wb > $(TITLE).unresolved.gene.intersection










