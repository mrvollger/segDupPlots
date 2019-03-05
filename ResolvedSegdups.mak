TITLE=asm
ADJPLOTS=$(TITLE).adjustingExtra
all: $(TITLE).mash \
	$(ASM).fai \
	$(TITLE).bed5 \
	$(TITLE).bed \
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
	$(TITLE).inter.dups \
	$(ADJPLOTS)Col/zeroTo250k.png

make_file := $(abspath $(lastword $(MAKEFILE_LIST)))
BASE_DIR := $(shell dirname $(make_file))

help:
	@echo "Make file location $(BASE_DIR)"
	@echo "Usage: make -f $(BASE_DIR)/ResolvedSegdups.mak ASM=assembly.fa TITLE=your.prefix"	
	@echo "     REF has to be updated to a no alts copy of ucsc's hg38"


# THESE NEED TO BE CHANGED BY THE USER
REF=/data/korens/devel/segDupPlots/ucsc.hg38.noalts.fasta

EUCHROMATIC=$(BASE_DIR)/ucsc.hg38.no_alts.fasta.euchromatic.bed
GENES=$(BASE_DIR)/ucsc.hg38.gene.locations.bed 
MAX_SEGDUPS=$(BASE_DIR)/ucsc.merged.max.segdups.bed
MEAN_SEGDUPS=$(BASE_DIR)/ucsc.merged.mean.segdups.bed
COL_SEGDUPS=$(BASE_DIR)/ucsc.collapsed.sorted.segdups

TITLE?=$(ASM)

$(ASM).fai: $(ASM)
	samtools faidx $(ASM)

$(TITLE).mash: $(ASM)
	mashmap -r $(REF) -q $(ASM) -t 8 -o $(TITLE).mash

$(TITLE).bed5: $(TITLE).mash
	# the 250 number is jsut a fake mapq
	awk '{ if ($$10 > 85) print  $$6"\t"$$8"\t"$$9"\t"$$1"\t"250  }' ${TITLE}.mash > $@
	#cat $< | grep -v "^@" | awk '{ if ($$5 > 20) print $$3"\t"$$4"\t"$$4+$$9"\t"$$1"\t"$$5 }'> $@

$(TITLE).bed: $(TITLE).mash
	# the 250 number is jsut a fake mapq
	awk '{ if ($$10 > 85) print  $$6"\t"$$8"\t"$$9"\t"$$1"\t"$$3"\t"$$4"\t"$$2"\t"$$10  }' ${TITLE}.mash | \
	   bedtools sort -i - >	$@
	#cat $< | grep -v "^@" | awk '{ if ($$5 > 20) print $$3"\t"$$4"\t"$$4+$$9"\t"$$1"\t"$$5 }'> $@



$(TITLE).bb: $(TITLE).bed5
	bedtools sort -i $(TITLE).bed5 > $(TITLE).bed5.sorted
	module load ucsc && bedToBigBed $(TITLE).bed5.sorted $(REF).fai $(TITLE).bb -type=bed5




$(TITLE).resolved: $(TITLE).bed5
	bedtools intersect -a $^ -b $(MAX_SEGDUPS) -wa -wb  \
		> $(TITLE).segdup_intersect
	$(BASE_DIR)/PrintResolvedSegdups.py $(TITLE).segdup_intersect \
		$(MAX_SEGDUPS) --resolved $(TITLE).resolved --unresolved $(TITLE).unresolved \
		--allseg $(TITLE).all --extra 50000

$(TITLE).mean.resolved: $(TITLE).bed5
	bedtools intersect -a $^ -b $(MEAN_SEGDUPS) -wa -wb \
		> $(TITLE).segdup_intersect.mean
	$(BASE_DIR)/PrintResolvedSegdups.py $(TITLE).segdup_intersect.mean \
		$(MEAN_SEGDUPS) --resolved $(TITLE).mean.resolved --unresolved $(TITLE).mean.unresolved \
		--allseg $(TITLE).mean.all --extra 50000

$(TITLE).col.resolved $(TITLE).col.all: $(TITLE).bed5
	bedtools intersect -a $^ -b $(COL_SEGDUPS) -wa -wb \
		> $(TITLE).segdup_intersect.col
	$(BASE_DIR)/PrintResolvedSegdups.py $(TITLE).segdup_intersect.col \
	   	$(MEAN_SEGDUPS) --resolved $(TITLE).col.resolved --unresolved $(TITLE).col.unresolved \
		--allseg $(TITLE).col.all --extra 50000





$(TITLE)-SegdupResolutionByPoint.pdf: $(TITLE).resolved
	Rscript $(BASE_DIR)/PlotResolved.R $(TITLE).resolved $(TITLE).unresolved  $(TITLE) $(TITLE)
	convert -density 150  $(TITLE)-SegdupResolutionByDensity.pdf $(TITLE)-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-SegdupResolutionByPoint.pdf $(TITLE)-SegdupResolutionByPoint.png


$(TITLE)-Mean-SegdupResolutionByDensity.pdf: $(TITLE).mean.resolved
	Rscript $(BASE_DIR)/PlotResolved.R $(TITLE).mean.resolved $(TITLE).mean.unresolved  $(TITLE)-Mean $(TITLE)-Mean
	convert -density 150  $(TITLE)-Mean-SegdupResolutionByDensity.pdf $(TITLE)-Mean-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-Mean-SegdupResolutionByPoint.pdf $(TITLE)-Mean-SegdupResolutionByPoint.png


$(TITLE)-Col-SegdupResolutionByDensity.pdf: $(TITLE).col.resolved
	Rscript $(BASE_DIR)/PlotResolved.R $(TITLE).col.resolved $(TITLE).col.unresolved  $(TITLE)-Col $(TITLE)-Col
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
	Rscript $(BASE_DIR)/SegDupPlots.R --dest $(ADJPLOTS) --all $(TITLE).all --euchro $(TITLE).euchromatic.all --gene $(TITLE).unresolved.gene.intersection


$(ADJPLOTS)Mean/zeroTo250k.png: $(TITLE).mean.all $(TITLE).euchromatic.mean.all $(TITLE).unresolved.gene.intersection
	mkdir -p $(ADJPLOTS)Mean 
	Rscript $(BASE_DIR)/SegDupPlots.R --dest $(ADJPLOTS)Mean --all $(TITLE).mean.all --euchro $(TITLE).euchromatic.mean.all --gene $(TITLE).unresolved.gene.intersection


$(ADJPLOTS)Col/zeroTo250k.png: $(TITLE).col.all $(TITLE).euchromatic.col.all $(TITLE).unresolved.gene.intersection
	mkdir -p $(ADJPLOTS)Col
	Rscript $(BASE_DIR)/SegDupPlots.R --dest $(ADJPLOTS)Col --all $(TITLE).col.all --euchro $(TITLE).euchromatic.col.all --gene $(TITLE).unresolved.gene.intersection


$(TITLE).unresolved.gene.intersection: $(TITLE).resolved
	bedtools intersect -a $(GENES) -b $(TITLE).unresolved -wb > $(TITLE).unresolved.gene.intersection


$(TITLE).inter.dups: $(TITLE).bed5	
	bedtools intersect -a $(MAX_SEGDUPS) -b $(TITLE).bed5 > $(TITLE).inter.dups


