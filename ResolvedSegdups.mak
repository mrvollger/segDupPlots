TITLE=asm
all: $(ASM).sam \
   $(ASM).fai \
   $(ASM).bed5 \
   $(ASM).resolved \
   $(TITLE)-SegdupResolutionByPoint.pdf \
   $(ASM).euchromatic.resolved \
   $(ASM).euchromatic.unresolved \
   $(TITLE)-Mean-SegdupResolutionByDensity.pdf \
   $(ASM).mean.resolved \
   $(ASM).bb

help:
	@echo "Usage: make -f CompareAssembly.mak ASM=assembly.fa "	
	@echo "  [REF=/path/to/ref]  The files REF.sa, and REF.contig_sizes (2nd col of .fai file) must exist."
	@echo "  [SEGDUPS=file.segdups]  The segdups file. Must have 'chrom start end identity' in the format"
	@echo "  [TITLE= name of plot title]  "

HOME=~mchaisso
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/indexes/GRCh38.fasta
MAX_SEGDUPS=/net/eichler/vol5/home/mchaisso/projects/SegDups/GRCh38/grch38_superdups.max_identity.bed
MEAN_SEGDUPS=/net/eichler/vol5/home/mchaisso/projects/SegDups/GRCh38/grch38_superdups.mean_identity.bed
EUCHROMATIC?=/net/eichler/vol5/home/mchaisso/projects/GenomeTracks/GRCh38/euchromatic.hg38.bed
BLASR=/net/eichler/vol21/projects/bac_assembly/nobackups/scripts/testResolvedSegDups/blasr/alignment/bin/blasr
TITLE?=$(ASM)

$(ASM).fai: $(ASM)
	samtools faidx $(ASM)

$(ASM).sam: $(ASM)
	blasr $(ASM) $(REF) -alignContigs -nproc 8 -minMapQV 30 -out $@ -sam

$(ASM).bed5: $(ASM).sam
	cat $< | grep -v "^@" | awk '{ if ($$5 > 20) print $$3"\t"$$4"\t"$$4+$$9"\t"$$1"\t"$$5 }'> $@

$(ASM).bb: $(ASM).bed5
	bedtools sort -i $(ASM).bed5 > $(ASM).bed5.sorted
	module load ucsc && bedToBigBed $(ASM).bed5.sorted $(REF).fai $(ASM).bb -type=bed5

$(ASM).resolved: $(ASM).bed5
	bedtools intersect -a $^ -b $(MAX_SEGDUPS) -wa -wb > $(ASM).segdup_intersect
	$(HOME)/projects/AssemblyAnalysis/scripts/PrintResolvedSegdups.py $(ASM).segdup_intersect $(MAX_SEGDUPS) --resolved $@ --unresolved $(ASM).unresolved  --extra 50000

$(ASM).mean.resolved: $(ASM).bed5
	bedtools intersect -a $^ -b $(MEAN_SEGDUPS) -wa -wb > $(ASM).segdup_intersect
	$(HOME)/projects/AssemblyAnalysis/scripts/PrintResolvedSegdups.py $(ASM).segdup_intersect $(MEAN_SEGDUPS) --resolved $@ --unresolved $(ASM).mean.unresolved  --extra 100000


$(TITLE)-SegdupResolutionByPoint.pdf: $(ASM).resolved
	Rscript $(HOME)/projects/AssemblyAnalysis/scripts/PlotResolved.R $(ASM).resolved $(ASM).unresolved  $(TITLE) $(TITLE)
	convert -density 150  $(TITLE)-SegdupResolutionByDensity.pdf $(TITLE)-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-SegdupResolutionByPoint.pdf $(TITLE)-SegdupResolutionByPoint.png


$(TITLE)-Mean-SegdupResolutionByDensity.pdf: $(ASM).mean.resolved
	Rscript $(HOME)/projects/AssemblyAnalysis/scripts/PlotResolved.R $(ASM).mean.resolved $(ASM).mean.unresolved  $(TITLE)-Mean $(TITLE)-Mean
	convert -density 150  $(TITLE)-Mean-SegdupResolutionByDensity.pdf $(TITLE)-Mean-SegdupResolutionByDensity.png
	convert -density 150  $(TITLE)-Mean-SegdupResolutionByPoint.pdf $(TITLE)-Mean-SegdupResolutionByPoint.png


$(ASM).euchromatic.resolved: $(ASM).resolved
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

$(ASM).euchromatic.unresolved: $(ASM).unresolved
	bedtools intersect -a $< -b $(EUCHROMATIC) > $@

