# segDupPlots

## Install
Mashmap must be installed and in your path. The best way to do this is through conda. 
```
conda install -c bioconda mashmap 
```
One line in the makefile ResolvedSegdups.mak must be changed as well. `REF=/net/...` must be updated
to point at your local download of ucsc's hg38 with no alts. 

## Run 
```
make -f {path.to.local.repo}/ResolvedSegdups.mak ASM={path.to.de.novo.assembly} TITLE={prefix.to.add.to.output}
```
Output plots can then be found in `adjustingExtra/` direcotry.


## Methods 
One of the output plots shows the percentage of segmental duplications that are "Resolved." Our definition of resolved is that for a SD to be resolved the assembly must continue into unique sequence on either side of the SD by at least some minimal extension. The percentage plot shows the fraction of resolved bases as the minimal extension is varied from 0 to 250 kbp.

The basic steps of identifying resolved vs unresolved duplications is as follows:
1. Map the de novo assembly to the human reference using Mashmap 2.0 defaults. 
2. Download the UCSCs annotated SD track and merge overlapping segmental duplications by their maximum percent identity. 
3. Intersect the de novo assembly track with the modified segmental duplication track.
4. Determine if and by how much the de novo assembly extends past segmental duplication blocks on either side.
5. Mark segmental duplications as resolved or unresolved based on whether the de novo assembly extends at least X kb into unique sequence on either side. 
6. Plot the percentage of segmental duplication bases resolved as a function of the minimal extension into unique sequence past a duplication block.


