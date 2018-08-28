# segDupPlots

# Run 
```
make -f ResolvedSegdups.mak ASM={path.to.de.novo.assembly} TITLE={prefix.to.add.to.output}
```

# Methods 
The percentage plot reflects the percentage of bases within segmental duplications that are “Resolved.” Our definition of resolved is nuanced but the basic idea is that for a segdup to be resolved the assembly must continue into unique sequence on either side of the segmental duplication by at least some minimal extension. The basic steps are as follows:
1. Map the de novo assembly to the human reference using Mashmap 2.0 defaults. 
2. Download the UCSCs annotated segdup track and merge overlapping segmental duplications by their maximum percent identity. 
3. Intersect the de novo assembly track with the modified segmental duplication track
4. Determine if and by how much the de novo assembly extends past segmental duplication blocks on either side.
5. Mark segmental duplications as resolved or unresolved based on whether the de novo assembly extends at least X kb into unique sequence on either side. (note a 50 kbp threshold was chosen for all the red and black dot plots showing resolved and unresolved duplications)
6. Plot the percentage of segmental duplication bases resolved as a function of the minimal extension into unique sequence past a duplication block.
7. The percentages reported in the text now represent the range of the percentage of resolved bases obtained by varying the minimal extension from 0-250 kbp. 

