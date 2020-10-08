## ZentUtils

<<<<<<< HEAD
=======
---

>>>>>>> c4987e96dd5894730bcfaca0c86df17bb50b4227
#### About
ZentUtils is a collection of functions designed to streamline analysis of genomic
intervals.

#### What does it do?
* Region sorting by coordinate or score
* Region expansion
* Retrieval of region sequences
* Plotting of region sequences as a color map

---

#### Example: plotting sequences around motif matches
Included with ZentUtils is a BED file of matches to motifs for the yeast Reb1
transcription factor. 

##### Import regions into ZentUtils

<<<<<<< HEAD
```R
=======
```
>>>>>>> c4987e96dd5894730bcfaca0c86df17bb50b4227
zent <- zentutils(system.file("extdata", "homer_reb1_badis_motif_scan.bed", package = "ZentUtils"), genome = "sacCer3")
```

##### Sort regions

<<<<<<< HEAD
```R
=======
```
>>>>>>> c4987e96dd5894730bcfaca0c86df17bb50b4227
zent <- sort_regions(zent, by = "score", decreasing = TRUE)
```

##### Expand regions

<<<<<<< HEAD
```R
=======
```
>>>>>>> c4987e96dd5894730bcfaca0c86df17bb50b4227
zent <- expand_regions(zent, length = 50)
```
