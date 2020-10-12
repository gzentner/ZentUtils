## ZentUtils

---

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
First, a 6-column BED file is imported into ZentUtils and converted to a GRanges object that is
stored in the ```regions``` slot of a new ZentUtils object.

```R
zent <- zentutils(system.file("extdata", "homer_reb1_badis_motif_scan.bed", package = "ZentUtils"), genome = "sacCer3")
```

##### Sort regions
Next, we sort the regions. Regions can be sorted by coordinate with ```by = "coord"``` or
by score (5th column of the BED file) with ```by = "score"```.

```R
zent <- sort_regions(zent, by = "score", decreasing = TRUE)

```

##### Expand regions
We then optionally expand the regions of interest to a specified size. For expansion,
```length/2``` bases are added to both ends of each region. Expanded regions are stored
in the ```expanded_regions``` slot of the ZentUtils object. Sequences that are out-of-bounds
following expansion are discarded.

```R
zent <- expand_regions(zent, length = 10)
```

##### Retrieve region sequences
Sequences for the imported or expanded regions are then obtained from a BSgenome object
and stored in the ```seqs``` slot of the ZentUtils object.

```R
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
zent <- get_seqs(zent, region_type = "expanded", genome = genome)
```

##### Generate sequence color map
Lastly, we generate a tile plot, wherein each base is represented by a specific
color. This allows efficient visual detection of specific sequences patterns
around features of interest.

```R
color_map(zent)
```

![Reb1 motif match color map](https://github.com/gzentner/ZentUtils/blob/master/inst/images/reb1_motif_colormap.png)
