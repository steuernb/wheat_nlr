# The bread wheat intracellular immune receptor repertoire


This repository is a collection of source code that supported the analyses of the NLR gene family in hexaploid wheat.

The code in this repository is to support reproducibility. It is not designed as stand-allone applications. 

The **NLR-Annotator**, which is also described in this paper can be found [here](https://github.com/steuernb/NLR-Annotator).


## Content

### `java/src/CalculateNLRGenes.java`
Methods to generate supplementary tables S5 and S6

Two files `data/ADR1_addon.txt` and `data/manual_inspection_nlrloci.txt` are input data for this code. These files cannot be gererated based on public data alone and are added here.

-----------

### `java/src/PreparePhylogenetics.java`
Methods to generate annotation files for [iTOL](https://itol.embl.de/shared/steuernb). Basis for Figures 2, S4, S5 and S6.

-----------
### `NLR_Expression_Level.java`
Methods to create Figure 3. 


-----------

### `java/src/PAMP_RnaSeq.java`
Prepare data for ternary plots (Figure 4).


### `R/expression_triangles.R`
Plot ternary plots for NLR expression



-----------

### support

`BioSequence.java`, `BioSequenceReader.java` and `FastaReader.java` are support classes that are used by methods in source code above.
