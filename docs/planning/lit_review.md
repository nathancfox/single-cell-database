# Literature Review

This covers everything I could find that has already been done to create an scRNA-seq database.

## General Patterns

The work mostly falls into one of three types:

1. **Metadata Databases**: Databases that only contain metadata and have done no work to organize,
   parse, or process the actual data itself.
2. **Reference Databases**: Databases that are used as a reference for a separate tool. Here, the
   database is not the main resource. It has been compiled to support a different tool. Unlike
   Metadata Databases, the data has actually parsed and organized.
3. **Full Databases**: Databases where datasets have been parsed, organized, and occasionally
   reprocessed entirely to create a consistent collection of comparable scRNA-seq data.

There are also a number of large datasets that self-host their own data instead of submitting it to
a repository. Other than these large self-contained projects, it looks like all other scRNA-seq
datasets will be found in a data repository like the following:

- GEO (Gene Expression Omnibus)
- ENA (European Nucleotide Archive)
- ArrayExpress
- EGA (European Genome-Phenome Archive)
- CNSA (CNGB (China National GeneBank) Sequence Archive)
- SRA (Sequence Read Archive)

Additionally, the raw data can come in a variety of formats:

- Raw reads (e.g. FASTQ files)
- Alignment Output (e.g. BAM files)
- Expression Matrices
    - UMIs vs NOT UMIs
    - Raw Counts vs RPKM vs FPKM vs TPM
    - Varying gene lists, both depending on the genome used for alignment/quantification and on the
      format of gene ID used (e.g. Ensembl ID vs Human Readable Name)

----

## Existing Work

### Metadata Databases

##### Valentine Svensson's Database

Valentine Svensson, a postdoc in Lior Pachter's lab at Cal Tech maintains probably the most
comprehensive scRNA-seq database in existence. This is purely metadata, he has not processed the
expression data at all. Extremely helpfully, he includes accession numbers across multiple
repositories. There is a lot of metadata, some that I don't want. The entire database is accessible
programmatically as a TSV file.

[https://www.biorxiv.org/content/10.1101/742304v2.full](https://www.biorxiv.org/content/10.1101/742304v2.full)</br>
[http://www.nxn.se/single-cell-studies](http://www.nxn.se/single-cell-studies)

### Reference Databases

##### scRNASeqDB

Web server run by Zhongming Zhao's lab at the University of Texas Health Science Center at Houston.
Source data is not directly available. The main purpose of the tool is to look at a gene's
expression profile across many many cells and datasets. You can technically also look at cell types,
bu these are so specific to dataset, that this isn't really effective. They basically just piped the
metadata from GEO into a SQL database, then used the source data as a backing tool for their web
app.

[https://www.biorxiv.org/content/10.1101/104810v1.full](https://www.biorxiv.org/content/10.1101/104810v1.full)</br>
[https://bioinfo.uth.edu/scrnaseqdb/](https://bioinfo.uth.edu/scrnaseqdb/)

##### Cell Blast

Web server run by Ge Gao's lab at Peking University. It's a tool similar to BLAST, but for finding
cells of similar expression profile to the cells that you submit. It uses 98 reference datasets to
accomplish its task. The datasets were manually curated and pre-processed. There aren't too many
details, just that they were processed for quality. All the datasets are available for download in
HDF5 format via an HTML table along with some basic metadata. I don't like that I don't what
pre-processing they did, but this could be quite helpful. Many of the original studies had
conditions or multiple datasets and these have been separated here.

[https://www.biorxiv.org/content/10.1101/587360v2.full](https://www.biorxiv.org/content/10.1101/587360v2.full)</br>
[https://cblast.gao-lab.org/download](https://cblast.gao-lab.org/download)

#### CytoTRACE

Web server run by Aaron Newman's lab at Stanford. It is a tool for predicting differentiation cell
state. You can upload scRNA-seq data and it will predict differentiation paths for the cells. It
relies on 42 public scRNA-seq datasets. The metadata for these datasets is available for download. I
don't think the actual data is available. Also unclear if they pre-processed the data at all.

[https://www.biorxiv.org/content/10.1101/649848v1.full](https://www.biorxiv.org/content/10.1101/649848v1.full)</br>
[https://cytotrace.stanford.edu/](https://cytotrace.stanford.edu/)

###### Cell Atlas Search

Similar to Cell Blast. It's a cell-to-expression profile matcher. You submit single cell expression
profiles and it returns the most similar single-cell or bulk-tissue matches. The reference data was
compiled from the 10X datasets, GEO, and Recount2. The backing reference data is not available.

[http://103.25.231.35:8097/cellatlassearch/](http://103.25.231.35:8097/cellatlassearch/)

### Full Databases

##### Hemberg Datasets

Martin Hemberg's lab at the Sanger Institute put together 14 human scRNA-seq datasets and 28 mouse
datasets. It looks like they did it pretty carefully and manually. They provide scripts to download
the datasets and do the parsing that they did. The metadata is processed with a Perl script that I
might steal. However, it's not in a database, just a website that you have to access manually.

[https://hemberg-lab.github.io/scRNA.seq.datasets/](https://hemberg-lab.github.io/scRNA.seq.datasets/)

##### SCPortalen

Database of ~70 human and mouse scRNA-seq datasets, run by the RIKEN institute. They've done a fair
amount of processing on the data and include PCA and t-SNE plots. It looks like it's all in FPKM?
The interface is very clunky and not very easily accessed. Links to some metadata, but again very
clunky. However, there is a place to batch download the datasets. Haven't tried it yet.

[https://academic.oup.com/nar/article/46/D1/D781/4555233](https://academic.oup.com/nar/article/46/D1/D781/4555233)</br>
[http://single-cell.clst.riken.jp/non_riken_data/study_meta_info_list.php](http://single-cell.clst.riken.jp/non_riken_data/study_meta_info_list.php)

##### EBI Single Cell Expression Atlas

Database of 132 single cell expression studies. Direct download links and some minor metadata. I
checked one dataset and it also came with metadata internally. Looks like the original data has a
repository ID (e.g. E-MTAB-7195). Last updated in October, not bad.

[https://www.ebi.ac.uk/gxa/sc/experiments](https://www.ebi.ac.uk/gxa/sc/experiments)

##### JingleBells

Database, intended for immunology scRNA-seq datasets. They don't provide expression matrices, just
universally reformatted BAM files. Useful, but potentially a pain. Could be useful if my QC
standards include re-aligning/re-quantifying. Lists 120 immune datasets and 182 non-immune, HOWEVER,
only 55 immune datasets and 7 non-immune datasets are processed and available. The rest only have
basic metadata.

[https://www.jimmunol.org/content/198/9/3375.long](https://www.jimmunol.org/content/198/9/3375.long)</br>
[http://jinglebells.bgu.ac.il/](http://jinglebells.bgu.ac.il/)

##### Conquer

scRNA-seq database covering 40 datasets from 2013-2017, put together by Charlotte Soneson while she
was in Mark Robinson's lab at the University of Zurich. She completely reprocessed using Salmon.
Part of the process included removing samples from the original data, so I doubt I could use her
data directly, but it's very well done and documented. I believe the data is available, aligned both
to transcriptomes and to genomes.

[https://www.biorxiv.org/content/10.1101/143289v1.full](https://www.biorxiv.org/content/10.1101/143289v1.full)</br>
[http://imlspenticton.uzh.ch:3838/conquer/](http://imlspenticton.uzh.ch:3838/conquer/)

##### Human Cell Atlas Data Portal

Database containing the parsed scRNA-seq studies that have been compiled/submitted for the HCA
project. I think it depends on the pipeline, but it generally looks like they try to start from raw
FASTQ files. They use a mix of programmatic and human checks, then process and store the
data/metadata in a standard way. They are actively maintaining, and have 31 studies available, with
13 that have direct data download available.

[https://data.humancellatlas.org/](https://data.humancellatlas.org/)

##### Broad Single Cell Portal

Database where study authors can upload their scRNA-seq data. I don't think most of the data was
compiled by the database maintainers. It seems mostly like a place for authors to submit if they
choose. Looks like they have 142 studies. The UI is pretty clunky, but I think there's API access if
you link your Google account?

[https://singlecell.broadinstitute.org/single_cell](https://singlecell.broadinstitute.org/single_cell)

##### 10X Single Cell Gene Expression Datasets

10X has roughly 100 scRNA-seq datasets generated to demonstrate their physical chemistry or their
Cell Ranger software. It's possible that some of these are duplicated, but I haven't checked yet.
Unclear if these are available for batch download, but it doesn't look too complicated if I have to
do it manually. Some of the datasets are also the ones backing a paper:

Zheng et al, "Massively parallel digital transcriptional profiling of single cells"

[https://support.10xgenomics.com/single-cell-gene-expression/datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets)</br>
[https://www.nature.com/articles/ncomms14049](https://www.nature.com/articles/ncomms14049)
