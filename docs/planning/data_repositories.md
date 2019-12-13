# Data Repositories

Description of the data repositories that I will be pulling data from and their relationships. Also
some basic preliminary stats on the data I'll be yanking.

### ArrayExpress

> ArrayExpress is one of the repositories recommended by major scientific journals to archive
> functional genomics data from microarray and sequencing platforms to support reproducible
> research. To serve this mission we facilitate submissions in compliance with Minimum Information
> About A Microarray Experiment (MIAME) and Minimum Information About A Sequencing Experiment
> (MINSEQE) guidelines. Experiments are submitted directly to ArrayExpress or are imported from the
> NCBI Gene Expression Omnibus database. For high-throughput sequencing experiments the raw data is
> brokered to the European Nucleotide Archive, while the experiment descriptions and processed data
> are archived in ArrayExpress.

A first search of the metadata returned:

146 Experiments submitted directly to ArrayExpress
3 Experiments imported from GEO

Of the 146, all but 2 are in Mage-TAB format. The other three are in ERAD. All 149 experiments came
from the search `exptype:"from single cells"`.

### European Nucleotide Archive (ENA)

> The European Nucleotide Archive (ENA) provides a comprehensive record of the world's nucleotide
> sequencing information, covering raw sequencing data, sequence assembly information, and
> functional annotation.

The ENA stores high-throughput sequencing data, analogously to SRA and DDBJ. All three are part of
the International Nucleotide Sequence Database Collaboration (INSDC). Data submitted to any of the
three organizations are shared among them. Looks like, for scRNA-seq experiments, will be almost
entirely FASTQ files with raw reads. No expression matrices. Maybe alignments? A first search of
the database with:

`rna seq ("scRNA-seq" OR "single cell" OR "single-cell")`

Returned 1413 results under "Study" and 1512 under "Study (Sequence)". I don't know what the
difference is, or if one is inclusive of the other.

### Gene Expression Omnibus (GEO)

> GEO is a public functional genomics data repository supporting MIAME-compliant data submissions.
> Array- and sequence-based data are accepted. Tools are provided to help users query and download
> experiments and curated gene expression profiles.

Likely to be the largest source of data for this project. Run by the NCBI. Individual datasets
generally are represented as GEO Series (GSEXXXXX). Series contain samples, which for scRNA-seq
typically represent single cells. GEO takes both raw data (FASTQ files), which it hosts in SRA, or
processed data like expression matrices, which are typically stored as supplementary files. Metadata
and ftp links to supplementary files are programmatically accessible via Entrez and BuildaURL
functionality. A first search of the database with:

`(("gse"[Entry Type]) AND "expression profiling by high throughput sequencing"[DataSet Type]) AND
("single cell" OR "single-cell" OR "scRNA-seq")`

returned 2601 results.

### Sequence Read Archive (SRA)

> Sequence Read Archive (SRA) makes biological sequence data available to the research community to
> enhance reproducibility and allow for new discoveries by comparing data sets. The SRA stores raw
> sequencing data and alignment information from high-throughpu sequencing platforms, including
> Roche 454 GS System, Illumina Genome Analyzer, Applied Biosystems SoLiD System, Helicos Heliscope,
> Complete Genomics, and Pacific Biosciences SMRT.

The SRA stores high-throughput sequencing data, analogously to ENA and DDBJ. All three are part of
the International Nucleotide Sequence Database Collaboration (INSDC). Data submitted to any of the
three organizations are shared among them.

A first search of the Studies section of the SRA: `scRNA-seq` returned 254 studies. 

### European Genome-Phenome Archive (EGA)

> The EGA provides a service for the permanent archiving and distribution of personally identifiable
> genetic and phenotypic data resulting from biomedical reserach projects.  Data at EGA was
> collected from individuals whose consent agreements authorise data release only for specific
> research use to bona fide researches. Strict protocols govern how information is managed, stored,
> and distributed by the EGA project.

Looks like it's only disease phenotypes, potentially only in humans. Their self reported stats show
about 95 single cell datasets.

### China National GeneBank (CNGB) Sequence Archive (CNSA)

> CNGB Sequence Archive (CNSA) is a convenient and fast online submission system for biological
> research projects, samples, experiments, runs, assemblies, variations, and other information data.
> Based on the International Nucleotide Sequence Database Collaboration (INSDC) standard and
> DataCite standard, accepting the submission of global scientific research sequencing data,
> information, and analysis result data, its data submission service can be used as a supplement to
> the literature publishing process to support early data sharing. CNSA is committed to the storage
> and sharing of biological sequencing data, information, and analysis result data, and is designed
> to provide global researchers with the most comprehensive data and information resources, enabling
> researches to access data more easily and facilitate data reuse.

Note that CNSA is not part of the INSDC, but it serves a similar purpose. Not sure about its
reputation. A quick first search of "Projects" with `scRNA-seq` returned 13,381 results. This seems
highly unlikely that there are that many usable datasets. This is likely a poor search query.
