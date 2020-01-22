# Single Cell Database Specification <a name="single_cell_database_specification"></a>

This document describes the current version of the single cell database (SCDB) specification,
including the schema of the data, the structure of the stored data, the process of retrieving data,
the process of storing data, and the future plans.

```
Version : WORKING 0.2
Date    : 2019-01-17
Author  : Nathan Fox
```

## Table of Contents <a name="table_of_contents"></a>

1. [Table of Contents](#table_of_contents)
2. [Database Conceptual Schema](#database_conceptual_schema)
    1. [Expression Data](#expression_data_schema)
    2. [Internal Universal Metadata](#internal_universal_metadata_schema)
        1. [Cell Specific Universal Metadata](#cell_specific_universal_metadata_schema)
        2. [Gene Specific Universal Metadata](#gene_specific_universal_metadata_schema)
    3. [Internal Author Annotated Metadata](#internal_author_annotated_metadata_schema)
    4. [External Metadata](#external_metadata_schema)
3. [Schema Implementation](#schema_implementation)
    1. [HDF5 File Structure](#hdf5_file_structure)
        1. [Expression Data](#expression_data_structure)
        2. [Internal Universal Metadata](#internal_universal_metadata_structure)
            1. [Cell-Specific Internal Universal Metadata](#cell_specific_internal_universal_metadata_structure)
            2. [Gene-Specific Internal Universal Metadata](#gene_specific_internal_universal_metadata_structure)
        3. [Internal Author Annotated Metadata](#internal_author_annotated_metadata_structure)
            1. [Cell-Specific Internal Author Annotated Metadata](#cell_specific_internal_author_annotated_metadata_structure)
            2. [Gene-Specific Internal Author Annotated Metadata](#gene_specific_internal_author_annotated_metadata_structure)
        4. [External Metadata](#external_metadata_structure)
4. [Appendix](#appendix)
    1. [Tissue List](#tissue_list)
    2. [Loom File Specification](#loom_file_specification)

## Database Conceptual Schema <a name="database_conceptual_schema"></a>

One entry in the SCDB is a single-cell RNA-seq dataset, composed of expression data for a
single set of cells and any cell/gene metadata to go with it. 

The SCDB stores 4 types of data:

1. Expression Data
2. Internal Universal Metadata
3. Internal Author-Annotated Metadata
4. External Metadata

### Expression Data <a name="expression_data_schema"></a>

Single-cell RNA-seq expression data is the main thing that the SCDB stores. These are gene x cell
expression matrices, where the value at `mat[i, j]` is the number of RNA molecules for gene `i` in
cell `j`. An entry may have multiple expression matrices, but they must all be the same dimensions
and have the same row/column labels. The count format of the values may be raw read counts or any
form of normalized values.

### Internal Universal Metadata <a name="internal_universal_metadata_schema"></a>

This is cell-specific or gene-specific metadata. These metadata follow a strict schema that must be
followed for every entry. This is to allow programmatic comparison of cell-specific or gene-specific
metadata across multiple datasets without having to reformat or accomodate different field names.

Here are the current schemata:

#### Cell Specific Universal Metadata <a name="cell_specific_universal_metadata_schema"></a>

All fields must exist for every entry. If a value is missing, it should be filled with "-1". If a
value does not fit the specification, it must be filled with "OTHER" (e.g. if a cell came from a
human with Klinefelter syndrome where they are XXY for sex chromosomes, or if a cell came from a
tissue not on the "Tissue List"). This also includes if a cell's information is known, but has
multiple values (e.g. if a cell came from a pooled sample of both male and female organisms). This
indicates that the data is available, but unusable, as opposed to data that is missing, but might
be filled in by inquiring with the authors.

All fields may optionally have a "description" attribute associated with them. This provides any
information needed to interpret the values of the field, not the generic description of the field
itself (e.g. For the "condition" field, an appropriate "description" would be "postnatal day" for
values {10, 10, 20, 20}, not "the experimental condition of the cell"). Additionally, the "batch"
field must have a "batch\_key" attribute associated with it. This is a string containing the names
of the batch variables included in this field. If this field has all "-1" missing values, this
string may be empty, but the attribute must exist anyway.

| Field | Description | Possible Values |
|:------|:------------|:----------------|
|CellID|Suggested by the loom file specification. Any unique ID identifying the cell within its entry.|Anything as long as every value is unique inside a single entry.|
|cluster|The most specific cell grouping reported by the authors. This should generally correspond to cell type. If this information is available as an annotated name, it should be used over a pure grouping.|Anything clearly distinguishing groups. (e.g. {1, 2, 3}). If an annotated name is available, it should be used over a pure grouping (e.g. {acinar, beta, vascular} over {1, 2, 3})|
|species|The species of origin. Must be the correctly capitalized scientific name. If there is additional information, such as a strain, it may be added after a "\|" (vertical bar) delimiter.|Any valid scientific name for a species in the format "Genus species". (e.g. "Mus musculus" or "Escherichia coli\|K-12")|
|tissue|The tissue of origin, taken from the "Tissue List", described in the Appendix. This field is not intended to parallel the "cluster" field. Rather, this is intended to be a high-level organizing field. This should be useful if multiple datasets are concatenated, or if an entry contains cells from multiple tissues.|Any value from the "Tissue List", described in the Appendix.|
|source\_organism|The specific organism of origin.|Anything clearly distinguishing different organisms. (e.g. "Mouse\_1")|
|sex|The sex of the organism of origin.|"M" for canonical cis-male and "F" for canonical cis-female.|
|condition|The experimental or disease condition of the cell. For example, if they came from multiple time points, disease states, or drug cohorts. If a cell is normal, healthy, and free of experimental condition, this value should be "NORMAL".|Anything clearly distinguishing multiple conditions. If a descriptive name is available, it should be used over a pure grouping (e.g. {"day0", "day10", "day20"}, not {"0, 1, 2"}. Normal, healthy cells free of experiment condition must be labeled "NORMAL".|
|batch|Any batch variable(s) not covered by "source\_organism", "sex", or "condition". Multiple remaining batch variables should be concatenated into a single string with a "\|" (vertical bar) delimiter. Additionally, this field must have a "batch\_key" string associated with it. This "batch\_key" must be the original variable names, also "\|" (vertical bar)-concatenated.|Any collection of batch variables, "\|"-concatenated into a single string. Each batch variable must be distinguishable within its own space, but there may be shared labels across multiple batch variables. (e.g. "1\|E10\|1" referring to plate 1, well E10, and sequencing run 1). This field must also store a "batch\_key" string containing the "\|"-concatenated variable names (e.g. "plate\|well\|sequencing\_run")|
|uuid|The Universally Unique Identifier (UUID) of the entry to which this cell belongs. This value should be identical for all cells from a single entry and is intended for when multiple entries are concatenated.|Any UUID in the database. (e.g. "6c582340-5f89-4ca1-a403-9a746f8e75b1")|

#### Gene Specific Universal Metadata <a name="gene_specific_universal_metadata_schema"></a>

An entry Must contain "Gene", "Accession", or both, but it may not omit both fields.

| Field | Description | Possible Values |
|:------|:------------|:----------------|
|Gene|Suggested by the loom file specification. Any human-readable gene ID. This field is not required to have unique values|Any human readable gene ID.|
|Accession|Suggested by the loom file specification. Any accession gene ID. This field must have unique values.|e.g. "ENSMUSG00000082108"|

### Internal Author Annotated Metadata <a name="internal_author_annotated_metadata_schema"></a>

This is also cell-specific and gene-specific metadata. This is an unedited dump of all the metadata
that the authors provided. This is to avoid losing any potential metadata, thus avoiding any need to
go back later and retrieve it. However, it is impossible to ensure that these metadata fields are
comparable across multiple entries. Different publications may provide different levels of
specificity, and building a consistent naming schema for every possible metadata field would not be
worth the effort. Not all of this metadata will be transformed/transformable into the internal
universal metadata, but it will still be stored in its original form if necessary.

No edits may be made to this metadata. It should be stored exactly as is. In general, cell-specific
metadata is more desired than gene\_specific metadata.

All fields may optionally have a "description" attribute associated with them. This provides any
information needed to interpret the values of the field, or a generic description of the field
itself (e.g. For a "well" field, an appropriate "description" would be "384-well plate well in
which the cell was processed").

### External Metadata <a name="external_metadata_schema"></a>

This is entry-specific metadata. These metadata describe each individual entry in the SCDB, not the
cells or genes inside each entry. If an entry does not have a value for a field, "-1" should be
stored. If an entry does not match the schema, or contains multiple values that are not allowed,
"OTHER" should be stored. If a field is Mandatory, a missing value may not be stored, but an
"OTHER" value is still allowed.

| Field | Mandatory | Description | Possible Values |
|:------|:----------|:------------|:----------------|
|species|Yes|The species of origin for the cells in this entry. Should match the format for the internal universal metadata field. Multiple species/strains may be included with a ";" (semicolon) delimiter.|e.g. "Homo sapiens" or "Mus musculus;Escherichia coli\|K-12"|
|tissue|No|The tissue of origin for the cells in this entry. Should match the format for the internal universal metadata field and be taken from the "Tissue List", described in the Appendix. Multiple tissues in a dataset may be included with a ";" (semicolon) delimiter.|e.g. "brain", "kidney;blood"|
|number\_of\_cells|Yes|The number of cells in this entry. Should be a non-negative integer.|e.g. "4028"|
|condition|No|Whether or not the cells in this entry fall under an experimental or disease condition. If any cells in this entry have anything other than "NORMAL" or "-1" in the internal universal metadata field "condition", this value should be "True".|"True" or "False"|
|date\_generated|Yes|The date that this entry was generated. Typically, the date on the publication reporting these data is used. Must be in YYYY-MM-DD format.|e.g. "2019-01-01"|
|count\_format|Yes|The format of the expression values in the expression matrix(ces) in this entry. Refers to the normalization state of the data. The number of values here must match the number of expression matrices in the entry. Multiple values should be ";" (semicolon) delimited.|Must be from this set: {"raw", "cpm", "tpm", "rpkm", "fpkm"}|
|umis|No|Whether or not Unique Molecular Identifiers (UMIs) were used in generating the data in this entry. |"True" or "False"|
|spikeins|No|Whether or not spike-in RNA molecules are reported in this entry. This is not an indication of whether or not they were used, only whether or not they are reported in the expression matrices.|"True" or "False"|
|technology|No|The chemical technology used to generate the expression data in this entry. Should be all lowercase and only alphanumeric. There is not a fixed set of values, but an effort must be made to keep a generally consistent format and not name the same technology multiple things.|e.g. "10xchromiumv2", "smartseq2", "celseq2", "fluidigmc1"|
|doi|No|The Digital Object Identifier (DOI) of the publication that presented the data in this entry. Must be in minimal format, not a URL.|e.g. "10.1038/s41467-019-13056-x"|
|accession|No|The accession number, if available, for the data in this entry. This could be a GEO Series ID, an ArrayExpress ID, or something similar. If the source database is not evident from the format of the accession ID, it may be indicated immediately after a "\|" (vertical bar) delimiter. Entries containing data from multiple accessions are discouraged, but allowed. Multiple entries should be ";" (semicolon) delimited.|e.g. "GSE108291", "E-MTAB-7195", "12345\|Obscure Database", "GSE10044;GSE10034"|
|date\_integrated|Yes|The date that this entry was entered into the SCDB. Must be in YYYY-MM-DD format.|e.g. "2019-01-01"|
|uuid|Yes|The Universally Unique Identifier (UUID) uniquely identifying this entry in the SCDB. Should be generated with the Python method `uuid.uuid4()`. Must be in the string format of a UUID: 12345678-1234-5678-1234-567812345678, where each of the 32 digits is 1 hexadecimal digit. This can be gotten with the Python code `str(uuid)` if `uuid` is a Python class `uuid.UUID`.|e.g. "6c582340-5f89-4ca1-a403-9a746f8e75b1"|
|file\_location|Yes|The full, absolute path to the folder holding this entry. Should not end with a "/" character.|e.g. "/data/single\_cell\_database/6c582340-5f89-4ca1-a403-9a746f8e75b1"|
|internal|Yes|Whether or not the data in this entry were generated by the Gillis Lab. This includes data generated by collaborators explicitly for us. This does not include data generated by others that we then began collaborating with after the fact.|"True" or "False"|

## Schema Implementation <a name="schema_implementation"></a>

This is a description of how the above SCDB schema is actually implemented.

The database is stored in a single folder on the Gillis Lab server, `tyrone`, at
`/data/single_cell_database/`. Inside that folder, there are 2 files and `n` folders, where `n` is
the number of entries in the database. There should be no other files or folders in this
root SCDB folder.

The 2 files are:

| File | Description |
|:-----|:------------|
|`external_metadata.tsv`|A tab-separated value text file holding the external metadata as a table. Fields are columns, in the order described above. Each row is a single entry. There are column headers, but no row labels.|
|`notes.tsv`|A tab-separated value text file holding notes about each entry, presumably if they require further attention. Each line has two values. The UUID and a variable-length string holding the note. This is a very loose file and not explicitly part of the schema.|

Each folder contains the expression matrices and internal metadata for a single entry. The folders
are named the UUID associated with the entry, found in the external metadata. So the output of an
`ls *` shell command when SCDB contains one entry could look like this:

```
$ pwd
/data/single_cell_database
$ ls *
external_metadata.tsv  notes.tsv

2d963abc-a990-45b1-a9ff-8d474eaa9649:
expr_mat.loom
```

Each folder must contain a file called `expr_mat.loom`. This file contains the expression matrices
and all internal metadata for that entry. It is an HDF5 file that conforms to the loom file
specification, version 3.0, as published by Sten Linnarsson's lab, albeit with a few small
additions that don't violate the loom file specification rules. The loom file specification is
reproduced in the appendix.

The basic summary is that an HDF5 file is organized and accessible as a file tree
inside the object. All items in an HDF5 file are either Groups or Datasets. Groups can hold other
Groups or Datasets. Datasets can hold an data array of a single datatype without labels. Both Groups
and Datasets can have small scalar or array pieces of data attached as Attributes.

### HDF5 File Structure <a name="hdf5_file_structure"></a>

This is an example of the HDF5 structure of an `expr_mat.loom` file. HDF5 Groups are designated
with a "/" (forward slash) at the end of the name. The mapping of this file structure to the SCDB
schema is described below. Throughout, the `expr_mat.loom` file will be referred to as `lfile`.

```
/
├── attrs/
│   └── LOOM_SPEC_VERSION
├── matrix
├── layers/
│   └── tpm
├── col_attrs/
│   ├── CellID
│   ├── cluster
│   ├── species
│   ├── tissue
│   ├── source_organism
│   ├── sex
│   ├── condition
│   ├── batch
│   └── uuid
├── col_graphs/
├── cell_author_annot/
│   ├── Age
│   ├── SampleID
│   └── Experimental_Day
├── row_attrs/
│   ├── Gene
│   └── Accession
├── row_graphs/
└── gene_author_annot/
    ├── highly_variable_gene
    └── astrocyte marker
```

#### Expression Data <a name="expression_data_structure"></a>

Each expression matrix is stored as a single HDF5 Dataset. Every expression matrix must be stored
with genes in the rows and cells in the columns. The raw read counts, if available, should
be stored in `lfile['/matrix']`. If not available, any other expression matrix may be stored in
`lfile['/matrix']`. The expression matrix stored in `lfile['/matrix']` MUST be the first one
described in the external metadata "count_format" field (e.g. if the "count_format" field contains
"tpm;cpm", then the "tpm" expression matrix must be stored at `lfile['/matrix']`. All other
expression matrices should be stored as HDF5 Datasets under `lfile['/layers/']`. They must be named
their corresponding value in the external metadata "count_format" field. However, if they are
designated as "OTHER" in the external metadata "count_format" field, they may be named anything
sensible that could be identified by reading the original publication, but is not one of the
pre-determined values (e.g. it could be named "norm_500" but not "cpm"). Expression matrices stored
under `lfile['/layers/']` cannot be named "matrix". Expression matrices May be named "counts" or
"logcounts", but this is HEAVILY discouraged because of potential confusing imports into
R SingleCellExperiment objects.

#### Internal Universal Metadata <a name="internal_universal_metadata_structure"></a>

##### Cell-Specific Internal Universal Metadata <a name="cell_specific_internal_universal_metadata_structure"></a>

The internal universal metadata for cells is stored in `lfile['/col_attrs/']`. Each field is stored
as a separate HDF5 Dataset, named the name of the field. Each Dataset must be a 1-D array with the
same length as the number of cells in the entry. If present, the field "description" attribute may
be stored as an HDF5 Attribute named "description" attached to the relevant HDF5 Dataset. It must be
a single scalar string containing the "description". The "batch\_key" attribute should be stored as
an HDF5 Attribute named "batch\_key" attached to the `lfile['/col_attrs/batch']` HDF5 Dataset.

##### Gene-Specific Internal Universal Metadata <a name="gene_specific_internal_universal_metadata_structure"></a>

The internal universal metadata for genes is stored in `lfile['/row_attrs/']`. Each field is stored
as a separate HDF5 Dataset, named the name of the field. Each Dataset must be a 1-D array with the
same length as the number of genes in the entry.

#### Internal Author Annotated Metadata <a name="internal_author_annotated_metadata_structure"></a>

##### Cell-Specific Internal Author Annotated Metadata <a name="cell_specific_internal_author_annotated_metadata_structure"></a>

**Note:** This is an addition to the loom file format.

The internal author annotated metadata for cells is stored in `lfile['/cell_author_annot/']`. Each
field is stored as a separate HDF5 Dataset, named the name of the field. Each Dataset must be a 1-D
array with the same length as the number of cells in the entry. If present, the field "description"
attribute may be stored as an HDF5 Attribute named "description" attached to the relevant HDF5
Dataset. It must be a single scalar string containing the "description".

##### Gene-Specific Internal Author Annotated Metadata <a name="gene_specific_internal_author_annotated_metadata_structure"></a>

**Note:** This is an addition to the loom file format.

The internal author annotated metadata for genes is stored in `lfile['/gene_author_annot/']`. Each
field is stored as a separate HDF5 Dataset, named the name of the field. Each Dataset must be a 1-D
array with the same length as the number of genes in the entry. If present, the field "description"
attribute may be stored as an HDF5 Attribute named "description" attached to the relevant HDF5
Dataset. It must be a single scalar string containing the "description".

#### External Metadata <a name="external_metadata_structure"></a>

External metadata is not stored in any `expr_mat.loom` files. It is stored in a TSV file as
described in the table above.

## Accessing the SCDB <a name="accessing_the_scdb"></a>

Access to the SCDB is provided in both Python and R as a set of access functions. Most are available
in both languages, but a few are language-specific. They are described below.

|Function Signature|Language|Description|
|:------------|:-------|:----------|
|`get_extern_md()`|Python, R|Get the external metadata as a dataframe.|
|`uuid_to_row(uuid, columns=None)`|Python, R|Get the row of the external metadata that corresponds to the UUID.|
|`get_loom_filename(uuid)`|Python, R|Get the full path to the `expr_mat.loom` file that corresponds to the UUID.|
|`get_h5_conn(uuid)`|Python, R|Get an open HDF5 connection to the requested `expr_mat.loom` file. h5py in Python and hdf5r in R.|
|`get_loom_conn(uuid)`|Python, R|Get an open loom connection to the requested `expr_mat.loom` file. loompy in Python and loomR in R.|
|`get_expr_mat(uuid, matrix='matrix')`|Python, R|Get an expression matrix as a matrix from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_expr_mat_names(uuid)`|Python, R|Get a list of all the expression matrix names from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_cell_univ(uuid, keep_missing=True)`|Python, R|Get the internal cell-specific universal metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_gene_univ(uuid, keep_missing=True)`|Python, R|Get the internal gene-specific universal metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_cell_author_annot(uuid)`|Python, R|Get the internal cell-specific author annotated metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_gene_author_annot(uuid)`|Python, R|Get the internal gene-specific author annotated metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_batch_key(uuid)`|Python, R|Get the "batch\_key" attribute for the internal cell-specific universal metadata field "batch" from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_cell_ids(uuid)`|Python, R|Get the cell IDs as a vector from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_gene_ids(uuid, accession=False)`|Python, R|Get the gene IDs as a vector from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_column_description(uuid, column, var='cell', metadata='universal')`|Python, R|Get the "description" attribute for the requested column from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_anndata(uuid, **kwargs)`|Python|Get the entire entry corresponding to the UUID as an AnnData object.|
|`get_sce(uuid, assay_for_matrix='counts', counts_assay=NULL, logcounts_assay=NULL`|R|Get the entire entry corresponding to the UUID as a SingleCellExperiment object.|

They all have documentation that will not be described here.

### Dependencies <a name="dependencies"></a>

The following non-standard libraries are dependencies for the access functions described above:

|Language|Dependencies|
|:-------|:-----------|
|Python|numpy, pandas, h5py, loompy, scanpy|
|R|hdf5r, loomR, SingleCellExperiment|

## Adding Entries to the SCDB <a name="adding_entries_to_the_scdb"></a>

**Note:** Regular users may not edit any files in the SCDB. This section is included as a matter
of record and to communicate how easily/quickly new entries may be added. All files in the SCDB are
read-only to all users except a local account called `scdb`. This is a local user on all `rugen`
nodes and on `tyrone`, but with a synchronized UID to enable permissions across NFS mounts.

Database management is all in Python 3 and relies on a library of functions written specifically for
this project.

Currently, adding a new entry happens in an interactive Python session. It is relatively trivial
and takes <10 minutes if the data is available in the following format:

1. Expression matrices as gene x cell numpy arrays.
2. Cell and gene IDs as numpy arrays.
3. Author annotated metadata as Pandas dataframes.

Creating the universal metadata is done manually in Pandas from the author annotated metadata.

There are single functions that accomplish the following tasks:

1. Given an expression matrix, cell IDs, and gene IDs, create a new loom file and store the data.
2. Given a Pandas dataframe, store it as universal or author annotated metadata for a given UUID.
3. Given a Pandas dataframe and a column list, construct a "batch" column and "batch\_key".
4. Given a string, set it as the "description" attribute for a given column.
5. Run an interactive command line application to set a new entry in the external metadata. Supports
   pre-filling fields (e.g. programmatically calculating the "date\_integrated") and provides the
   capacity to restart, cancel, or get field descriptions on the fly.

Right now, all other architecture is built to parse GEO deposited datasets into the above format.
To that end, there are single functions that accomplish the following tasks.

1. Given a GEO Series ID, generate a new UUID, create a new folder in the database, and download the
   SOFT metadata file, as well as all supplementary files. The supplementary files typically include
   the expression data.
2. Process the output of a 10X Cell Ranger pipeline into an expression matrix, cell IDs, and gene
   IDs numpy arrays.
3. Process a CSV-like file into an expression matrix, cell IDs, and gene IDs numpy arrays.
4. Parse a SOFT metadata file for single cell characteristics and store in a Pandas dataframe.

Together, this means that adding a new GEO Series-sourced entry typically takes ~20 minutes and
looks like this:

1. Start an interactive Python session.
2. Set the GSE ID.
3. Run a function to create a new entry and download/unzip all files.
4. Examine the downloaded files in a separate shell using vim and other Unix utilities. Make any
   edits necessary to feed the files to my parsers (e.g. renaming files, checking for headers).
5. Run a function to process the files into an expression matrix, cell ID list, and a gene ID list.
6. Run a function to create a new loom file and store the above in it.
7. If available, repeat 5 and 6 to parse and store other expression matrices.
8. Run 2 functions to parse the SOFT file for cell-specific author annotated metadata into a
   dataframe.
9. Check the original publication for supplementary data files containing internal cell-specific
   and gene-specific author annotated metadata. These are often Excel files that have to be manually
   downloaded on my local computer, exported into TSV files, and pushed to my `rugen` node.
10. Run a function to import any author annotated metadata into Pandas dataframes.
11. Run a function to store author annotated metadata.
12. Manually construct the internal universal metadata dataframes using Pandas and numpy. This
    usually looks like reformatting certain columns from the author annotated metadata.
13. Run a function to construct the "batch" field and "batch\_key" attribute.
14. Run a function to store the internal universal metadata.
15. Run a function to store internal metadata field "description" attributes if necessary.
16. Run a function to store the external metadata.
17. Add a `notes.tsv` entry if there is missing data to be requested from the authors. Send an
    email if necessary.

About two thirds of this process is straightforward and accessible to anyone with light
computational experience. The other third is the manual curation, where speed and accuracy heavily
depends on the programmer's skill with Pandas and with Unix utilities. Additionally, there is a
large amount of error-handling code, but it is not comprehensive and definitely needs to be
revisited. There is very little validation code, i.e. code that validates adherence to the schema
described above.

## Future Directions <a name="future_directions"></a>

This is an informal list of features, functionality, and code that should be written in future
versions of the SCDB and this specification.

* A semi-universally comparable language for cell types that is more specific than the Tissue List.
  For example, this would enable programmatic comparison of many datasets of the same tissue. It
  would be interesting to see if cell types within a tissue are replicable across datasets.
* Input validation code. This would ensure adherence with the above schema. For example, right now
  any dataframes submitted as cell-specific universal metadata are checked to verify that they are
  the right shape and that all the columns are correctly named. But the "sex" field is not verified
  to ensure that all its values are "M", "F", "-1" or "OTHER".
* Maintenace validation code. This is code that runs regularly to ensure that data corruption hasn't
  occurred. The external metadata could include a "checksum" field that holds a hashsum value for
  the corresponding `expr_mat.loom` file. That way if anything changes unexpectedly in the file,
  e.g. from corrupted data or an accidental overwrite, an alert is raised. There could be a check to
  make sure that all entries in the external metadata have corresponding folders. There could be
  checks to ensure that the external metadata fields that come from the internal metadata agree.
* Infrastructure to speed up adding entries from non GEO Series sources. This should not be written
  until needed, but will likely be needed eventually.

## Appendix <a name="appendix"></a>

### Tissue List <a name="tissue_list"></a>

This is a compromise for the need to have universally comparable tissue labels across entries.
Originally, this was going to be a subtree of the Cell Ontology (everything below "organ" using
both "is\_a" and "part\_of" relationships), but that turned out to be more complicated than was
worth it at this stage of development. So instead, Nathan manually read the subtree and pulled out
his best approximation of a rank-equal list of tissues. This will be used until a better solution is
created. Fields that use this Tissue List still have the option to store "OTHER". All terms should
be used verbatim. Categories are not eligible as terms, and are only presented here for clarity.

|Tissue Category|Tissue Terms|
|:--------------|:-----------|
|Nervous|brain, spinal cord, peripheral nervous system|
|Respiratory|lung, trachea|
|Digestive|tongue, esophagus, stomach, liver, gall bladder, pancreas, small intestine, large intestine|
|Cardiovascular|heart, vasculature, blood|
|Immune|spleen, thymus, lymph node|
|Urinary|urethra, kidney, urinary bladder|
|Musculoskeletal|ligament, tendon, skeletal muscle, smooth muscle, bone|
|Gland|thyroid gland|
|Reproductive|uterus, placenta, ovary, clitoris, vagina, testis, penis|
|Sensory|ear, eye, nose|
|Skin|skin|

### Loom File Specification <a name="loom_file_specification"></a>

The following is copied directly from Sten Linnarsson's lab website:

**<< BEGINNING OF LOOM FILE SPECIFICATION >>**

#### Loom file format specs

##### Versions

This specification defines the loom file format version `3.0.0`.

##### Introduction

The Loom file format is designed to efficiently hold large omics datasets. Typically, such data
takes the form of a large matrix of numbers, along with metadata for the rows and columns. For
example, single-cell RNA-seq data consists of expression measurements for all genes (rows) in a
large number of cells (columns), along with metadata for genes (e.g. `Chromosome`, `Strand`,
`Location`, `Name`), and for cells (e.g. `Species`, `Sex`, `Strain`, `GFP positive`).

We designed Loom files to represent such datasets in a way that treats rows and columns the same.
You may want to cluster both genes and cells, you may want to perform PCA on both of them, and
filter based on quality controls. SQL databases and other data storage solutions almost always treat
data as a table, not a matrix, and makes it very hard to add arbitrary metadata to rows and columns.
In contrast, Loom makes this very easy.

Furthermore, current and future datasets can have tens of thousands of rows (genes) and hundreds of
thousands of columns (cells). We designed Loom for efficient access to arbitrary rows and columns.

The annotated matrix format lends itself to very natural representation of common analysis tasks.
For example, the result of a clustering algorithm can be stored simply as another attribute that
gives the cluster ID for each cell. Dimensionality reduction such as PCA or t-SNE, similarly, can be
stored as two attributes giving the projection coordinates of each cell.

Finally, we recognize the importance of graph-based analyses of such datasets. Loom supports graphs
of both the rows (e.g. genes) and the columns (e.g. cells), and multiple graphs can be stored each
file.

##### HDF5 concepts

The Loom format is based on HDF5, a standard for storing large numerical datasets. Quoting from
h5py.org:

> An HDF5 file is a container for two kinds of objects: datasets which are array-like collections
> of data, and groups, which are folder-like containers that hold datasets and other groups. The
> most fundamental thing to remember when using h5py is: *Groups work like dictionaries, and
> datasets work like NumPy arrays*.

A valid Loom file is simply an HDF5 file that contains specific *groups* containing the main matrix
as well as row and column attributes. Because of this, Loom files can be created and read by any
language that supports HDF5, including Python, R, MATLAB, Mathematica, C, C++, Java, and Ruby.

##### Standards

The official MIME media type for a loom file is `application/vnd.loom`, approved by IANA. You should
use this media type when requesting a file using HTTP, or if you create a server that offers Loom
files for download.

See the IANA record for application/vnd.loom for important security considerations.

The Loom file format is designated as a standard format for gene expression matrices by the Global
Alliance for Genomics and Health (GA4GH) standards body, as part of the rnaget API. The rnaget API
documentation contains good examples of how to use the media type specification.

##### Specification

A valid Loom file conforms to the following:

###### Main matrix and layers

* There MUST be a single HDF5 dataset at `/matrix`, of dimensions (N, M)
* There can OPTIONALLY be an HDF5 group `/layers` containing additional matrices (called "layers")
* Each additional layer MUST have the same (N, M) shape
* Each layer can have a different data type, compression, chunking, etc.

###### Global attributes

* There MUST be an HDF5 group `/attrs` containing global attributes.
* There MUST be an HDF5 dataset `/attrs/LOOM_SPEC_VERSION` with the value `3.0.0`

Global attributes apply semantically to the whole file, not any specific part of it. Such attributes
are stored in the HDF5 group `/attrs` and can be any valid scalar or multidimensional datatype.

As of Loom file format v3.0.0, only one global attribute is mandatory: the `LOOM_SPEC_VERSION`
attribute, which is a string value giving the loom file spec version that was followed in creating
the file. See top of this document for the current version of the spec.

Note: previous versions of the loom file format stored global attributes as HDF5 attributes on the
root `/` group. However, such attributes are size-limited, which caused problems for some
applications.

###### Row and column attributes

* There MUST be a group `/row_attrs`
* There can OPTIONALLY be one or more datasets at `/row_attrs/{name}` whose first dimension has
  length N
* There MUST be a group `/col_attrs`
* There can OPTIONALLY be one or more datasets at `/col_attrs/{name}` whose first dimension has
  length M

The datasets under `/row_attrs` should be semantically interpreted as row attributes, with one value
per row of the main matrix, and in the same order. Therefore, all datasets under this group must
be arrays with exactly N elements, where N is the number of rows in the main matrix.

The datasets under `/col_attrs` should be semantically interpreted as column attributes, with one
value per column of the main matrix, and in the same order. Therefore, all datasets under this
group must be arrays with exactly M elements, where M is the number of columns in the main matrix.

###### Row and column sparse graphs

* There MUST be a group `/col_graphs`
* There can OPTIONALLY be one or more groups at `/col_graphs/{name}`
* Under each `/col_graphs/{name}` group, there MUST be three one-dimensional datasets called `a`
  (integer), `b` (integer) and `w` (float). These should be interpreted as a sparse graph in
  coordinate list format. The lengths of the three datasets MUST be equal, which defines the number
  of edges in the graph. Note that the number of columns in the dataset defines the vertices, so an
  unconnected vertex is one that has no entry in `a or b`.
* There MUST be a group `/row_graphs`
* There can OPTIONALLY be one or more groups at `/row_graphs/{name}`
* Under each `/row_graphs/{name}` group, there MUST be three one-dimensional datasets called `a`
  (integer), `b` (integer) and `w` (float). These should be interpreted as a sparse graph in
  coordinate list format. The lengths of the three datasets MUST be equal, which defines the number
  of edges in the graph. Note that the number of rows in the dataset defines the vertices, so an
  unconnected vertex is one that has no entry in `a` or `b`.
* Vertex indexing is zero-based. When an entry in `a` or `b` is zero, this denotes the first column
  in the matrix. If there are N columns, then vertices are numbered from 0 to N - 1.

##### Datatypes

The main matrix and additional layers MUST be two-dimensional arrays of one of these numeric types:
`int8`, `int16`, `int32`, `int64`, `uint8`, `uint16`, `uint32`, `uint64`, `float16`, `float32` and
`float64`. Each layer can have its own datatype.

Row and column attributes are multidimensional arrays whose first dimension matches the
corresponding main matrix dimension. The elements MUST be of one of the numeric datatypes `int8`,
`int16`, `int32`, `int64`, `uint8`, `uint16`, `uint32`, `uint64`, `float16`, `float32` and `float64`
or fixed-length ASCII strings.

Global attributes are scalars or multidimensional arrays of any shape, whose elements are any of the
numeric datatypes `int8`, `int16`, `int32`, `int64`, `uint8`, `uint16`, `uint32`, `uint64`,
`float16`, `float32` and `float64` or fixed-length ASCII strings.

Starting with v3.0.0 of the spec, all strings are stored as variable-length UTF-8 encoded.

Note: in previous version, strings were stored as fixed-length null-padded 7-bit ASCII. Unicode
characters outside 7-bit ASCII were stored using XML entity encoding, to ensure maximum
compatibility. Strings were decoded when read and encoded when written.

##### Example

Here's an example of the structure of a valid Loom file:

|Group|Type|Description|
|:----|:---|:----------|
|/attrs/|(subgroup)|Global attributes|
|/matrix|float32[N,M] or uint16[N,M]|Main matrix of N rows and M columns|
|/layers/|(subgroup)|Subgroup of additional matrix layers|
|/row\_attrs/|(subgroup)|Subgroup of all row attributes|
|/row\_attrs/Name|string[N]|Row attribute "Name" of type string|
|/row\_attrs/Species|string|Row attribute "Species" of type string|
|/col\_attrs/|(subgroup)|Subgroup of all column attributes|
|/col\_attrs/CellID|float64[M]|Column attribute "CellID" of type float64|
|/col\_graphs/|(subgroup)|Subgroup of all column graphs|
|/col\_graphs/KNN|(subgroup)|A column graph "KNN"|
|/col\_graphs/KNN/a|int32[E]|Vector of edge 'from' vertices|
|/col\_graphs/KNN/b|int32[E]|Vector of edge 'to' vertices|
|/col\_graphs/KNN/w|float64[E]|Vector of edge weights|
|/row\_graphs/|(subgroup)|Subgroup of all row graphs|

###### Backwards compatibility

Loom v3.0.0 introduces two major backwards-incompatible changes (global attributes and
variable-length strings; see above).

A compliant Loom reader MUST check the LOOM\_SPEC\_VERSION and treat files consistently with their
spec. For example, when writing a global attribute, the writer MUST write only to the `/attrs` group
if `LOOM_SPEC_VERSION` is `3.0.0` or higher. The writer MUST write the HDF5 attributes on the root
`/` group if `LOOM_SPEC_VERSION` is lower than `3.0.0` or if it does not exist. This is to preserve
a consistent format for legacy files.

##### Conventions

In order to maximize interoperability of Loom files between analysis pipelines, we suggest adhering
to the following conventions.

**Note:** This document is work in progress and subject to change! You can follow the
[discussion here](https://github.com/linnarsson-lab/loompy/issues/19).

###### Single analysis per file

Each Loom file stores a single analysis. If you want to try two different ways of clustering, you
store the results separately.

This convention simplifies things a lot, because we only need to keep track of one set of cluster
labels (for example). It also means that we don't need to store the relationship between
different attributes (e.g. that *this clustering* was done using *that PCA*). Such relationships
must be stored external to the file.

###### Orientation

* Columns represent cells or aggregates of cells
* Rows represent genes

Loom files can grow along the column axis, but not the row axis, so this makes sense.

###### Attribute naming conventions

We propose that algorithms always accept an argument specifying the name of each attribute it will
use, with the defaults set as listed below. For example, an algorithm that needs a unique gene
accession string would take an argument `accession_attr="Accession"`, and a clustering algorithm
that generates cluster labels would take an argument `cluster_id_attr="ClusterID"`. In this way,
Loom files that conform to the conventions below would work without fuss, but files that don't
could still be made to work by supplying the non-standard attribute names.

###### Column attributes

`CellID` a string label unique to each cell (preferably globally unique across all datasets)

`Valid` integers 1 or 0, indicating cells that are considered valid after some QC step

`ClusterID` an integer label with values in [0, n\_clusters].

`ClusterName` a string label with arbitrary values representing cluster names. If both `ClusterID`
and `ClusterName` are present, they should correspond `:`.

`Outliers` an integer label 1 or 0, indicating cells that are outliers relative to the clusters.

Any attribute that is an M-by-Y matrix (where M is the number of columns and Y is the dimensionality
of the embedding) can be used to store e.g. a PCA or t-SNE dimensionality reduction.

###### Row attributes

`Gene` a string human-readable gene name, not necessarily unique

`Accession` a string, unique within the file (e.g. an ENSEMBL accession etc.)

`Selected` integers 1 or 0, indicating genes that were selected by some previous step.

`Valid` integers 1 or 0, indicating genes that are considered valid after some QC step.

**<< END OF LOOM FILE SPECIFICATION >>**
