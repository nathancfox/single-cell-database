# SCDB Tutorial <a name=scdb_tutorial></a>

```
Version : WORKING 0.3
Date    :
Author  : Nathan Fox
```

## Table of Contents <a name=table_of_contents></a>

1. [Tutorial](#tutorial)
2. [Function Index](#function_index)
3. [Using a Different Version](#using_a_different_version)

## Tutorial <a name="tutorial"></a>

This is a tutorial for using the Single Cell Database (SCDB). The SCDB is a database of
scRNA-seq datasets. Before we get started, we should load the access code.

### Setup Code

#### R

The R access code is currently not available as a formal package. Instead, you'll have to load
all the functions into your current environment.

To use these access functions, you need to install the following dependencies:

* hdf5r
* loomR
* SingleCellExperiment
* httr

```
> source("/data/single_cell_database/src/current/scdb/access.R")
```

**NOTE:** If you are on a non-`rugen` node, you may have to change the path to reflect your
correct path to the `tyrone` NFS mount on your system.

#### Python

To use these access functions, you need to install the following dependencies:

* numpy
* pandas
* h5py
* loompy
* scanpy
* requests

The Python access code is importable as a package. However, you first have to update your
`PYTHONPATH` environment variable so that the Python interpreter can find it. To do this,
first run this bash command:

```
$ source /data/single_cell_database/src/setup -v=current
```

**NOTE:** If you are on a non-`rugen` node, you may have to change the path to reflect your
correct path to the `tyrone` NFS mount on your system.

To avoid having to run the setup script in the future, you can add this line to your `.bashrc`
file.

Then, start a Python REPL and import the package:

```
>>> import scdb
```

### Database Layout

The database entries are split into two parts:

1. External Metadata
2. Expression Data and Internal Metadata

The external metadata is stored in a TSV file as a regular table. It has the following columns:

| Field | Mandatory | Description | Possible Values |
|:------|:----------|:------------|:----------------|
|title|No|The title of the publication that reported these data.|All characters must be lowercase alphanumeric, " "(whitespace), "\_"(underscore), or "-"(hyphen)|
|authors|No|The author(s) of the publication that reported these data.|A ", "(comma whitespace)-delimited list of the authors. No capitalization or punctuation rules inside an author's name. All author names must be in the format "GIVEN\_NAME FAMILY\_NAME". The first author should be first, but all other order is not guaranteed.|
|abstract|No|The abstract of the publication that reported these data.|Must be a single string with no line breaks or leading/trailing whitespace.|
|species|Yes|The species of origin for the cells in this entry. Should match the format for the internal universal metadata field. Multiple species/strains may be included with a ";" (semicolon) delimiter.|e.g. "Homo sapiens" or "Mus musculus;Escherichia coli\|K-12"|
|tissue|No|The tissue of origin for the cells in this entry. Should match the format for the internal universal metadata field and be taken from the "Tissue List", described in the Appendix. Multiple tissues in a dataset may be included with a ";" (semicolon) delimiter.|e.g. "brain", "kidney;blood"|
|number\_of\_cells|Yes|The number of cells in this entry. Should be a non-negative integer.|e.g. "4028"|
|author\_clusters|Yes|Whether or not the author(s) provided cell groupings or clusters. This should be "True" even if the clusters are unannotated.|"True" or "False"|
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
|gillis\_lab|Yes|Whether or not the data in this entry were generated by the Gillis Lab. This includes data generated by collaborators explicitly for us. This does not include data generated by others that we then began collaborating with after the fact.|"True" or "False"|

**NOTE:** All columns can have OTHER as a value if the truth does not match the acceptable
values. Additionally, non-mandatory columns can have "-1" as a value if the information is
missing.

This external metadata can be gotten as a dataframe with the following method:

#### R

```
> emdf <- get_extern_md()
> str(emdf)
'data.frame':   17 obs. of 18 variables:
 $ title          : chr  "variation in activity state axonal projection..."
 $ authors        : chr  "Maxime Chevée, Johanna De Jong Robertson, ..."
 $ abstract       : chr  "Single-cell RNA sequencing has generated catalogs..."
 $ species        : Factor w/ 1 level "Mus musculus": 1 1 ...
 $ tissue         : Factor w/ 3 levels "brain","skin",..: 1 1 ...
 $ number_of_cells: int  1023 5454 2303 24185 3186 1504 250 2669 ...
 $ author_clusters: Factor w/ 2 levels "False","True": 2 2 ...
 $ condition      : Factor w/ 2 levels "False","True": 2 2 ...
 $ date_generated : Date, format: "2018-01-09" "2018-01-15" ...
 $ count_format   : Factor w/ 6 levels "fpkm", "OTHER",..: 2 3 ...
 $ umis           : Factor w/ 2 levels "False","True": 1 2 ...
 $ spikeins       : Factor w/ 3 levels "False","OTHER"",..: 2 1
 $ technology     : Factor w/ 6 levels "10chromiumv2",..: 6 2
 $ doi            : chr  "10.1016/j.celrep.2017.12.046" "10.1038/s41593-017-0056-2" ...
 $ accession      : chr  "GSE107632" "GSE95315" "GSE95752" "GSE104323" ...
 $ date_integrated: Date, format: "2020-01-02" "2020-01-03" ...
 $ uuid           : chr  "ec344eaf-bde2-4555-aa42-163d4b373914" ...
 $ gillis_lab     : Factor w/ 1 level "False": 1 1 ...
```

**NOTE:** The columns that look like they should be booleans are left as factors to accomodate
potential values of "OTHER" or "-1".

#### Python

```
>>> emdf = scdb.get_extern_md()
>>> emdf.iloc[0]
title               variation in activity state axonal projection ...
authors             Maxime Chevée, Johanna De Jong Robertson, Gabr...
abstract            Single-cell RNA sequencing has generated catal...
species                                                  Mus musculus
tissue                                                          brain
number_of_cells                                                  1023
author_clusters                                                  True
condition                                                        True
date_generated                                    2018-01-09 00:00:00
count_format                                                    OTHER
umis                                                            False
spikeins                                                        OTHER
technology                                                  smartseq2
doi                                      10.1016/j.celrep.2017.12.046
accession                                                   GSE107632
date_integrated                                   2020-01-02 00:00:00
uuid                             ec344eaf-bde2-4555-aa42-163d4b373914
gillis_lab                                                      False
Name: 0, dtype: object
```

Each entry in the database has a row in this table. All other access functions take UUIDs
as lookups. You can use the dataset-specific metadata in this table to identify the datasets
you want to use, then use their UUIDs to retrieve them.

Additionally, you also have a lookup function, `uuid_to_row()` to extract a dataset's external
metadata using its UUID.

#### R

```
> str(uuid_to_row(emdf[["uuid"]][1]))
'data.frame':   1 obs. of 18 variables:
 $ title          : chr  "variation in activity state axonal projection..."
 $ authors        : chr  "Maxime Chevée, Johanna De Jong Robertson, ..."
 $ abstract       : chr  "Single-cell RNA sequencing has generated catalogs..."
 $ species        : Factor w/ 1 level "Mus musculus": 1
 $ tissue         : Factor w/ 3 levels "brain","skin",..: 1
 $ number_of_cells: int  1023
 $ author_clusters: Factor w/ 2 levels "False","True": 2
 $ condition      : Factor w/ 2 levels "False","True": 2
 $ date_generated : Date, format: "2018-01-09"
 $ count_format   : Factor w/ 6 levels "fpkm", "OTHER",..: 2
 $ umis           : Factor w/ 2 levels "False","True": 1
 $ spikeins       : Factor w/ 3 levels "False","OTHER"",..: 2
 $ technology     : Factor w/ 6 levels "10chromiumv2",..: 6
 $ doi            : chr  "10.1016/j.celrep.2017.12.046"
 $ accession      : chr  "GSE107632"
 $ date_integrated: Date, format: "2020-01-02"
 $ uuid           : chr  "ec344eaf-bde2-4555-aa42-163d4b373914"
 $ gillis_lab     : Factor w/ 1 level "False": 1
```

#### Python
```
>>> scdb.uuid_to_row(emdf.loc[0, 'uuid'])
title               variation in activity state axonal projection ...
authors             Maxime Chevée, Johanna De Jong Robertson, Gabr...
abstract            Single-cell RNA sequencing has generated catal...
species                                                  Mus musculus
tissue                                                          brain
number_of_cells                                                  1023
author_clusters                                                  True
condition                                                        True
date_generated                                    2018-01-09 00:00:00
count_format                                                    OTHER
umis                                                            False
spikeins                                                        OTHER
technology                                                  smartseq2
doi                                      10.1016/j.celrep.2017.12.046
accession                                                   GSE107632
date_integrated                                   2020-01-02 00:00:00
uuid                             ec344eaf-bde2-4555-aa42-163d4b373914
gillis_lab                                                      False
Name: 0, dtype: object
```

For the purposes of this tutorial, let's use a small mouse brain dataset.

#### R

```
dset_id <- dplyr::filter(emdf, number_of_cells < 3000
                             & species == "Mus musculus"
                             & tissue == "brain")[4, "uuid"]
> dset_id
[1] "19e00811-d167-45ce-a9f9-67142e213195"
```

#### Python

```
>>> dset_id = emdf[(emdf['number_of_cells'] < 3000)
                   & (emdf['species'].str.match('Mus musculus'))
                   & (emdf['tissue'].str.match('brain')), 'uuid'].iloc[3]
>>> dset_id
'19e00811-d167-45ce-a9f9-67142e213195'
```

Everything except the external metadata is stored in an individual HDF5 file. This HDF5 file
conforms to the Loom file format, albeit with a few additions. None of the additions should
invalidate the file for any external functions expecting to receive a standard Loom file.

HDF5 files are internally structured like a file tree. Each item can either be a Group
(like a file directory) or a Dataset (like a file). Each Dataset is an n-dimensional array
of a single datatype. Additionally, all Groups and Datasets can have one or more Attributes.
These are small pieces of data attached to any Group or Dataset. The Loom file has the
following internal structure:

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

`/matrix` is a 2-D Dataset that holds the first expression matrix. This is a genes x cells matrix
of transcript counts. If there are other expression matrices of the **same** data, but in a
different normalization state, they are stored under the `/layers` Group. The Dataset under
`/matrix` will always be the first one in the list in the "count_format" external metadata
column. All other layers will be named their corresponding "count_format" value unless they
are "OTHER". In this case, they will be named something short and descriptive.

`/col_attrs` holds "internal universal cell-specific metadata". This is cell-specific (as 
opposed to dataset specific) metadata that is universally named and parseable across multiple
entries in the database. Each Dataset in this Group stores the data for a column in a table.
The fields are described below. Like the non-mandatory external metadata columns, all columns
can have an "OTHER" or "-1" value. `CellID` is mandatory. All others may be missing.

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

`/row_attrs` holds "internal universal *gene*-specific" metadata. This is the parallel of
`/col_attrs` but for genes instead of cells. Each dataset must have either `Accession` and/or
`Gene`, but cannot be missing both.

| Field | Description | Possible Values |
|:------|:------------|:----------------|
|Accession|Suggested by the loom file specification. Any accession gene ID. This field must have unique values.|e.g. "ENSMUSG00000082108"|
|Gene|Suggested by the loom file specification. Any human-readable gene ID. This field is not required to have unique values|Any human readable gene ID.|

`/cell_author_annot` is a Group similar to `/col_attrs`. However, it holds "internal
author-annotated cell-specific metadata". This is a raw dump of all the cell-specific
metadata that the authors provided with this dataset. There is no consistency in these metadata
between multiple entries in the database. The intention of this Group is to preserve all
metadata that the author(s) provided, especially metadata that could not be parsed into
the "universal metadata". It does not have any predetermined fields because it is just a
dump of all provided cell-specific metadata. The Datasets in the HDF5 structure above
are examples.

`/gene_author_annot` is the analogue of `/cell_author_annot` for genes.

`/col_graphs` and `/row_graphs` are Groups required by the Loom file format, but the SCDB
does not use them.

#### Attributes

There are a few Attributes used by the SCDB.

|Attribute|Description|
|:--------|:----------|
|`all_missing`|All internal metadata columns have an `all_missing` Attribute that is True if all values in that column are "-1" (missing) and False otherwise. This enables constant-time lookup checks for missing internal metadata columns, allowing faster filtering for desired datasets.|
|`batch_key`|The `batch` column in the internal universal cell-specific metadata is a concatenation of 1 or more columns. The `/col_attrs/batch` Dataset has a `batch_key` Attribute that stores the names of the concatenated columns.|
|`column_order`|HDF5 files do not remember order of creation/addition for Datasets inside a Group. The `/cell_author_annot` and `/gene_author_annot` Groups both have a `column_order` attribute that saves the order of columns in the table constructed from each group.|
|`description`|Additionally, any internal metadata Dataset can have a `description` Attribute that explains the values. For example, the `condition` column in the internal universal cell-specific metadata might only contain the value "14.5". Because the values are not self-explanatory, the `description` attribute for this column might contain the string "embryonic day age". These descriptions are not guaranteed.|

### Accessing Datasets

Once you have the UUIDs of the datasets you want to analyze, you can access the data in a number of ways.

You can get an HDF5 or Loom file connection object for any dataset. This is not encouraged,
as it requires a detailed knowledge of HDF5 to be effective, but low-level access functions
are provided:

#### R

```
> hfile <- get_h5_conn(dset_id)

WARNING: R and Python interpret binary matrices
in transposed ways. R uses column-major order and
Python uses row-major order. Because the database
entries are created in Python, this means that all
n-D HDF5 datasets (where n > 1) accessed manually
from an HDF5 or loom connection will be returned
transposed, i.e. cells x genes. However, cell-specific
internal universal metadata will still be in the
"col_attrs" HDF5 group and vice versa for genes.

To silence this warning, pass the function arg:
    "warning = FALSE"

> hfile
Class: H5File
Filename: /data/single_cell_database/database/ec344eaf-bde2-4555-aa42-163d4b373914/expr_mat.loom
Access type: H5F_ACC_RDONLY
Attributes: last_modified
Listing:
              name    obj_type dataset.dims dataset.type_class
             attrs   H5I_GROUP         <NA>               <NA>
 cell_author_annot   H5I_GROUP         <NA>               <NA>
         col_attrs   H5I_GROUP         <NA>               <NA>
        col_graphs   H5I_GROUP         <NA>               <NA>
 gene_author_annot   H5I_GROUP         <NA>               <NA>
            layers   H5I_GROUP         <NA>               <NA>
            matrix H5I_DATASET 1023 x 31749        H5T_INTEGER
         row_attrs   H5I_GROUP         <NA>               <NA>
        row_graphs   H5I_GROUP         <NA>               <NA>
> lfile <- get_loom_conn(dset_id, warning=FALSE)
```

**NOTE:** This obnoxious warning is there because R and Python have an incompatibility when
handling HDF5 files. In binary form, R reads and writes matrices in row-major order (walking
down rows and changing columns). However, Python reads and writes matrices in column-major
order (walking down columns and changing rows). This means that if R is used to read a
matrix stored in binary form by Python, R will read a transposed matrix relative to how Python
originally stored it. All R access functions handle this behind the scenes, but low-level
R access functions will throw this warning unless you deliberately mute it to prevent any
tricky bugs when you retrieve expression matrices by hand.

#### Python

```
>>> hfile = scdb.get_h5_conn(dset_id)
>>> hfile.keys()
<KeysViewHDF5 ['attrs', 'cell_author_annot', 'col_attrs', 'col_graphs', 'gene_author_annot', 'layers', 'matrix', 'row_attrs', 'row_graphs']>
>>> lfile = scdb.get_loom_conn(dset_id)
```

If you want to read datasets using these low-level connections, see the hdf5r/loomR or the
h5py/loompy packages respectively. For now, let's close the file connections and use the
more accessible functions.

#### R

```
> hfile$close_all()
> lfile$close()
```

#### Python

```
>>> hfile.close()
>>> lfile.close()
```

All elements of the HDF5 file can be accessed via getter methods. Let's see what expression
matrices are stored in our chosen dataset.

#### R

```
> get_expr_names(dset_id)
[1] "matrix" "rpkm"
```

#### Python

```
>>> scdb.get_expr_mat_names(dset_id)
array(['matrix', 'rpkm'], dtype='<U6')
```

There are two available forms of the expression data. We can double check the external metadata
to see what is in the "matrix" layer.

#### R

```
> as.character(uuid_to_row(dset_id)[["count_format"]])
[1] "row;rpkm"
```

#### Python

```
>>> scdb.uuid_to_row(dset_id)['count_format']
'raw;rpkm'
```

So now we know that the "matrix" layer holds raw counts and the "rpkm" layer holds the same
data, but normalized to Reads Per Kilobase Million. Let's get the raw counts and look at
the first few values.

#### R

```
> mat <- get_expr_mat(dset_id, matrix = "matrix")
> mat[1:10, 1:10]
10 x 10 sparse Matrix of class "dgCMatrix"
   [[ suppressing 10 column names 'C1-101-A10_CGAGGCTG-GCGTAAGA_L008_R1_all', ...]]

ENSMUSG00000000001   . .    .   1    .    .  57  .    .    .
ENSMUSG00000000003   . .    .   .    .    .   .  .    .    .
ENSMUSG00000000028   . . 1439  31    .    1   .  .    .  626
ENSMUSG00000000031 189 1   84 151 2231 2161 128  1    . 3504
ENSMUSG00000000037   . .    .   .    .    .   .  .    .    .
ENSMUSG00000000049   . .    .   .    .    .   .  .    .    .
ENSMUSG00000000056 193 .    .   1    .    .   . 38    .    .
ENSMUSG00000000058   . .    .   .    .    .   .  .    .    .
ENSMUSG00000000078  13 .    .   .    .    .  20  .    .    .
ENSMUSG00000000085   . .    2  22    .    .   .  . 1932    .
```

#### Python

```
>>> mat = scdb.get_expr_mat(dset_id, matrix = "matrix")
>>> mat
<37310x2669 sparse matrix of type '<class 'numpy.int64'>'
        with 8301426 stored elements in Compressed Sparse Column format>
>>> mat[:10, :10].toarray()
array([[   0,    0,    0,    1,    0,    0,   57,    0,    0,    0],
       [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
       [   0,    0, 1439,   31,    0,    1,    0,    0,    0,  626],
       [ 189,    1,   84,  151, 2231, 2161,  128,    1,    0, 3504],
       [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
       [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
       [ 193,    0,    0,    1,    0,    0,    0,   38,    0,    0],
       [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
       [  13,    0,    0,    0,    0,    0,   20,    0,    0,    0],
       [   0,    0,    2,   22,    0,    0,    0,    0, 1932,    0]])
```

These are `Matrix::dgCMatrix` and `scipy.sparse.csc_matrix` sparse matrices respectively.

You can also look at the cell and gene IDs.

#### R

```
> cell_ids <- get_cell_ids(dset_id)
> cell_ids[1:5]
[1] "C1-101-A10_CGAGGCTG-GCGTAAGA_L008_R1_all"
[2] "C1-101-A1_TAAGGCGA-GCGTAAGA_L008_R1_all"
[3] "C1-101-A4_TCCTGAGC-GCGTAAGA_L008_R1_all"
[4] "C1-101-A5_GGACTCCT-GCGTAAGA_L008_R1_all"
[5] "C1-101-A6_TAGGCATG-GCGTAAGA_L008_R1_all"
> gene_names <- get_gene_ids(dset_id, accession=FALSE)
> gene_names[1:5]
[1] "GNAI3" "PBSN"  "CDC45" "H19"   "SCML2"
> gene_ids <- get_gene_ids(dset_id, accession=TRUE)
> gene_ids[1:5]
[1] "ENSMUSG00000000001" "ENSMUSG00000000003" "ENSMUSG00000000028"
[1] "ENSMUSG00000000031" "ENSMUSG00000000037"
```

#### Python

```
>>> cell_ids = scdb.get_cell_ids(dset_id)
array(['C1-101-A10_CGAGGCTG-GCGTAAGA_L008_R1_all'
       'C1-101-A1_TAAGGCGA-GCGTAAGA_L008_R1_all'
       'C1-101-A4_TCCTGAGC-GCGTAAGA_L008_R1_all'
       'C1-101-A5_GGACTCCT-GCGTAAGA_L008_R1_all'
       'C1-101-A6_TAGGCATG-GCGTAAGA_L008_R1_all'], dtype=object)
>>> gene_names = scdb.get_gene_ids(dset_id, accession = False)
>>> gene_names[:5]
array(['GNAI3', 'PBSN', 'CDC45', 'H19', 'SCML2'], dtype=object)
>>> gene_ids = scdb.get_gene_ids(dset_id, accession = True)
>>> gene_ids[:5]
array(['ENSMUSG00000000001', 'ENSMUSG00000000003', 'ENSMUSG00000000028'
       'ENSMUSG00000000031', 'ENSMUSG00000000037'], dtype=object)
```

The way that gene IDs are handled is a bit complicated. As you read above, every dataset must
have `Accession` and/or `Gene` in its internal universal gene-specific metadata. This means
that any dataset may have one or the other or both. Although either field can be retrieved
with the above functions, whenever rownames or other unique gene IDs are needed, the first
valid option from the following list is used.

1. `Accession`, if available.
2. `Gene`, if available **and** does not contain duplicates.
3. An integer index.

Any of the internal metadata tables can also be accessed as dataframes. Additionally, you
can choose to automatically drop columns that have all missing values from the universal
metadata.

#### R

```
cell_univ <- get_cell_univ(dset_id, keep_missing = False)
cell_aa <- get_cell_author_annot(dset_id)
gene_univ <- get_gene_univ(dset_id)
gene_aa <- get_gene_author_annot(dset_id)
```

#### Python

```
cell_univ = scdb.get_cell_univ(dset_id, keep_missing = False)
cell_aa = scdb.get_cell_author_annot(dset_id)
gene_univ = scdb.get_gene_univ(dset_id)
gene_aa = scdb.get_gene_author_annot(dset_id)
```

For cell-specific internal metadata, the rownames of the returned dataframes are taken from
the `CellID` column and the `CellID` column is dropped. For gene-specific metadata, the
rownames of the returned dataframes are chosen using the above ruleset, but no columns
are dropped (unless you pass `keep_missing = False`). If a table does not exist (e.g. if
the author did not provide any gene-specific metadata other than gene IDs), then an
empty dataframe is returned.

Additionally, the `description` and `all_missing` attributes can be retrieved for any column.

#### R

```
> get_column_allmissing(dset_id, "Gene", var = "gene", metadata = "universal")
[1] FALSE
> get_column_description(dset_id, "cluster", var="cell", metadata="author_annot")
No description available!
NULL
```

#### Python

```
>>> scdb.get_column_allmissing(dset_id, "Gene", var = "gene", metadata = "universal")
False
>>> scdb.get_column_description(dset_id, "cluster", var="cell", metadata="author_annot")
No description available!
```

Finally, each dataset can be retrieved as a popular format for each language. In R, each
dataset can be retrieved as a `SingleCellExperiment` object. The universal metadata is
stored in `rowData` and `colData`. The rownames are assigned based on the ruleset above.
The colnames are assigned from the `Cell ID` column. The author-annotated metadata and
the `batch_key` Attribute are placed in the `S4Vectors::metadata()` list. Assay naming
is somewhat complicated. This is an excerpt from the function documentation.

> The loom file has a main matrix in the HDF5 dataset 'matrix'.
The rest of the matrices, if available, are under the
HDF5 group 'layers', as datasets of their own. SingleCellExperiment
objects typically have the raw counts stored under the assay "counts",
and the log-normalized counts stored under the assay "logcounts".
Other matrices can be named anything. The `assay_for_matrix`
argument should hold the name of the assay that the 'matrix' dataset
should be stored in. Then, the `counts_assay` and
`logcounts_assay` arguments can hold the names of HDF5 layers
that should be stored in "counts" or "logcounts" assays respectively.
For example: if the loom file has a raw reads matrix in 'matrix'
and a normalized matrix in 'norm' and a 3rd matrix in 'misc', the
appropriate function call would be:
```
get_sce(UUID,
        assay_for_matrix = "counts",
        counts_assay = NULL,
        logcounts_assay = "norm")}
```
> However, if the 3rd matrix was in 'matrix' and the raw reads and
normalized reads were in 'raw' and 'norm', the appropriate function
call would be:
```
get_sce(UUID,
        assay_for_matrix = "misc",
        counts_assay = "raw",
        logcounts_assay = "norm")}
```

There is also an R function that attempts to naively merge two `SingleCellExperiment` objects
returned by the `get_sce()` function. This is also a bit complicated, and is best explained by
another function documentation excerpt.

```
Merge SingleCellExperiment objects returned from the database.

Attempts to merge 2 SingleCellExperiment objects returned from
get_sce(). Will reject the merge if any incompatibilities
exist.

@param sce1 SingleCellExperiment object. The first
  SingleCellExperiment object to be merged. Should be an
  object returned from get_sce().
@param sce2 SingleCellExperiment object. The second
  SingleCellExperiment object to be merged. Should be an
  object returned from get_sce().
@param keep_author_annot Logical vector of length 1. If TRUE,
  the internal author-annotated metadata from each SCE will
  be labeled with cell_id_prefix and saved in the metadata
  of the new SCE. Otherwise, they will be dropped.
@param min_common_genes Numeric vector of length 1. Must
  be an non-negative integer. The two SCE objects must
  have this many valid gene IDs in common. Otherwise, the
  merge will fail.
@param cell_id_prefix Can be a single integer from the
  following list c(4, 8, 12, 16, 20, 24, 28, 32) or
  a character vector of length 2 with 2 different prefixes
  for the cell labels from sce1 and sce2
  respectively. If it is an integer, that number of digits
  from the UUIDs are used as a prefix (does not include the
  hyphens). If the number is too small to produce unique
  prefixes, it is increased by 4 until it works.
@return A new SingleCellExperiment object containing all
  cells from both objects (sce1 first). It only
  retains genes that were guaranteed to be unique inside
  each SCE object and in common between the two objects.
  Cell IDs are prefixed with cell_id_prefix. The
  colData() from each is also merged. The names of
  the author-annotated internal metadata dataframe, if
  kept are prefixed with the cell_id_prefix. The
  function attempts to combine the batch_keys into a 
  vector, but may not do this correctly if sce1 or sce2
  are the product of a previous merge.
```

#### R

```
> sce <- get_sce(dset_id, assay_for_matrix = "counts", logcounts_assay = "rpkm")
> sce
class: SingleCellExperiment
dim: 37310 2669
metadata(3): batch_key cell_author_annot gene_author_annot
assays(2): counts logcounts
rownames(37310): ENSMUSG00000000001 ENS MUSG0000000003 ...
  ENSMUSG00000093373 ENSMUSG00000093374
rowData names(2): Accession Gene
colnames(2669): C1-101-A10_CGAGGCTG-GCGTAAGA_L008_R1_all
  C1-101-A1_TAAGGCGA-GCGTAAGA_L008_R1_all ...
  C1-190-2-H8_CAGAGAGG-CTAAGCCT_L001_R1_all
  C1-190-2-H9_GCTACGCT-CTAAGCCT_L001_R1_all
colDAta names(i): cluster species ... batch uuid
reducedDimNames(0):
spikeNames(0):
altExpNames(0):
> dset_id2 <- dplyr::filter(emdf, number+of_cells < 3000
                                  & species == "Mus musculus"
                                  & tissue == "brain")[8, "uuid"]
> uuid_to_row(dset_id_2)[["count_format"]]
[1] raw
Levels: fpkm OTHER raw raw;OTHER raw;rpkm rpkm
> sce2 <- get_sce(dset_id_2, assay_for_matrix = "counts")
> sce2
class: SingleCellExperiment
dim: 27998 1700
metadata(3): batch_key cell_author_annot gene_author_annot
assays(2): counts
rownames(27998): ENSMUSG00000051951 ENSMUSG00000089699 ...
  ENSMUSG00000096730 ENSMUSG00000095742
rowData names(2): Accession Gene
colnames(250): GSM2677817_AAAGATGAGACCTAGG-1
  GSM2677817_AAAGATGCATCCTAGA-1 ... GSM2677819_TTTGTCAAGGTGACCA-1
  GSM2677819_TTTGTCATCATGTCTT-1
colData names(8): cluster species ... batch uuid
reducedDimNames(0):
spikeNames(0):
altExpNames(0):
> merged_sce <- merge_sce(sce, sce2,
                          keep_author_annot = TRUE,
                          min_common_genes = 15000,
                          cell_id_prefix = 4)
Warning message:
In merge_sce(sce, sce2) : Dropped unmergeable assay "logcounts" from sce1!
> merged_sce
class: SingleCellExperiment
dim: 22876 4369
metadata(5): 19e0__cell_author_annot 19e0__gene_author_annot
  fea1__cell_author_annot fea1__gene_author_annot batch_key
assays(1): counts
rownames(22876): ENSMUSG00000000001 ENSMUSG00000000003 ...
  ENSMUSG00000092622 ENSMUSG00000092625
rowData names(2): Accession Gene
colnames(4369): 19e0__C1-101-A10_CGAGGCTG-GCGTAAGA_L008_R1_all
  19e0__C1-101-A1_TAAGGCGA-GCGTAAGA_L008_R1_all ...
  fea1__GSM2677819_TTTGTCAAGGTGACCA-1
  fea1__GSM2677819_TTTGTCATCATGTCTT-1
colData names(8): cluster species ... batch uuid
reducedDimNames(0):
spikeNames(0):
altExpNames(0):
```

In Python, each dataset can be retrieved as an `anndata` object. This is also a bit complex.
Here is the function docstring.

```
Get an AnnData object from a dataset.

Given a UUID, loads the associated loom file
completely into memory into a scanpy AnnData object.
Unlike most other access functions, this function
can edit the schema. To facilitate pipelining
into other tools, the user can choose to concatenate
the universal and author_annot internal metadata
into the adata.var and adata.obs dataframes instead
of leaving the author_annot in adata.uns.

Args:
    uuid: String. UUID of the desired dataset
    keep_missing: String. Must be 'cells', 'genes', 'both',
        or 'none'. Indicates which internal universal metadata
        columns to keep. If 'cells, then cell-specific internal
        universal metadata columns that are all_missing will
        still be kept, but gene-specific missing columns will
        be dropped. Vice versa for 'genes'. If 'both', then
        all columns will be kept. If 'none', then missing columns
        will be dropped from both cell-specific and gene-specific
        internal universal metadata. Does not apply to
        author-annotated metadata.
    combine: String. Must be 'cells', 'genes', 'both',
        or 'none'. Indicates which internal metadata to
        combine. If 'cells', the cell-specific internal
        metadata will be combined into a single dataframe,
        stored under adata.obs. Universal columns will have
        the prefix "scdb_" added to their names. Vice versa
        for 'genes'. If 'both', then both will be combined.
        If 'none', then neither will be combined and the
        author_annot internal metadata will be left in
        adata.uns.
    **kwargs: Keyword arguments passed to
        scanpy.read_loom()

Returns:
    AnnData object holding the dataset from the
    desired loom file.

Raises: None
```

#### Python

```
>>> adata = scdb.get_anndata(dset_id,
                             keep_missing = 'both',
                             combine = 'none')
>>> adata
AnnData object with n_obs x n_vars = 2669 x 37310
    obs: 'cluster', 'species', 'tissue', 'source_organism', 'sex', 'condition', batch', 'uuid'
    var: 'Accession', 'Gene'
    uns: 'batch_key', 'gene_author_annot', 'cell_author_annot'
    layers: 'rpkm'
```

## Function Index <a name=function_index></a>

Detailed documentation does not exist in document form yet. All functions have in-code
documentation.

In R, these are roxygen function docstrings. However, because these have not been built into a
formal R package, these must be accessed by reading the source code file that was originally
imported. Again, this is `/data/single_cell_database/src/current/scdb/access.R`

In Python, all functions have Google-style docstrings. These can be conveniently accessed
from the interactive Python REPL. In the Python REPL, this is `help(scdb.get_anndata)` and in the
iPython REPL, this is `scdb.get_anndata?`.

|Function Signature|Language|Description|
|:------------|:-------|:----------|
|`get_extern_md()`|Python, R|Get the external metadata as a dataframe.|
|`uuid_to_row(uuid, columns=None)`|Python, R|Get the row of the external metadata that corresponds to the UUID.|
|`get_h5_conn(uuid)`|Python, R|Get an open HDF5 connection to the requested `expr_mat.loom` file. h5py in Python and hdf5r in R.|
|`get_loom_conn(uuid)`|Python, R|Get an open loom connection to the requested `expr_mat.loom` file. loompy in Python and loomR in R.|
|`get_expr_mat(uuid, matrix='matrix')`|Python, R|Get an expression matrix as a sparse matrix from the `expr_mat.loom` file that corresponds to the UUID. Note that this function has a small difference in the Python vs the R implementation.|
|`get_expr_mat_names(uuid)`|Python, R|Get a list of all the expression matrix names from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_cell_univ(uuid, keep_missing=True)`|Python, R|Get the internal cell-specific universal metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_gene_univ(uuid, keep_missing=True)`|Python, R|Get the internal gene-specific universal metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_cell_author_annot(uuid)`|Python, R|Get the internal cell-specific author annotated metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_gene_author_annot(uuid)`|Python, R|Get the internal gene-specific author annotated metadata as a dataframe from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_batch_key(uuid)`|Python, R|Get the "batch\_key" attribute for the internal cell-specific universal metadata field "batch" from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_cell_ids(uuid)`|Python, R|Get the cell IDs as a vector from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_gene_ids(uuid, accession=True)`|Python, R|Get the gene IDs as a vector from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_column_allmissing(uuid, column, var='cell', metadata='universal')`|Python, R|Get the "all_missing" attribute for the requested column from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_column_description(uuid, column, var='cell', metadata='universal')`|Python, R|Get the "description" attribute for the requested column from the `expr_mat.loom` file that corresponds to the UUID.|
|`get_note(uuid)`|Python, R|Get the note from the "notes.tsv" file for the dataset that corresponds to the UUID.|
|`get_bibtex(uuid)`|Python, R|Get a BibTeX citation for the entry corresponding to the UUID.|
|`export_bibtex(uuids, outfile, append=False)`|Python, R|Export a collection of BibTeX citations to a file.|
|`get_anndata(uuid, keep_missing='both', combine='none', **kwargs)`|Python|Get the entire entry corresponding to the UUID as an AnnData object.|
|`get_sce(uuid, assay_for_matrix='counts', counts_assay=NULL, logcounts_assay=NULL`|R|Get the entire entry corresponding to the UUID as a SingleCellExperiment object.|
|`merge_sce(sce1, sce2, keep_author_annot=FALSE, min_common_genes=15000, cell_id_prefix=4)`|R|Attempt to merge 2 SingleCellExperiment objects returned from `get_sce()`.|

## Using a Different Version <a name=using_a_different_version></a>

These access functions are versioned. By using the `current` version in the `source()` (R) and setup
script (Python), you automatically loaded the last stable version. You can get a list of all
available versions with the setup script on the command line.

```
$ source /data/single_cell_database/src/setup -l
current -> v0.3
devel
v0.2
v0.3
```

The `devel` version is the development version. It is unstable and may or may not be working. It can

If you want to use a different version of the code, simply change the version in the R `source()`
call or the Python setup script:

#### R

```
> source("/data/single_cell_database/src/VERSION/scdb/access.R")
```

#### Python

```
$ source /data/single_cell_database/src/setup -v=VERSION
```

**NOTE:** Versioning may be abandoned in the future. It is not very helpful because older versions
are sometimes not future-compatible. If the schema of the database changes (e.g. if HDF5 Groups get
renamed), then the entire database has to be retroactively changed. Because only the access code is versioned and not the actual database, older access code that worked fine in the past will not work for newer versions of the database. Thus, old versions are not guaranteed to be
working even if they were stable in the past.
