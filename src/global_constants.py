"""A collection of global constants.

This is a collection of global constants that reflect the schema
of the database implementation. These SHOULD NOT be edited in any
code throughout the entire project. This is intended as a reference
for all other Python code in the project.

NOTE: These constants are Python-readable, but are not directly
referenced by the R access functions. Some of them are mirrored
manually in the access.R script. They should always be up-to-date
with the constants here.

Descriptions
============

get_PATH_TO_DATABASE():
    A method returning the full, absolute path to the root folder
    of the database. It is only a function because tyrone is
    mounted as different names on different servers.
get_PATH_TO_METADATA():
    A method returning the full, absolute path to the external
    metadata file. It is only a function because tyrone is
    mounted as different names on different servers.
_EM_COLUMN_DESCRIPTIONS:
    A dict holding printable descriptions for each of the columns
    in the external metadata schema.
_EM_COLUMN_INDEX:
    A dict holding the 0-based index for each column in the external
    metadata schema, indicating column order.
_EM_COLUMN_MANDATORY:
    A dict holding a boolean for each column in the external metadata
    schema, indicating whether or not the column is mandatory.
_IMU_CELL_COLUMN_DESCRIPTIONS:
    A dict holding printable descriptions for each of the columns
    in the internal universal cell-specific metadata schema.
_IMU_CELL_COLUMN_INDEX:
    A dict holding the 0-based index for each column in the internal
    universal cell-specific metadata schema, indicating column order.
_IMU_CELL_COLUMN_MANDATORY:
    A dict holding a boolean for each column in the internal universal
    cell-specific metadata schema, indicating whether or not the
    column is mandatory.
_IMU_GENE_COLUMN_DESCRIPTIONS:
    A dict holding printable descriptions for each of the columns
    in the internal universal gene-specific metadata schema.
_IMU_GENE_COLUMN_INDEX:
    A dict holding the 0-based index for each column in the internal
    universal gene-specific metadata schema, indicating column order.
_IMU_GENE_COLUMN_MANDATORY:
    A dict holding a boolean for each column in the internal universal
    gene-specific metadata schema, indicating whether or not the
    column is mandatory.
_TISSUE_LIST:
    A dict of lists, holding the list of tissues used in the tissue
    metadata fields. Each list is named an organ system category
    and holds the tissue names relevant to it.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
from socket import getfqdn
# These 2 are mirrored in the access.R file.
# If you ever change these, MAKE SURE TO UPDATE
# THE ONES IN THE access.R FILE!
def get_PATH_TO_DATABASE():
    if getfqdn() == 'dactyl.cshl.edu':
        return('/tyronedata/single_cell_database/database')
    else:
        return('/data/single_cell_database/database')
def get_PATH_TO_METADATA():
    if getfqdn() == 'dactyl.cshl.edu':
        return('/tyronedata/single_cell_database/database/external_metadata.tsv')
    else:
        return('/data/single_cell_database/database/external_metadata.tsv')

_EM_COLUMN_DESCRIPTIONS = {
    'title':
        (
            '01. Title\n'
            '    The title of the publication that published the data.\n'
            '    e.g. \"Hedgehog stimulates hair follicle neogenesis by\n'
            '          creating inductive dermis during murine skin\n'
            '          wound healing.\"'
        ),
    'authors':
        (
            '02. Authors\n'
            '    The authors of the publication that published the data.\n'
            '    This should be a single string containing a list of authors,\n'
            '    delimited by \", \". The first author must be first, but\n'
            '    the rest may be in any order. The author names must be in\n'
            '    the format: \"GIVEN_NAME FAMILY NAME\"\n'
            '    e.g. \"Chae Ho Lim, Qi Sun, Karan Ratti, Mayumi Ito'
        ),
    'abstract':
        (
            '03. Abstract\n'
            '    The abstract of the publication that published the data.\n'
            '    This should be a single string with no line breaks or\n'
            '    leading/trailing whitespace.'
        ),
    'species': 
        (
            '04. Species\n'
            '    The species that the cells originated from. Should be a\n'
            '    two word scientific name in standard format. If there is\n'
            '    additional information, such as a string, it may be added\n'
            '    after a \"|\" delimiter. You can have multiple species/strains\n'
            '    in a dataset. These should be \";\" delimited.\n'
            '    e.g. \"Homo sapiens\", \"Mus musculus;Escherichia coli|K-12\"'
        ),
    'tissue':
        (
            '05. Tissue\n'
            '    The tissue/organ that the cells originated from. Must be from\n'
            '    the TISSUE_LIST. Enter \"LIST\" to view the options. Multiple\n'
            '    tissues in a dataset can be \";\" delimited.\n'
            '    e.g. \"brain\", \"kidney\", \"blood\"'
        ),
    'number_of_cells':
        (
            '06. Number of Cells\n'
            '    The number of cells in this dataset. Should just be a\n'
            '    simple, non-negative integer.\n'
            '    e.g. \"5402\"'
        ),
    'author_clusters':
        (
            '07. Author Clusters\n'
            '    Whether or not the author provided cell clusters. Must be\n'
            '    one of these two values: {\"True\", \"False\"}. Should be\n'
            '    \"True\" if cell clusters are provided by the author. They do\n'
            '    not have to be annotated, as long as some type of cell\n'
            '    grouping is provided. This can be \"True\", even if all cells\n'
            '    are in the same group, as long as this aligns with the\n'
            '    author\'s report. Should be \"False\" otherwise.'
        ),
    'condition':
        (
            '08. Condition\n'
            '    Whether or not these cells fall under a disease condition\n'
            '    or an experimental condition. Must be one of these two\n'
            '    values: {\"True\", \"False\"}. Should be False if all values\n'
            '    for condition in the internal cell-specific universal metadata\n'
            '    are \"NORMAL\" or \"-1\". Should be \"True\" otherwise.\n'
            '    This could be \"False\" if the organism is unhealthy, if the\n'
            '    cells came from a tumor, if the organism or cells were subject\n'
            '    to an experimental condition other than control, etc.'
        ),
    'date_generated':
        (
            '09. Date Generated\n'
            '    The date that these data were created. Typically, the date\n'
            '    of the publication reporting this dataset is used.\n'
            '    Must be in YYYY-MM-DD format.\n'
            '    e.g. \"2019-01-01\"'
        ),
    'count_format':
        (
            '10. Count Format\n'
            '    The format of the expression values in the actual data.\n'
            '    Refers to the normalization state of the data. There can be\n'
            '    multiple values here, delimited by \";\". The first value\n'
            '    MUST correspond to the expression matrix under lfile[\'matrix\'].\n'
            '    If there are multiple matrices under lfile[\'layers\'], the\n'
            '    HDF5 datasets should be named their corresponding \'count_format\'\n'
            '    value. Must be one of the following values:\n'
            '        {\n'
            '         \'raw\', \'cpm\', \'tpm\',\n'
            '         \'rpkm\', \'fpkm\'\n'
            '        }\n'
            '    Enter \"LIST\" to view these options.'
            '    e.g. raw;tpm;OTHER'
        ),
    'umis':
        (
            '11. UMIs\n'
            '    Whether or not the read counts were quantified using UMIs\n'
            '    or not. Must be one of these two values: {\"True\", \"False\"}.'
        ),
    'spikeins':
        (
            '12. Spike-ins\n'
            '    Whether or not spike-in transcripts were included.\n'
            '    Must be one of these two values: {\"True\", \"False\"}.'
        ),
    'technology':
        (
            '13. Technology\n'
            '    The protocol used to capture and sequence the reads.\n'
            '    Should be all lowercase and only alphanumeric.\n'
            '    Enter \"LIST\" to view existing values.\n'
            '    e.g. \"10xchromiumv2\", \"smartseq2\"'
        ),
    'doi':
        (
            '14. DOI\n'
            '    The DOI for the publication that presented these data.\n'
            '    Must be in minimal format\n'
            '    e.g. \"10.1038/s41467-019-13056-x\"'
        ),
    'accession':
        (
            '15. Accession\n'
            '    The accession number, if available, for the data. This\n'
            '    can be a GEO Series ID, an ArrayExpress ID, or something\n'
            '    else similar.\n'
            '    e.g. \"GSE108291\", \"E-MTAB-7195\"'
        ),
    'date_integrated':
        (
            '16. Date Integrated\n'
            '    The date that this dataset was entered into this database.\n'
            '    Must be in YYYY-MM-DD format.\n'
            '    e.g. \"2019-01-01\"'
        ),
    'uuid':
        (
            '17. UUID\n'
            '    The UUID associated with this entry in the database.\n'
            '    Should be generated with the python method: uuid.uuid4()\n'
            '    Must be in the string form of a UUID:\n'
            '        12345678-1234-5678-1234-567812345678\n'
            '    where each of the 32 digits is one hexadecimal digit.\n'
            '    This can be gotten with str(uuid), if uuid is a\n'
            '    python class uuid.UUID'
        ),
    'gillis_lab':
        (
            '18. Gillis Lab\n'
            '    Whether or not the dataset was scraped from a repository.\n'
            '    False if it is internally generated data.\n'
            '    Must be one of these two values: {\"True\", \"False\"}.'
        )
}
_EM_COLUMN_INDEX = {
    'title'          : 0,
    'authors'        : 1,
    'abstract'       : 2,
    'species'        : 3,
    'tissue'         : 4,
    'number_of_cells': 5,
    'author_clusters': 6,
    'condition'      : 7,
    'date_generated' : 8,
    'count_format'   : 9,
    'umis'           : 10,
    'spikeins'       : 11,
    'technology'     : 12,
    'doi'            : 13,
    'accession'      : 14,
    'date_integrated': 15,
    'uuid'           : 16,
    'gillis_lab'     : 17
}
_EM_COLUMN_MANDATORY = {
    'title'          : False,
    'authors'        : False,
    'abstract'       : False,
    'species'        : True,
    'tissue'         : False,
    'number_of_cells': True, 
    'author_clusters': True,
    'condition'      : False,
    'date_generated' : True,
    'count_format'   : True,
    'umis'           : False,
    'spikeins'       : False,
    'technology'     : False,
    'doi'            : False,
    'accession'      : False,
    'date_integrated': True,
    'uuid'           : True,
    'gillis_lab'     : True
}

_IMU_CELL_COLUMN_DESCRIPTIONS = {
    'cluster':
        (
            '01. Cluster\n'
            '    The most specific cell grouping provided by the authors.\n'
            '    These should generally correspond to cell type, but this\n'
            '    is a gray area and is up for debate.\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'species':
        (
            '02. Species\n'
            '    The species of origin for the cells. Must be the correctly\n'
            '    capitalized scientific name. If there is additional\n'
            '    information, such as a strain, it may be added after a\n'
            '    \"|\" (vertical bar) delimiter. e.g. \"Escherichia coli|K-12\"\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'tissue':
        (
            '03. Tissue\n'
            '    The tissue of origin, taken from the Tissue List. This should\n'
            '    not be a parallel of cluster. this field is intended to be a more\n'
            '    high-level organizing field, especially if multiple datasets are\n'
            '    concatenated or if a dataset took cells from multiple tissues.\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'source_organism':
        (
            '04. Source Organism\n'
            '    The specific organism of origin. For example, if cells\n'
            '    from multiple mice were pooled, the original mouse ID\n'
            '    should be stored here.\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'sex':
        (
            '05. Sex\n'
            '    The sex of the organism of origin. \"M\" for cis-canonical\n'
            '    male and \"F\" for cis-canonical female.\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'condition':
        (
            '06. Condition\n'
            '    The experimental or disease condition of the cell. If they are\n'
            '    healthy, untreated cells, this value should be \"NORMAL\".\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'batch':
        (
            '07. Batch\n'
            '    Any batch variable not covered by source_organism, sex, or\n'
            '    condition. If there are multiple batch variables remaining,\n'
            '    they should be \"|\"-concatenated to form all possible\n'
            '    combinations.\n'
            '\n'
            '    Missing values are designated with \"-1\". If a cell\'s value\n'
            '    is not missing, but is invalid, it is designated with \"OTHER\"'
        ),
    'uuid':
        (
            '08. UUID\n'
            '    The UUID of the dataset entry the cell is associated with.\n'
            '    Should be identical for all cells inside a single dataset\n'
            '    entry. This field is intended for when a user concatenates\n'
            '    cells from multiple datasets.'
        )
}

# This is mirrored in the access.R file.
# If you ever change this, MAKE SURE TO UPDATE
# THE ONE IN THE access.R FILE!
_IMU_CELL_COLUMN_INDEX = {
    'cluster'        : 0,
    'species'        : 1,
    'tissue'         : 2,
    'source_organism': 3,
    'sex'            : 4,
    'condition'      : 5,
    'batch'          : 6,
    'uuid'           : 7
}

_IMU_CELL_COLUMN_MANDATORY = {
    'cluster'        : False,
    'species'        : False,
    'tissue'         : False,
    'source_organism': False,
    'sex'            : False,
    'condition'      : False,
    'batch'          : False,
    'uuid'           : True
}

_IMU_GENE_COLUMN_DESCRIPTIONS = {
    'Accession':
        (
            '01. Accession\n'
            '    Any standardized gene ID. These are typically accession numbers\n'
            '    from a gene database, e.g. an ENSEMBL gene ID. These must be\n'
            '    unique among all genes in the dataset.'
        ),
    'Gene':
        (
            '02. Gene\n'
            '    Any human-readable, non-standardized gene ID. These are typically\n'
            '    human-readable "common names", e.g. TP53. These do not have to\n'
            '    be unique among genes in a dataset.'
        )
}
# This is mirrored in the access.R file.
# If you ever change this, MAKE SURE TO UPDATE
# THE ONE IN THE access.R FILE!
_IMU_GENE_COLUMN_INDEX = {
    'Accession' : 0,
    'Gene'      : 1
}
_IMU_GENE_COLUMN_MANDATORY = {
    'Accession' : True,
    'Gene'      : True
}

_TISSUE_LIST = {
    'Nervous':
        [
            'brain',
            'spinal cord',
            'peripheral nervous system'
        ],
    'Respiratory':
        [
            'lung',
            'trachea'
        ],
    'Digestive':
        [
            'tongue',
            'esophagus',
            'stomach',
            'liver',
            'gall bladder',
            'pancreas',
            'small intestine',
            'large intestine'
        ],
    'Cardiovascular':
        [
            'heart',
            'vasculature',
            'blood'
        ],
    'Immune':
        [
            'spleen',
            'thymus',
            'lymph node'
        ],
    'Urinary':
        [
            'urethra',
            'kidney',
            'urinary bladder'
        ],
    'Musculoskeletal':
        [
            'ligament',
            'tendon',
            'skeletal muscle',
            'smooth muscle',
            'bone'
        ],
    'Gland':
        [
            'thyroid gland'
        ],
    'Reproductive':
        [
            'uterus',
            'placenta',
            'ovary',
            'clitoris',
            'vagina',
            'testis',
            'penis'
        ],
    'Sensory':
        [
            'ear',
            'eye',
            'nose'
        ],
    'Skin':
        [
            'skin'
        ]
}