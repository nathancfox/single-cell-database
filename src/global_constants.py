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

_PATH_TO_DATABASE:
    The full, absolute path to the root folder of the database.
_PATH_TO_METADATA:
    The full, absolute path to the external metadata file.
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
# These 2 are mirrored in the access.R file.
# If you ever change these, MAKE SURE TO UPDATE
# THE ONES IN THE access.R FILE!
_PATH_TO_DATABASE = '/data/single_cell_database'
_PATH_TO_METADATA = '/data/single_cell_database/external_metadata.tsv'

_EM_COLUMN_DESCRIPTIONS = {
    'species': 
        (
            '01. Species\n'
            '    The species that the cells originated from. Should be a\n'
            '    two word scientific name in standard format. If there is\n'
            '    additional information, such as a string, it may be added\n'
            '    after a \"|\" delimiter. You can have multiple species/strains\n'
            '    in a dataset. These should be \";\" delimited.'
            '    e.g. \"Homo sapiens\", \"Mus musculus;Escherichia coli|K-12\"'
        ),
    'tissue':
        (
            '02. Tissue\n'
            '    The tissue/organ that the cells originated from. Must be from\n'
            '    the TISSUE_LIST. Enter \"LIST\" to view the options. Multiple\n'
            '    tissues in a dataset can be \";\" delimited.\n'
            '    e.g. \"brain\", \"kidney\", \"blood\"'
        ),
    'number_of_cells':
        (
            '03. Number of Cells\n'
            '    The number of cells in this dataset. Should just be a\n'
            '    simple, non-negative integer.\n'
            '    e.g. \"5402\"'
        ),
    'condition':
        (
            '04. Condition\n'
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
            '05. Date Generated\n'
            '    The date that these data were created. Typically, the date\n'
            '    of the publication reporting this dataset is used.\n'
            '    Must be in YYYY-MM-DD format.\n'
            '    e.g. \"2019-01-01\"'
        ),
    'count_format':
        (
            '06. Count Format\n'
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
            '07. UMIs\n'
            '    Whether or not the read counts were quantified using UMIs\n'
            '    or not. Must be one of these two values: {\"True\", \"False\"}.'
        ),
    'spikeins':
        (
            '08. Spike-ins\n'
            '    Whether or not spike-in transcripts were included.\n'
            '    Must be one of these two values: {\"True\", \"False\"}.'
        ),
    'technology':
        (
            '09. Technology\n'
            '    The protocol used to capture and sequence the reads.\n'
            '    Should be all lowercase and only alphanumeric.\n'
            '    Enter \"LIST\" to view existing values.\n'
            '    e.g. \"10xchromiumv2\", \"smartseq2\"'
        ),
    'doi':
        (
            '10. DOI\n'
            '    The DOI for the publication that presented these data.\n'
            '    Must be in minimal format\n'
            '    e.g. \"10.1038/s41467-019-13056-x\"'
        ),
    'accession':
        (
            '11. Accession\n'
            '    The accession number, if available, for the data. This\n'
            '    can be a GEO Series ID, an ArrayExpress ID, or something\n'
            '    else similar.\n'
            '    e.g. \"GSE108291\", \"E-MTAB-7195\"'
        ),
    'date_integrated':
        (
            '12. Date Integrated\n'
            '    The date that this dataset was entered into this database.\n'
            '    Must be in YYYY-MM-DD format.\n'
            '    e.g. \"2019-01-01\"'
        ),
    'uuid':
        (
            '13. UUID\n'
            '    The UUID associated with this entry in the database.\n'
            '    Should be generated with the python method: uuid.uuid4()\n'
            '    Must be in the string form of a UUID:\n'
            '        12345678-1234-5678-1234-567812345678\n'
            '    where each of the 32 digits is one hexadecimal digit.\n'
            '    This can be gotten with str(uuid), if uuid is a\n'
            '    python class uuid.UUID'
        ),
    'file_location':
        (
            '14. File location\n'
            '    A full, absolute path to the folder holding the loom\n'
            '    file and all internal metadata for this dataset.\n'
            '    Must end with a \'/\'.\n'
            '    e.g. \"/data/single_cell_database/12345678-1234-5678-1234-567812345678/\"'
        ),
    'internal':
        (
            '15. Internal\n'
            '    Whether or not the dataset was scraped from a repository.\n'
            '    False if it is internally generated data.\n'
            '    Must be one of these two values: {\"True\", \"False\"}.'
        )
}
_EM_COLUMN_INDEX = {
    'species'        : 0,
    'tissue'         : 1,
    'number_of_cells': 2,
    'condition'      : 3,
    'date_generated' : 4,
    'count_format'   : 5,
    'umis'           : 6,
    'spikeins'       : 7,
    'technology'     : 8,
    'doi'            : 9,
    'accession'      : 10,
    'date_integrated': 11,
    'uuid'           : 12,
    'file_location'  : 13,
    'internal'       : 14
}
_EM_COLUMN_MANDATORY = {
    'species'        : True,
    'tissue'         : False,
    'number_of_cells': True, 
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
    'file_location'  : True,
    'internal'       : True
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

_IMU_GENE_COLUMN_DESCRIPTIONS = {}
_IMU_GENE_COLUMN_INDEX = {}
_IMU_GENE_COLUMN_MANDATORY = {}

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