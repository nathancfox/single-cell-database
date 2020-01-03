# These 2 are mirrored in the access.R file.
# If you ever change these, MAKE SURE TO UPDATE
# THE ONES IN THE access.R FILE!
_PATH_TO_DATABASE = '/data/single_cell_database'
_PATH_TO_METADATA = '/data/single_cell_database/external_metadata.tsv'
_COLUMN_DESCRIPTIONS = {
    'species': 
        (
            '01. Species\n'
            '    The species that the cells originated from. Should be a\n'
            '    two word scientific name in standard format. e.g. \"Homo sapiens\"'
        ),
    'organ':
        (
            '02. Organ\n'
            '    The organ that the cells originated from. Should be all\n'
            '    lowercase and as simple and high-level as possible.\n'
            '    e.g. \"brain\", \"kidney\", \"blood\"'
        ),
    'number_of_cells':
        (
            '03. Number of Cells\n'
            '    The number of cells in this dataset. Should just be a\n'
            '    simple integer. e.g. \"5402\"'
        ),
    'condition':
        (
            '04. Condition\n'
            '    Whether or not these cells fall under a disease condition\n'
            '    or an experimental condition. Must be one of these two\n'
            '    values: {\"True\", \"False\"}. Should be True if these are\n'
            '    assumed to be normal cells from a healthy organisms.\n'
            '    Should be False otherwise. This could be False if the\n'
            '    organism is unhealthy, if the cells came from a tumor,\n'
            '    if the organism or cells were subject to an experimental\n'
            '    condition other than control.'
        ),
    'date_generated':
        (
            '05. Date Generated\n'
            '    The date that these data were created. Often, the date\n'
            '    of the publication reporting this dataset is used.\n'
            '    Must be in YYYY-MM-DD format. e.g. \"2019-01-01\"'
        ),
    'count_format':
        (
            '06. Count Format\n'
            '    The format of the expression values in the actual data.\n'
            '    Are the read counts raw? Normalized? TPM? Should be\n'
            '    one of the following:\n'
            '        {\n'
            '         \'raw\', \'cpm\', \'tpm\', \'rpkm\',\n'
            '         \'fpkm\', \'other: DESCRIPTION\'\n'
            '        }\n'
            '    where DESCRIPTION can be anything that adequately describes\n'
            '    the format of the read counts. If multiple versions are\n'
            '    available, the one stored in \'matrix\' should be first.'
            '    They should be comma-delimited.'
            '    e.g. raw_tpm_other: normalized to 500,000'
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
            '    e.g. \"10xchromiumv2\", \"smartseq2\", \"celseq\"'
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
            '    The date that these data were entered into this database.\n'
            '    Must be in YYYY-MM-DD format. e.g. \"2019-01-01\"'
        ),
    'uuid':
        (
            '13. UUID\n'
            '    The UUID associated with this entry in the database.\n'
            '    Should be generated with the python method: uuid.uuid4()\n'
            '    Must be in the string form of a UUID:\n'
            '        12345678-1234-5678-1234-567812345678\n'
            '    where each of the 32 digits is one hexadecimal digit.\n'
            '    This can be gotten with str(uuid), if uuid is a \n'
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
_COLUMN_INDEX = {
    'species'        : 0,
    'organ'          : 1,
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
_COLUMN_MANDATORY = {
    'species'        : True,
    'organ'          : False,
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