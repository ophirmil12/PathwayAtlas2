from os.path import join as pjoin

#  CLUSTER CONFIG (MACHINE DEPENDENT)
BASE_P = "/cs/labs/dina/ophirmil12/PathwayAtlas2"

#  PATHS
DATA_P = pjoin(BASE_P, 'data')                                      # The main data folder

#           CBIO
CBIO_RAW_P = pjoin(DATA_P, 'cbio_raw')                              # Raw downloaded data from cBioPortal

CBIO_MUTATION_STUDIES = pjoin(CBIO_RAW_P, 'mutation_studies')       # Raw downloaded mutations studies

CBIO_P = pjoin(DATA_P, 'cbio')                                      # Processed data from cBioPortal

CBIO_MUTATION_STUDIES_WITH_SEQUENCES = pjoin(CBIO_P, 'mutation_studies_with_sequences')  # Studies mutations with added sequences
CBIO_CANCER_MUTATIONS = pjoin(CBIO_P, 'cancer_mutations')           # Merged mutations studies by cancer


#           KEGG
KEGG_RAW_P = pjoin(DATA_P, 'kegg_raw')                              # Raw downloaded data from KEGG
KEGG_P = pjoin(DATA_P, 'kegg')                                      # Processed data from KEGG

KEGG_GENES_P = pjoin(KEGG_RAW_P, 'genes')                           # Gene objects (pickles)

KEGG_PATHWAY_OBJECTS_P = pjoin(KEGG_P, 'pathway_dicts')             # Dictionaries of a pathway`s genes to CSV file name
KEGG_PATHWAY_SCORES_P = pjoin(KEGG_P, 'pathway_snvs')               # All snvs scoring for entire pathway
KEGG_GENE_SCORES_P = pjoin(KEGG_P, 'gene_snvs')                     # All snvs scoring for single gene
KEGG_PATHWAY_METADATA_P = pjoin(KEGG_P, 'pathway_id_to_metadata.pickle')    # Mapping from pathway KEGG ID to the pathway's metadata

KEGG_PATHWAY_METADATA_P = pjoin(KEGG_P, 'pathway_id_to_metadata.pickle')    # Mapping from pathway KEGG ID to the pathway's metadata


#           ESM
ESM_EMBEDDINGS_P = pjoin(DATA_P, 'esm_1b_emb')                      # Embeddings for all sequences

#           DISORDER PREDICT
DISORDER_PRED_P = pjoin(DATA_P, 'disorder_pred')                    # Predictions of disorder of all sequences


#           CLINVAR
CLINVAR_P = pjoin(DATA_P, 'clinvar')                                # The ClinVar things folder
CLINVAR_MODELS_P = pjoin(CLINVAR_P, 'clinvar_models.pickle')    # The trained models
CLINVAR_DATA_TABLE_P = pjoin(CLINVAR_P, 'clinvar_data.csv')         # The data of ClinVar


#           RESULTS
RESULTS_P = pjoin(BASE_P, 'results')                                # The basic results (textual/csv)
RESULTS_DISTANCES_P = pjoin(RESULTS_P, 'distances')                 # The calculated bg-cancer distances


#           PLOTS
PLOTS_P = pjoin(RESULTS_P, 'plots')                                 # plots







#  ESM
ESM1B_MODEL = 'esm1b_t33_650M_UR50S'
REP_LAYERS = [33]

ESM_MAX_LENGTH = 1020
SLIDING_WINDOW_STRIDE = 300


"""
pip install fair-esm
import esm
import os
_, alphabet = esm.pretrained.load_model_and_alphabet("esm1b_t33_650M_UR50S")
print([(i, alphabet.get_tok(i)) for i in range(4, 24)])
['<cls>', '<pad>', '<eos>', '<unk>', 'L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K', 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C', 'X', 'B', 'U', 'Z', 'O', '.', '-', '<null_1>', '<mask>']
[(4, 'L'), (5, 'A'), (6, 'G'), (7, 'V'), (8, 'S'), (9, 'E'), (10, 'R'), (11, 'T'), (12, 'I'), (13, 'D'), (14, 'P'), (15, 'K'), (16, 'Q'), (17, 'N'), (18, 'F'), (19, 'Y'), (20, 'M'), (21, 'H'), (22, 'W'), (23, 'C')]
"""

ESM_1B_LOGITS_INDEX_TO_AA = {
    0: '<cls>', 1: '<pad>', 2: '<eos>', 3: '<unk>',
    4: 'L', 5: 'A', 6: 'G', 7: 'V', 8: 'S', 9: 'E', 10: 'R', 11: 'T', 12: 'I', 13: 'D',
    14: 'P', 15: 'K', 16: 'Q', 17: 'N', 18: 'F', 19: 'Y', 20: 'M', 21: 'H', 22: 'W', 23: 'C',
    24: 'X', 25: 'B', 26: 'U', 27: 'Z', 28: 'O',
    29: '.', 30: '-',
    31: '<null_1>', 32: '<mask>'
}







#  META PREDICT (SEQ -> DISORDER PREDICTION)
def V3_version_letters(sequence: str) -> str:
    return sequence.replace('B', 'N').replace('U', 'C').replace('X', 'G').replace('Z', 'Q')

DISORDERED_THRESHOLD = 0.7







#  CLINVAR













#  BIOLOGY
STOP_CODONS = ['tag', 'taa', 'tga']
CODON_TRANSLATOR = {'ata': 'I', 'atc': 'I', 'att': 'I', 'atg': 'M', 'aca': 'T',
                    'acc': 'T', 'acg': 'T', 'act': 'T', 'aac': 'N', 'aat': 'N',
                    'aaa': 'K', 'aag': 'K', 'agc': 'S', 'agt': 'S', 'aga': 'R',
                    'agg': 'R', 'cta': 'L', 'ctc': 'L', 'ctg': 'L', 'ctt': 'L',
                    'cca': 'P', 'ccc': 'P', 'ccg': 'P', 'cct': 'P', 'cac': 'H',
                    'cat': 'H', 'caa': 'Q', 'cag': 'Q', 'cga': 'R', 'cgc': 'R',
                    'cgg': 'R', 'cgt': 'R', 'gta': 'V', 'gtc': 'V', 'gtg': 'V',
                    'gtt': 'V', 'gca': 'A', 'gcc': 'A', 'gcg': 'A', 'gct': 'A',
                    'gac': 'D', 'gat': 'D', 'gaa': 'E', 'gag': 'E', 'gga': 'G',
                    'ggc': 'G', 'ggg': 'G', 'ggt': 'G', 'tca': 'S', 'tcc': 'S',
                    'tcg': 'S', 'tct': 'S', 'ttc': 'F', 'ttt': 'F', 'tta': 'L',
                    'ttg': 'L', 'tac': 'Y', 'tat': 'Y', 'taa': '_', 'tag': '_',
                    'tgc': 'C', 'tgt': 'C', 'tga': '_', 'tgg': 'W'}
STOP_AA = '_'
CODON_LENGTH = 3
NA_CHANGE = {'a': 'tcg', 't': 'acg', 'c': 'gta', 'g': 'cta'}










#  PSSM
UNIFY_PSSM = {
            "A>G": 1 / 12, "A>C": 1 / 12, "A>T": 1 / 12,
            "G>A": 1 / 12, "G>C": 1 / 12, "G>T": 1 / 12,
            "T>A": 1 / 12, "T>C": 1 / 12, "T>G": 1 / 12,
            "C>A": 1 / 12, "C>G": 1 / 12, "C>T": 1 / 12
        }

# https://academic.oup.com/nargab/article/3/3/lqab079/6371225#supplementary-data
MICHAL_HN1_PSSM = {
    "A>G": 0.000176772262074062, "A>C": 0.0000347746137285225, "A>T": 0.0000186740902005341,
    "G>A": 0.000243930190507464, "G>C": 0.0000455742821842972, "G>T": 0.0000323348824129854,
    "C>A": 0.0000410441152096065, "C>G": 0.0000510887825100302, "C>T": 0.0002600933180129,
    "T>A": 0.0000205111808477017, "T>G": 0.0000291510207312577, "T>C": 0.000195170174131106
}
# "A>A": 0.999769779033996, "G>G": 0.999678160644895, "C>C": 0.999647773784267, , "T>T": 0.99975516762429







#  CBIOPORTAL API
VERBOSE = {'critical': 0, 'program_warning': 1, 'program_progress': 1,
           'thread_warnings': 2, 'thread_progress': 3, 'raw_warnings': 3}
CBIO_BASE_URL = 'https://www.cbioportal.org/api'
CBIO_API_URL = CBIO_BASE_URL + '/v2/api-docs'
MISSENSE_MUTATION = 'Missense_Mutation'
MUTATION_STUDY_COLUMNS = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Protein', 'Variant']
DUPLICATE_EXCLUSION_COLUMNS = MUTATION_STUDY_COLUMNS + ['PatientKey']
STUDY_COLUMNS = MUTATION_STUDY_COLUMNS + ['PatientId', 'PatientKey', 'SampleId', 'StudyId', 'RefDNA']



#  UNIPROT API
Q_UID_PROT_ALL = "fields=&gene&format=tsv&query={}+AND+organism_id:9606"
UNIP_QUERY_URL = "https://rest.uniprot.org/uniprotkb/search?"
Q_UNIP_ALL_ISOFORMS = UNIP_QUERY_URL + "&format=fasta&query=" \
                                          "(accession:{}+AND+is_isoform:true)+OR+(accession:{}+AND+is_isoform:false)"
VARIATION_REGEX = r"([A-Z])(\d+)([A-Z])"



#  KEGG API
KEGG_API_URL = 'https://rest.kegg.jp'
COMMAND_TYPES = ['link', 'list', 'conv', 'get']
KEGG_LIST_COMMAND = KEGG_API_URL + '/list/{}/{}'
KEGG_LINK_COMMAND = KEGG_API_URL + '/link/{}/{}'
KEGG_CONV_COMMAND = KEGG_API_URL + '/conv/{}/{}'
KEGG_GET_COMMAND = KEGG_API_URL + '/get/' + '{}/{}'

KEGG_GENE_SEQ_URL = f'{KEGG_API_URL}/get/' + '{}/{}'
KEGG_HOMO_SAPIENS = 'hsa'
KEGG_PATHWAY_PREFIX = 'hsa'
KEGG_MODULE_PREFIX = 'M'

KEGG_MAX_IDS = 10
NETWORK_TYPES = ('module', 'pathway')
EMPTY_SET = {''}
EMPTY_LIST = ['']
GENE_SEPERATOR = '///\n'

GENE_DATA = {'kegg_id': None, 'uniprot_id': None, 'aa_seq': '', 'na_seq': '',
             'chr': None, 'start': None, 'end': None, 'coding_type': None, 'ref_names': None}

KEG_POSITION_RE = "(\d+)\.{2}(\d+)"

KEGG_COL = 'KeggId'





# QUERIES
UIDS_COL_IDX = 0
REVIEWED_COL_IDX = 2
GENE_NAME_COL_IDX = 4
UNIP_REVIEWED = 'reviewed'
UNIP_QUERY_URL = "https://rest.uniprot.org/uniprotkb/search?"

#  REQUESTS AND OS CONSTANTS
TIMEOUT = 20.0
WAIT_TIME = 1.0
RETRIES = 7
RETRY_STATUS_LIST = [429, 500, 502, 503, 504, 403, 400]
DEFAULT_HEADER = "https://"
WORKERS = None  # use default amount of CPUs
KEGG_API_RECOMMENDED_WORKERS = 6

#  ERRORS
NETWORK_TYPE_ERROR = f'Network type must be one of: {", ".join(NETWORK_TYPES)}'
NETWORK_ID_ERROR = f'KEGG id must be of a KEGG module or KEGG pathway'
LOAD_OBJ_ERROR = 'Data missing or invalid for {}. ' \
                 '\nDelete instance from DB and recreate the object'
CON_ERR_FUS = "Connection Error in fetch_uniport_sequences while fetching isoforms for {}\nURL: "
CON_ERR_GENERAL = "Connection Error in {} on protein {}"
CON_ERR_UFN = "Connection Error in uid_from_name failed to fetch Uniprot IDs for protein {}"





#  PROTEIN COVERAGE
COVERAGE_PERCENTAGE_THRESHOLD = 1  # 1%
ABSOLUTE_COUNT_THRESHOLD = 10  # 10 mutations





#  PAPER COLORING
import matplotlib.pyplot as plt
from cycler import cycler

MY_PALETTE = [
    '#CB7673',
    '#447D68',
    '#EC9D58',
    '#476067',
    '#F3B8BA',
    '#67383E',
    '#A5AE77',
    '#9F403A',
    '#787A91',
    '#B7CADB',
    '#917FB3'
]
def set_paper_palette():
    plt.rcParams['axes.prop_cycle'] = cycler(color=MY_PALETTE)




