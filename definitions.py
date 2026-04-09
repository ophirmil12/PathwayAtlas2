from os.path import join as pjoin

#  CLUSTER CONFIG (MACHINE DEPENDENT)
# TODO: make a .env file with the BASE_P, and use the dotenv library (needed also for the R script)
BASE_P = "/cs/labs/dina/ophirmil12/PathwayAtlas2"

#  PATHS
DATA_P = pjoin(BASE_P, 'data')                                      # The main data folder

#           CBIO
CBIO_RAW_P = pjoin(DATA_P, 'cbio_raw')                              # Raw downloaded data from cBioPortal

CBIO_MUTATION_STUDIES = pjoin(CBIO_RAW_P, 'mutation_studies')       # Raw downloaded mutations studies
CBIO_PATIENT_CLINICAL_STUDIES_P = pjoin(CBIO_RAW_P, 'patients_clinical_data')

CBIO_P = pjoin(DATA_P, 'cbio')                                      # Processed data from cBioPortal

CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P = pjoin(CBIO_P, 'mutation_studies_with_sequences')  # Studies mutations with added sequences
CBIO_CANCER_MUTATIONS_P = pjoin(CBIO_P, 'cancer_mutations')           # Merged mutations studies by cancer
CBIO_CANCER_MUTATIONS_NO_DUPS_P = CBIO_CANCER_MUTATIONS_P + "_no_dups"


#           KEGG
KEGG_RAW_P = pjoin(DATA_P, 'kegg_raw')                              # Raw downloaded data from KEGG
KEGG_P = pjoin(DATA_P, 'kegg')                                      # Processed data from KEGG

KEGG_GENES_P = pjoin(KEGG_RAW_P, 'genes')                           # Gene objects (pickles)

KEGG_PATHWAY_SCORES_P = pjoin(KEGG_P, 'pathway_snvs')               # All snvs scoring for entire pathway
KEGG_GENE_SCORES_P = pjoin(KEGG_P, 'gene_snvs')                     # All snvs scoring for single gene
KEGG_PATHWAY_METADATA_FILE = pjoin(KEGG_P, 'pathway_id_to_metadata.pickle')    # Mapping from pathway KEGG ID to the pathway's metadata

KEGG_PATHWAY_CLUSTERING_P = pjoin(KEGG_P, 'pathway_clustering')     # Clustering of pathways based on their gene content

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
RESULTS_DISTANCES_NO_DUPS_P = RESULTS_DISTANCES_P + "_no_dups"
CANCER_PATIENT_SURVIVAL_P = pjoin(RESULTS_P, 'cancer_patient_survival')     # The patient survival data for each cancer


#           PLOTS
PLOTS_P = pjoin(BASE_P, 'plots')                                 # plots
KAPLAN_MEIER_P = pjoin(PLOTS_P, 'kaplan_meier')                 # Kaplan-Meier plots


#           FIGURES
FIGURES_P = pjoin(BASE_P, 'figures')





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



CANCER_FULLNAME = {
    "acc": "Adrenocortical Carcinoma",
    "aml": "Acute Myeloid Leukemia",
    "blca": "Bladder Urothelial Carcinoma",
    "brca": "Breast Invasive Carcinoma",
    "ccrcc": "Clear Cell Renal Cell Carcinoma",
    "cesc": "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
    "chol": "Cholangiocarcinoma",
    "chrcc": "Chromophobe Renal Cell Carcinoma",
    "coad": "Colon Adenocarcinoma",
    "coadread": "Colorectal Adenocarcinoma",
    "difg": "Diffuse Glioma",
    "dlbclnos": "Diffuse Large B-Cell Lymphoma, Not Otherwise Specified",
    "egc": "Esophagogastric Cancer",
    "esca": "Esophageal Carcinoma",
    "hcc": "Hepatocellular Carcinoma",
    "hgsoc": "High-Grade Serous Ovarian Carcinoma",
    "hnsc": "Head and Neck Squamous Cell Carcinoma",
    "luad": "Lung Adenocarcinoma",
    "lusc": "Lung Squamous Cell Carcinoma",
    "mixed": "Mixed Cancer Types",
    "mnet": "Metastatic Neuroendocrine Tumor",
    "nsclc": "Non-Small Cell Lung Cancer",
    "nsgct": "Non-Seminomatous Germ Cell Tumor",
    "paad": "Pancreatic Adenocarcinoma",
    "pan_cancer": "Pan-Cancer",
    "plmeso": "Pleural Mesothelioma",
    "prad": "Prostate Adenocarcinoma",
    "prcc": "Papillary Renal Cell Carcinoma",
    "read": "Rectum Adenocarcinoma",
    "skcm": "Skin Cutaneous Melanoma",
    "soft_tissue": "Soft Tissue Sarcoma",
    "stad": "Stomach Adenocarcinoma",
    "testis": "Testicular Germ Cell Tumors",
    "thpa": "Thyroid Papillary Carcinoma",
    "thym": "Thymoma",
    "ucec": "Uterine Corpus Endometrial Carcinoma",
    "ucs": "Uterine Carcinosarcoma",
    "um": "Uveal Melanoma"
}

CANCER_STUDIES_DICT = {'chol': ['chol_tcga', 'chol_tcga_gdc', 'chol_tcga_pan_can_atlas_2018'],
                  'skcm': ['skcm_tcga_gdc', 'skcm_tcga_pub_2015', 'skcm_tcga', 'skcm_tcga_pan_can_atlas_2018'],
                  'acc': ['acc_tcga_gdc', 'acc_tcga', 'acc_tcga_pan_can_atlas_2018'],
                  'coadread': ['coadread_tcga', 'coadread_tcga_pub', 'coadread_tcga_pan_can_atlas_2018'],
                  'lusc': ['lusc_tcga', 'lusc_tcga_gdc', 'lusc_tcga_pub', 'lusc_tcga_pan_can_atlas_2018'],
                  'ccrcc': ['kirc_tcga', 'kirc_tcga_pub', 'kirc_tcga_pan_can_atlas_2018', 'ccrcc_tcga_gdc'],
                  'read': ['read_tcga_gdc'], 'testis': ['tgct_tcga_pan_can_atlas_2018'],
                  'aml': ['aml_tcga_gdc', 'laml_tcga', 'laml_tcga_pub', 'laml_tcga_pan_can_atlas_2018'],
                  'mixed': ['mixed_msk_tcga_2021'],
                  'luad': ['luad_tcga', 'luad_tcga_gdc', 'luad_tcga_pub', 'luad_tcga_pan_can_atlas_2018'],
                  'coad': ['coad_tcga_gdc'], 'nsclc': ['nsclc_tcga_broad_2016'],
                       'nsgct': ['nsgct_tcga_gdc', 'tgct_tcga'],
                       'chrcc': ['chrcc_tcga_gdc', 'kich_tcga_pub', 'kich_tcga', 'kich_tcga_pan_can_atlas_2018'],
                       'prad': ['prad_tcga_pub', 'prad_tcga', 'prad_tcga_gdc', 'prad_tcga_pan_can_atlas_2018'],
                       'thpa': ['thpa_tcga_gdc', 'thca_tcga_pub', 'thca_tcga', 'thca_tcga_pan_can_atlas_2018'],
                       'ucs': ['ucs_tcga', 'ucs_tcga_pan_can_atlas_2018', 'ucs_tcga_gdc'],
                       'egc': ['stes_tcga_pub'],
                       'soft_tissue': ['sarc_tcga_pub', 'pcpg_tcga_pub', 'pcpg_tcga_pan_can_atlas_2018',
                                  'sarc_tcga', 'sarc_tcga_pan_can_atlas_2018', 'soft_tissue_tcga_gdc'],
                       'blca': ['blca_msk_tcga_2020', 'blca_tcga_pub_2017', 'blca_tcga', 'blca_tcga_gdc', 'blca_tcga_pub',
                           'blca_tcga_pan_can_atlas_2018'],
                       'plmeso': ['meso_tcga', 'meso_tcga_pan_can_atlas_2018', 'plmeso_tcga_gdc'],
                       'paad': ['paad_tcga', 'paad_tcga_gdc', 'paad_tcga_pan_can_atlas_2018'],
                       'um': ['um_tcga_gdc', 'uvm_tcga', 'uvm_tcga_pan_can_atlas_2018'],
                       'difg': ['lgg_tcga', 'lgg_tcga_pan_can_atlas_2018', 'difg_tcga_gdc', 'gbm_tcga_pub2013',
                           'gbm_tcga_pub', 'gbm_tcga_gdc', 'gbm_tcga', 'gbm_tcga_pan_can_atlas_2018', 'lgggbm_tcga_pub'],
                       'esca': ['esca_tcga_gdc', 'esca_tcga_pan_can_atlas_2018', 'esca_tcga'],
                       'thym': ['thym_tcga_gdc', 'thym_tcga', 'thym_tcga_pan_can_atlas_2018'],
                       'cesc': ['cesc_tcga_gdc', 'cesc_tcga_pan_can_atlas_2018', 'cesc_tcga'],
                       'brca': ['brca_tcga_pub2015', 'brca_tcga', 'brca_tcga_pub', 'brca_tcga_pan_can_atlas_2018',
                           'brca_tcga_gdc'],
                       'hnsc': ['hnsc_tcga', 'hnsc_tcga_gdc', 'hnsc_tcga_pub', 'hnsc_tcga_pan_can_atlas_2018'],
                       'dlbclnos': ['dlbc_tcga_pan_can_atlas_2018', 'dlbclnos_tcga_gdc', 'dlbc_tcga'],
                       'ucec': ['ucec_tcga_gdc', 'ucec_tcga', 'ucec_tcga_pub', 'ucec_tcga_pan_can_atlas_2018'],
                       'hgsoc': ['hgsoc_tcga_gdc', 'ov_tcga', 'ov_tcga_pub', 'ov_tcga_pan_can_atlas_2018'],
                       'stad': ['stad_tcga_gdc', 'stad_tcga', 'stad_tcga_pub', 'stad_tcga_pan_can_atlas_2018'],
                       'prcc': ['kirp_tcga', 'kirp_tcga_pan_can_atlas_2018', 'prcc_tcga_gdc'],
                       'mnet': ['mnet_tcga_gdc', 'pcpg_tcga'],
                       'hcc': ['hcc_tcga_gdc', 'lihc_tcga', 'lihc_tcga_pan_can_atlas_2018']}

CANCER_SUBTYPES = {"renal": ["ccrcc", "chrcc", "prcc"],
                    "lung": ["luad", "lusc", "nsclc"],
                    "colorectal": ["coad", "coadread", "read"],
                    "glioma": ["difg", "dlbclnos"],
                    "liver": ["hcc", "chol"],
                    "uterine": ["ucec", "ucs"],
                    "testicular": ["testis", "nsgct"]
                    }
CANCER_FULLNAME = {
    "acc": "Adrenocortical Carcinoma",
    "aml": "Acute Myeloid Leukemia",
    "blca": "Bladder Urothelial Carcinoma",
    "brca": "Breast Invasive Carcinoma",
    "ccrcc": "Clear Cell Renal Cell Carcinoma",
    "cesc": "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
    "chol": "Cholangiocarcinoma",
    "chrcc": "Chromophobe Renal Cell Carcinoma",
    "coad": "Colon Adenocarcinoma",
    "coadread": "Colorectal Adenocarcinoma",
    "difg": "Diffuse Glioma",
    "dlbclnos": "Diffuse Large B-Cell Lymphoma, Not Otherwise Specified",
    "egc": "Esophagogastric Cancer",
    "esca": "Esophageal Carcinoma",
    "hcc": "Hepatocellular Carcinoma",
    "hgsoc": "High-Grade Serous Ovarian Carcinoma",
    "hnsc": "Head and Neck Squamous Cell Carcinoma",
    "luad": "Lung Adenocarcinoma",
    "lusc": "Lung Squamous Cell Carcinoma",
    "mixed": "Mixed Cancer Types",
    "mnet": "Metastatic Neuroendocrine Tumor",
    "nsclc": "Non-Small Cell Lung Cancer",
    "nsgct": "Non-Seminomatous Germ Cell Tumor",
    "paad": "Pancreatic Adenocarcinoma",
    "pan_cancer": "Pan-Cancer",
    "plmeso": "Pleural Mesothelioma",
    "prad": "Prostate Adenocarcinoma",
    "prcc": "Papillary Renal Cell Carcinoma",
    "read": "Rectum Adenocarcinoma",
    "skcm": "Skin Cutaneous Melanoma",
    "soft_tissue": "Soft Tissue Sarcoma",
    "stad": "Stomach Adenocarcinoma",
    "testis": "Testicular Germ Cell Tumors",
    "thpa": "Thyroid Papillary Carcinoma",
    "thym": "Thymoma",
    "ucec": "Uterine Corpus Endometrial Carcinoma",
    "ucs": "Uterine Carcinosarcoma",
    "um": "Uveal Melanoma"
}





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

NUMBER_OF_BINS = 100





#  CBIOPORTAL API
VERBOSE = {'critical': 0, 'program_warning': 1, 'program_progress': 1,
           'thread_warnings': 2, 'thread_progress': 3, 'raw_warnings': 3}
CBIO_BASE_URL = 'https://www.cbioportal.org/api'
CBIO_API_URL = CBIO_BASE_URL + '/v2/api-docs'
MISSENSE_MUTATION = 'Missense_Mutation'
#MUTATION_STUDY_COLUMNS = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Protein', 'Variant']
MUTATION_STUDY_COLUMNS = ['Chr', 'Ref', 'Alt', 'Protein', 'Variant']

#DUPLICATE_EXCLUSION_COLUMNS = MUTATION_STUDY_COLUMNS + ['PatientKey']
DUPLICATE_EXCLUSION_COLUMNS = MUTATION_STUDY_COLUMNS + ['PatientId']

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
COVERAGE_PERCENTAGE_THRESHOLD = 1    # 1% of the gene is covered
ABSOLUTE_COUNT_THRESHOLD = 10        # 10 mutations recorded for the gene
PATHWAY_COVERAGE_THRESHOLD = 40      # 40% of the genes in the pathway are covered (see p6A)

PATHWAY_JACCARD_SIMILARITY_THRESHOLD = 0.25      # If two pathways have Jaccard similarity of PATHWAY_JACCARD_SIMILARITY_THRESHOLD or higher, we consider them similar




#  PAPER COLORING
import matplotlib.pyplot as plt
from cycler import cycler

COLOR_MAP = {
    'pathogenic': "#CB7673",      # Mauve
    'benign': "#447D68",          # Green
    'non-significant': "#EC9D58",  # Orange
    'dark-blue': "#5b7d87",
    'pink': "#F3B8BA",
    'dark-red': "#67383E",
    'light-green': "#A5AE77",
    'red': "#9F403A",
    'grey': "#787A91",
    'light-blue': "#B7CADB",
    'significant': "#917FB3"  # Purple
}
def set_paper_palette():
    plt.rcParams['axes.prop_cycle'] = cycler(color=list(COLOR_MAP.values()))
