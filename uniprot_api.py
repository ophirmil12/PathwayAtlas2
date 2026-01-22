# The UniProt API class

import warnings
from Bio import SeqIO
import requests
from requests.adapters import HTTPAdapter, Retry
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool as Pool
import re
import tempfile
import Retry

Q_UID_PROT_ALL = "fields=&gene&format=tsv&query={}+AND+organism_id:9606"
UNIP_QUERY_URL = "https://rest.uniprot.org/uniprotkb/search?"
Q_UNIP_ALL_ISOFORMS = UNIP_QUERY_URL + "&format=fasta&query=" \
                                          "(accession:{}+AND+is_isoform:true)+OR+(accession:{}+AND+is_isoform:false)"
VARIATION_REGEX = r"([A-Z])(\d+)([A-Z])"

# QUERIES
UIDS_COL_IDX = 0
REVIEWED_COL_IDX = 2
GENE_NAME_COL_IDX = 4
UNIP_REVIEWED = 'reviewed'
UNIP_QUERY_URL = "https://rest.uniprot.org/uniprotkb/search?"
VERBOSE = {'critical': 0, 'program_warning': 1, 'program_progress': 1,
           'thread_warnings': 2, 'thread_progress': 3, 'raw_warnings': 3}

#  REQUESTS AND OS CONSTANTS
TIMEOUT = 20.0
WAIT_TIME = 1.0
RETRIES = 7
RETRY_STATUS_LIST = [429, 500, 502, 503, 504, 403, 400]
DEFAULT_HEADER = "https://"
WORKERS = None  # use default amount of CPUs
KEGG_API_RECOMMENDED_WORKERS = 6

#  ERRORS
NETWORK_ID_ERROR = f'KEGG id must be of a KEGG module or KEGG pathway'
LOAD_OBJ_ERROR = 'Data missing or invalid for {}. ' \
                 '\nDelete instance from DB and recreate the object'
CON_ERR_FUS = "Connection Error in fetch_uniport_sequences while fetching isoforms for {}\nURL: "
CON_ERR_GENERAL = "Connection Error in {} on protein {}"
CON_ERR_UFN = "Connection Error in uid_from_name failed to fetch Uniprot IDs for protein {}"


def print_if(verbose: object, thr: object, text: object) -> object:
    """
    print text if verbose > thr
    :param verbose: int
    :param thr: int
    :param text: str
    :return:
    """
    if verbose >= thr:
        print(text)

def process_fastas(text):
    """
    process multiple fasta sequences
    :return: {id: sequence}
    """
    temp = tempfile.TemporaryFile(mode='w+t')
    temp.writelines(text)
    temp.seek(0)
    ret = {seq_record.id: str(seq_record.seq) for seq_record in SeqIO.parse(temp, "fasta")}
    temp.close()
    return ret

def read_in_chunks(array, chunk_size):
    for i in range(0, len(array), chunk_size):
        yield array[i:i + chunk_size]


def create_session(header, retries=5, wait_time=0.5, status_forcelist=None):
    """
    Creates a session using pagination
    :param header: str url header session eill apply to
    :param retries: int number of retries on failure
    :param wait_time: float time (sec) between attempts
    :param status_forcelist: list HTTP status codes that we should force a retry on
    :return: requests session
    """
    s = requests.Session()
    retries = Retry(total=retries,
                    backoff_factor=wait_time,
                    status_forcelist=status_forcelist)

    s.mount(header, HTTPAdapter(max_retries=retries))
    return s


def extract_accession(uid):
    if "|" in uid:
        return uid.split("|")[1]  # get the middle part
    return uid


def safe_get_request(session, url, timeout=TIMEOUT, warning_msg='connection failed', return_on_failure=None):
    """
    creates a user friendly request raises warning on ConnectionError but will not crush
    verbose_level = 3 will return raw Error massage in warning
    :param session: requests session obj
    :param url: str url to query
    :param timeout: float max time to wait for response
    :param warning_msg: str msg to display on failure
    :param return_on_failure: value to return upon exception
    :return: response
    """
    try:
        r = session.get(url, timeout=timeout)
    except requests.exceptions.ConnectionError as e:
        warnings.warn(warning_msg)
        return return_on_failure
    return r


def multiprocess_task(tasks, target, workers=None, callback=lambda x: x):
    """
    :param tasks: list of iterables
    :param target: callable
    :param workers: optional int number of CPUs that will be used otherwise maximum available
    :param callback: callable
    :return:
    """
    workers = workers if workers else cpu_count()
    with Pool(workers) as p:
        try:
            for status in p.starmap(target, tasks):
                callback(status)
        except Exception as e:
            print(f"Multiprocessing task failed: {e}")


def warn_if(verbose, thr, text):
    """
    print text if verbose > thr
    :param verbose: int
    :param thr: int
    :param text: str
    :return:
    """
    if verbose >= thr:
        warnings.warn(text)


def safe_post_request(session, url, timeout, verbose_level, warning_msg='connection failed', return_on_failure=None,
                      warning_thr=VERBOSE['thread_warnings'], raw_err_thr=VERBOSE['raw_warnings']):
    """
    creates a user friendly request raises warning on ConnectionError but will not crush
    verbose_level = 3 will return raw Error massage in warning
    :param session: requests session obj
    :param url: str url to query
    :param timeout: float max time to wait for response
    :param verbose_level: int
    :param warning_msg: str msg to display on failure
    :param return_on_failure: value to return upon exception
    :param raw_err_thr: int threshold to print raw error messages
    :param warning_thr: int threshold to print warning messages
    :return: response
    """
    try:
        r = session.post(url, timeout=timeout)
    except requests.exceptions.ConnectionError as e:
        warn_if(verbose_level, warning_thr, warning_msg)
        warn_if(verbose_level, raw_err_thr, f"{e}")
        return return_on_failure
    return r



class UniprotApi:
    """
    This class is responsible to connect to online DBs and retrieve information
    """

    def __init__(self, verbose_level=1):
        self._v = verbose_level

    def fetch_uniport_sequences(self, uid):
        """
        Retrieve all known isoforms from uniprot
        :param uid: Uniprot id only primary name
        :return: {uid_iso_index: sequence}
        """
        print_if(self._v, VERBOSE['thread_progress'], "Retrieving isoforms from Uniprot...")
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        url = Q_UNIP_ALL_ISOFORMS.format(uid, uid)
        response = safe_post_request(s, url, TIMEOUT, self._v, CON_ERR_FUS.format(uid, url))
        if not response.ok:
            print_if(self._v, VERBOSE['thread_warnings'], CON_ERR_GENERAL.format('fetch_uniprot_sequences', uid))
            return {}
        if response.text == '':
            print_if(self._v, VERBOSE['thread_warnings'], f"no sequences found for {uid}")
            return {}
        return process_fastas(response.text)

    def expand_isoforms(self, ref_name, ref_mut=None, reviewed=True):
        """
        expand protein isoforms using all relevant Uniprot accession
        this will not override the default protein isoforms
        :param ref_name: protein name
        :param ref_mut: Mutation obj if given will search for isoform with the given mutation
        :return: {uid_iso_index: seq}
        """
        uids = self.uid_from_name(ref_name)['reviewed'] if reviewed else self.uid_from_name(ref_name)['all_enteries']
        isoforms = {}

        if isinstance(ref_mut, str):
            match = re.match(VARIATION_REGEX, ref_mut)
            if not match:
                raise ValueError(f"Invalid mutation format: {ref_mut}")
            orig_aa, loc, mut_aa = match.groups()
            idx = int(loc) - 1
        else:
            idx = ref_mut.loc - 1
            orig_aa = ref_mut.origAA


        # Search for isoforms that have a fitting sequence to the mutation
        for uid in uids:
            res = self.fetch_uniport_sequences(uid)
            if ref_mut:
                for iso, seq in res.items():
                    if idx >= len(seq):
                        continue
                    if seq[idx] == orig_aa:
                        return {iso: seq}
            isoforms = {**isoforms, **res}
        return isoforms if not ref_mut else {}

    def uid_from_name(self, ref_name):
        """
        return uniprot-id given a protein ref_name
        :param ref_name: protein name
        :return: {'reviewed': [...], 'non_reviewed': [...]}
        """
        ret = {'reviewed': [], 'non_reviewed': [], 'main_entery': [], 'all_enteries': [], 'aliases': []}
        query = UNIP_QUERY_URL + Q_UID_PROT_ALL.format(ref_name)
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        r = safe_get_request(s, query, TIMEOUT, self._v, CON_ERR_UFN.format(ref_name))
        if not r:
            return ret
        if r.text == '':
            return ret
        ret = self._process_uid_query(r.text)
        return ret

    @staticmethod
    def _process_uid_query(data):
        ret = {'reviewed': [], 'non_reviewed': [], 'main_entery': [], 'all_enteries': [], 'aliases': []}
        rows = data.split('\n')[1:-1]  # first is header last is blank
        if not rows:
            return ret
        main_entry = rows[0].split('\t')[UIDS_COL_IDX]  # first entry is considered main
        ret['main_entery'].append(main_entry)
        for row in rows:
            values = row.split('\t')
            entry, reviewed, gene = values[UIDS_COL_IDX], values[REVIEWED_COL_IDX], values[GENE_NAME_COL_IDX].split(" ")
            ret['all_enteries'].append(entry)
            ret['aliases'] += gene
            if reviewed == UNIP_REVIEWED:
                ret['reviewed'].append(entry)
            else:
                ret['non_reviewed'].append(entry)
        return ret