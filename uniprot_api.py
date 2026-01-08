# The UniProt API class



import re

from utils import *



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