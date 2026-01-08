# The KEGG API class, KEGG Network, and KEGG Gene



from tqdm import tqdm
from typing import List
import pandas as pd
import requests
import copy
import re
import os


from definitions import *
from utils import *



class KeggApi:
    """api for kegg"""

    def __init__(self):
        """
        constructor for Kegg
        """
        self.api = create_session(KEGG_API_URL, retries=RETRIES, wait_time=WAIT_TIME,
                                  status_forcelist=RETRY_STATUS_LIST)

    @staticmethod
    def format_multiple_genes(genes):
        return '+'.join(genes)

    def _kegg_command(self, command, *params, verbose=True):
        """
        Kegg general command
        :param database: pathway | brite | module | ko | <org> | vg | vp | ag | genome | compound |
             glycan | reaction | rclass | enzyme | network | variant | disease |
             drug | dgroup | disease_ja | drug_ja | dgroup_ja | compound_ja
        :param species: str
        :param verbose: bool if true request will be printed
        :return:
        """
        query = command.format(*params)
        if verbose:
            #print(command.format(*params))
            pass
        data = safe_get_request(self.api, query)
        if not data:
            raise ConnectionError(f'querry: {query} --> Failed')
        if not data.ok:
            raise ConnectionError(f'querry: {query} --> Failed with code {data.status_code}')
        return data.text

    def kegg_command(self, command_type, *params):
        """
        Factory for kegg commands
        :param command_type: one of list | link | conv | get
        :param params: parameters for query
        :return:
        """
        assert command_type in COMMAND_TYPES, 'command_type must be one of link | list | conv | get'
        if command_type == 'list':
            return self._kegg_command(KEGG_LIST_COMMAND, *params)
        elif command_type == 'link':
            return self._kegg_command(KEGG_LINK_COMMAND, *params)
        elif command_type == 'conv':
            return self._kegg_command(KEGG_CONV_COMMAND, *params)
        elif command_type == 'get':
            return self._kegg_command(KEGG_GET_COMMAND, *params)

    def _chunk_request(self, kegg_ids, kegg_command, *params):
        """
        efficiently process request of multiple entries into chunks
        :param kegg_command: kegg command type one of list | link | conv | get
        :param kegg_ids: str or list of kegg ids
        :param params: optional extra parameters for command
        :return: raw data of all kegg_ids with the given the kegg_command
        """
        all_data = ''
        if isinstance(kegg_ids, str):
            kegg_ids = [kegg_ids]
        #  kegg takes 10 values at a time
        chunks = [kegg_ids[i:i + KEGG_MAX_IDS] for i in range(0, len(kegg_ids), KEGG_MAX_IDS)]
        for chunk in chunks:
            query = self.format_multiple_genes(chunk)
            data = self.kegg_command(kegg_command, query, *params)
            all_data += data
        return all_data

    def get_all_genes(self, species=KEGG_HOMO_SAPIENS):
        """
        :param species: str
        :return: list of Kegg <species> genes
        """
        data = self.kegg_command('list', species, '')
        return self._process_response(data)

    def get_all_pathways(self, species=KEGG_HOMO_SAPIENS):
        """
        :param species: str
        :return: list of Kegg <species> pathways
        """
        data = self.kegg_command('list', 'pathway', species)
        return self._process_response(data)

    def get_all_modules(self):
        """
        :return: list of Kegg <species> pathways
        """
        data = self.kegg_command('list', 'module', '')
        return self._process_response(data)

    def module_orthologs(self, module_id):
        """
        list of ortholog involved in kegg module
        :param module_id:
        :return:
        """
        data = self.kegg_command('link', 'ko', module_id)
        orthologs = {ortholog.split(':')[1] for ortholog in self._process_response(data, return_as_set=True)}
        return orthologs

    def ortholog_genes(self, orthologs, species=KEGG_HOMO_SAPIENS):
        """
        :type orthologs: set or string
        :param species:
        :return:
        """
        if isinstance(orthologs, (set, list)):
            orthologs = self.format_multiple_genes(orthologs)
        data = self.kegg_command('link', species, orthologs)
        return self._process_response(data, return_as_set=True)

    @staticmethod
    def _process_response(pathway_data, return_as_list=False, return_as_set=False):
        if return_as_set:  # correction to support older version
            return_as_list = True
        res = {} if not return_as_list else []
        rows = pathway_data.split('\n')
        for row in rows[:-1]:
            items = row.split('\t')
            kegg_id = items[0]
            desc = items[1] if len(items) >= 2 else ''
            if return_as_list:
                res.append(desc)
            else:
                res[kegg_id] = desc
        if return_as_set:
            res = set(res)
        return res

    def get_pathway_info(self, pathway_id):
        """
        retrieves raw pathway information in xml format
        :param pathway_id: str - kegg pathway id
        :return:
        """
        return self.kegg_command('get', pathway_id, 'kgml')

    def get_pathway_gene_list(self, pathway_id):
        """
        retrieves set of genes in a pathway
        :param pathway_id: str kegg pathway id
        :return:
        """
        data = self.kegg_command('link', KEGG_HOMO_SAPIENS, pathway_id)
        genes = self._process_response(data, return_as_set=True)
        return genes

    def get_module_gene_list(self, module_id):
        """
        retrieves set of genes in a module
        :param module_id: str kegg pathway id
        :return:
        """
        orthologs = self.module_orthologs(module_id)
        return self.ortholog_genes(orthologs)

    def get_gene_list(self, kegg_id):
        """
        given general kegg_id retrieves all genes in network
        :param kegg_id:
        :return: set list of kegg gene ids
        """
        if kegg_id.startswith(KEGG_PATHWAY_PREFIX):
            return self.get_pathway_gene_list(kegg_id)
        elif kegg_id.startswith(KEGG_MODULE_PREFIX):
            return self.get_module_gene_list(kegg_id)
        else:
            raise ValueError(NETWORK_ID_ERROR)

    def convert_gene_names(self, genes, database='uniprot'):
        """
        convert gene names between datasets
        :param genes: str or list of gene names
        :param database: one of genes | ncbi-geneid | ncbi-proteinid | uniprot
        :return: dict containing primary entry and secondary entries
        """
        res = {'primary': None, 'secondary': set()}
        if isinstance(genes, (set, list)):
            genes = self.format_multiple_genes(genes)
        data = self.kegg_command('conv', database, genes)
        genes = self._process_response(data, return_as_set=False, return_as_list=True)
        if (not genes) or (genes == EMPTY_LIST):
            return res
        #  first instance is saved as primary
        res['primary'] = genes[0].split(':')[1]
        for gene in genes[1:]:
            res['secondary'].add(gene.split(':')[1])
        return res

    def gene_seq(self, genes, seq_type='aaseq'):
        """
        retrieves amino-acid or dna sequence of gene
        :param genes: str or list of gene names
        :param seq_type: one of aaseq | ntseq
        :return: dict {kegg_id : seq}
        """
        data = self._chunk_request(genes, 'get', seq_type)
        return self._gene_seq_to_dict(data)

    def genes_info(self, genes):
        """
        :param genes: str or list of gene names
        :return: dict {gene_id : dict gene details}
        Note: in rare cases some genes will be skipped
        """
        data = self._chunk_request(genes, 'get', '')
        return self._process_genes_info(data)

    @staticmethod
    def _gene_seq_to_dict(data):
        res = {}
        rows = data.split('\n')[:-1]
        gene_name = ''
        for row in rows:
            if row.startswith('>'):  # header
                gene_name = row.split(' ')[0][1:]
                if gene_name not in res:
                    res[gene_name] = ''
                continue
            else:
                res[gene_name] += row
        return res

    def _process_genes_info(self, data):
        genes_dict = {}
        genes = data.split(GENE_SEPERATOR)[:-1]
        for gene in genes:
            data = self._process_single_gene(gene)
            if not data:
                continue
            genes_dict = genes_dict | data
        return genes_dict

    @staticmethod
    def _process_single_gene(data):
        # define default gene data
        gene_data = copy.deepcopy(GENE_DATA)
        kegg_id = None
        aa_seq_flag, na_seq_flag = False, False
        aa_seq_len, na_seq_len = None, None

        # process each row in the gene data
        for row in data.split('\n'):
            if not row.strip():
                continue
            values = row.split()
            # while flags are True len(values) must be 1
            if len(values) > 1:
                aa_seq_flag, na_seq_flag = False, False
            # append sequences until flag is False
            if aa_seq_flag:
                gene_data['aa_seq'] += values[0]
                continue
            if na_seq_flag:
                gene_data['na_seq'] += values[0]
                continue
            try:
                title = values[0]
            except IndexError:  # happens in rare cases [empty lines, corrupted data...]
                print(f"[KEGG PARSE ERROR] Malformed line in gene block: {row}")
                continue
            # process each title
            if title == 'ENTRY':
                kegg_id = f'hsa:{values[1]}'
                gene_data['coding_type'] = values[2]
            if title == 'SYMBOL':
                gene_data['ref_names'] = values[1:]
            if title == 'POSITION':
                line_data = values[1].split(':')
                gene_data['chr'] = line_data[0]
                try:
                    res = re.search(KEG_POSITION_RE, line_data[1])
                    gene_data['start'], gene_data['end'] = int(res.groups()[0]), int(res.groups()[1])
                except IndexError:
                    gene_data['chr'] = values[1]
            if title.startswith('UniProt'):
                ids = values[1:]
                gene_data['uniprot_id'] = { 'primary': ids[0] if ids else None,
                                            'secondary': set(ids[1:]) if len(ids) > 1 else set()  }
            if title == 'AASEQ':
                aa_seq_flag = True
                aa_seq_len = int(values[1])
                continue
            if title == 'NTSEQ':
                na_seq_flag = True
                na_seq_len = int(values[1])
                continue
        # check if all required fields are present
        if gene_data['aa_seq']:
            assert len(gene_data['aa_seq']) == aa_seq_len, 'aa_seq length does not match ' + kegg_id
        if gene_data['na_seq']:
            assert len(gene_data['na_seq']) == na_seq_len, 'na_seq length does not match ' + kegg_id
        assert kegg_id, 'Unable to find Kegg id'
        gene_data['kegg_id'] = kegg_id
        return {kegg_id: gene_data}

    @staticmethod
    def hugo_to_kegg_hsa(hugo):
        url = f"https://rest.kegg.jp/find/genes/{hugo}"
        r = requests.get(url)
        r.raise_for_status()
        hsa_ids = []
        for line in r.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2 and parts[0].startswith("hsa:"):
                hsa_ids.append(parts[0])
        return hsa_ids



class KeggNetwork:
    """Object to represent KEGG Modules or Pathways using precomputed gene SNV files."""

    def __init__(self, kegg_id, network_type):
        self.id, self.type = kegg_id, network_type.lower()
        self._dict_path = pjoin(KEGG_PATHWAY_OBJECTS_P, f"{self.id}.pickle")

        assert self.type in NETWORK_TYPES, NETWORK_TYPE_ERROR

        if not os.path.exists(self._dict_path):
            raise FileNotFoundError(f"Missing gene SNV mapping for pathway/module: {self.id}")

        self.gene_snv_map: dict = pd.read_pickle(self._dict_path)  # {kegg_id: path_to_snv_csv}
        self.gene_names_list = list(self.gene_snv_map.keys())

    def __len__(self):
        return len(self.gene_names_list)

    @property
    def num_genes(self) -> int:
        """
        Returns the total number of genes in the network.
        This is an alias for len(self).
        """
        return len(self.gene_names_list)

    def total_length(self) -> int:
        """
        Calculates the total length of the amino acid sequences for all genes
        in the network.
        :return: int - sum of all gene AA sequence lengths.
        """
        total_len = 0
        for gene in self.genes:
            try:
                # KeggGene has a __len__ method that returns len(self.aa_seq)
                # KeggGene's initialization handles loading/fetching the sequence data.
                total_len += len(gene)
            except Exception as e:
                # Handles cases where KeggGene might be missing sequence data
                warnings.warn(f"[Warning] Could not get length for gene {gene.kegg_id}: {e}")
                continue
        return total_len

    @property
    def genes(self):
        """Yield KeggGene instances for all genes in the network."""
        for gene_id in self.gene_names_list:
            yield KeggGene(gene_id)

    def all_snvs(self, outpath='', index=True):
        """
        Load precomputed SNVs for all genes in the pathway.
        :param index: bool, include index in final CSV
        :param outpath: optional output CSV path (defaults to KEGG_PATHWAY_MUTATIONS_PATH/<id>.csv)
        :return: DataFrame of concatenated SNVs
        """
        collector = []
        for gene_id in tqdm(self.gene_names_list, desc=f"Reading SNVs for {self.id}", unit="gene"):
            snv_file = self.gene_snv_map.get(gene_id)
            if snv_file and os.path.exists(snv_file):
                try:
                    df = pd.read_csv(snv_file)
                    if not df.empty:
                        collector.append(df)
                except Exception as e:
                    print(f"[Error] Reading {snv_file}: {e}")
            else:
                print(f"[Warning] SNV file not found for {gene_id}")

        if not collector:
            print(f"[Warning] No SNVs found for pathway {self.id}")
            return pd.DataFrame(columns=FAMANALYSIS_COLUMNS)

        all_snvs = pd.concat(collector, ignore_index=True)

        if not outpath:
            outpath = pjoin(KEGG_GENE_SCORES_P, f"{self.id}.csv")

        all_snvs.to_csv(outpath, index=index)
        return all_snvs

    def get_pathway_scores(self):
        """Loads the precomputed background mutation scores from file"""
        scores_path = str(pjoin(KEGG_PATHWAY_SCORES_P, self.id + ".csv"))
        df = pd.read_csv(scores_path)
        return df

    def get_genes_lengths(self) -> List[tuple[str, int]]:
        """Returns a list of tuples of (gene name, length)"""
        return [(gene.kegg_id, len(gene)) for gene in self.genes]

    def get_genes_num_cancer_samples(self, cancer_file: str, cancer_folder: str = CBIO_P)\
            -> List[tuple[str, int]] or None:
        """Returns a list of tuples of (gene name, number of mutations recorded for this cancer)"""
        cancer_path = pjoin(cancer_folder, cancer_file)

        if not os.path.exists(cancer_path):
            warnings.warn(f"Cancer file not found: {cancer_path}")
            # Return list of tuples with count 0 for all genes in the network
            return None

        try:
            # Load only the required KEGG ID column
            df = pd.read_csv(cancer_path, usecols=[KEGG_COL])
        except ValueError as e:
            # Handle case where KEGG_COL might not exist in the CSV (based on the inspiration code's check)
            warnings.warn(f"Error loading {KEGG_COL} from {cancer_file}: {e}")
            return None

        def is_in_pathway(kegg_id_cell: str) -> bool:
            """Checks if any Kegg ID in the cell is part of the network's gene list."""
            return any(gid in self.gene_names_list for gid in str(kegg_id_cell).split(','))

        # Filter the DataFrame to only include mutations relevant to the network
        # The relevant column is KEGG_COL, which may contain multiple IDs separated by ','
        mask = df[KEGG_COL].astype(str).apply(is_in_pathway)
        pathway_cancer_mutations = df[mask].copy()

        # Initialize mutation counts for all genes in the network to zero
        gene_mutation_counts = {gene_id: 0 for gene_id in self.gene_names_list}

        # Count mutations per gene in the network
        for _, row in pathway_cancer_mutations.iterrows():
            # Get all KEGG IDs associated with this mutation
            kegg_ids_in_cell = str(row[KEGG_COL]).split(',')

            # Count the mutation for every gene in the network that is present in the cell
            for gid in kegg_ids_in_cell:
                if gid in gene_mutation_counts:
                    gene_mutation_counts[gid] += 1
                    break

        # Convert the dictionary to the required list of tuples format: (gene name, count)
        result_list = list(gene_mutation_counts.items())

        return result_list

    def get_genes_is_covered(self, cancer_file: str, cancer_folder: str = CBIO_P)\
            -> List[tuple[str, bool]] or None:
        """Calculates for each gene in the network whether it is "covered"
        by the given cancer dataset based on a threshold rule."""

        genes_num_cancer_samples = self.get_genes_num_cancer_samples(
            cancer_file, cancer_folder=cancer_folder)

        if not genes_num_cancer_samples:
            return None

        # Convert list of tuples to dictionary for quick lookup: {kegg_id: count}
        gene_counts = dict(genes_num_cancer_samples)

        # Get gene lengths
        genes_lengths = self.get_genes_lengths()

        # Convert list of tuples to dictionary for quick lookup: {kegg_id: length}
        gene_lengths = dict(genes_lengths)

        def is_covered(gene_id: str) -> bool:
            """Applies the coverage logic for a single gene."""
            gene_num_samples = gene_counts.get(gene_id, 0)
            gene_length = gene_lengths.get(gene_id, 0)  # Use 0 length for safety

            # First condition: Absolute count threshold
            if gene_num_samples >= ABSOLUTE_COUNT_THRESHOLD:
                return True

            # Second condition: Mutation density threshold
            if gene_length > 0:
                # Calculate the percentage: (count / length) * 100
                coverage_percentage = (gene_num_samples / gene_length) * 100
                if coverage_percentage >= COVERAGE_PERCENTAGE_THRESHOLD:
                    return True

            return False

        return [(gene_id, is_covered(gene_id)) for gene_id in self.gene_names_list]

    def get_coverage_percentage(self, cancer_file: str, cancer_folder: str = CBIO_P)\
            -> float or None:
        genes_coverage = self.get_genes_is_covered(cancer_file, cancer_folder=cancer_folder)
        if not genes_coverage:
            return None

        # Get the total number of genes in the network
        total_genes = len(genes_coverage)  # Same as len(self) or self.num_genes

        # Count the number of covered genes (where the tuple's second element is True)
        covered_genes_count = sum(1 for _, is_covered in genes_coverage if is_covered)

        # Calculate the percentage
        if total_genes == 0:
            return 0.0  # Avoid division by zero, though genes_coverage should be non-empty if not None

        return covered_genes_count / total_genes



class KeggGene:
    """Gene instance for KEGG pathway"""
    """
    CDS: 20536
    miRNA: 1913
    ncRNA: 1454
    rRNA: 761
    tRNA: 22
    """

    def __init__(self, kegg_id, redownload=False):
        """Constructor for Protein"""
        self.kegg_id, self.uniprot_id, self.ref_names = kegg_id, None, None
        self.na_seq, self.aa_seq, self.chr, self.start, self.end = None, None, None, None, None
        self.coding_type = None
        self._dir_name = kegg_id.replace(':', '_')
        self._directory = pjoin(KEGG_GENES_P, self._dir_name + '.pickle')
        self.gc_content = None

        if redownload:
            self._create_new_instance(kegg_id)
        if os.path.exists(self._directory):
            load_obj(self, self._directory, name=kegg_id)

    def __len__(self):
        return len(self.aa_seq)

    def _create_new_instance(self, kegg_id):
        kegg_api = KeggApi()

        # 1. Fetch the general info (this contains coding_type, ref_names, chr, etc.)
        info = kegg_api.genes_info(kegg_id)[kegg_id]

        self.uniprot_id = info['uniprot_id']
        self.ref_names = info['ref_names']
        self.chr = info['chr']
        self.start = info['start']
        self.end = info['end']
        self.coding_type = info['coding_type']

        # 2. Fetch sequences
        if self.coding_type == "CDS":
            self.aa_seq = kegg_api.gene_seq(kegg_id, 'aaseq').get(kegg_id)
        self.na_seq = kegg_api.gene_seq(kegg_id, 'ntseq').get(kegg_id)

        # 3. Calculate GC content and save
        self.gc_content = self.get_gc_content()
        save_obj(self, self._directory)

    def create_from_dict(self, data):
        self.__dict__ = data
        self._dir_name = self.kegg_id.replace(':', '_')
        self._directory = pjoin(KEGG_GENES_P, self._dir_name + '.pickle')
        save_obj(self, self._directory)

    def get_gc_content(self):
        if not isinstance(self.na_seq, str) or not self.na_seq:
            return None
        return (self.na_seq.upper().count('G') + self.na_seq.upper().count('C')) / len(self.na_seq)

    @property
    def uid(self):
        """
        primary uniprot id
        :return: str
        """
        if self.uniprot_id and 'primary' in self.uniprot_id:
            return self.uniprot_id['primary']
        return None

    @property
    def alias_uid(self):
        """
        alias uniprot ids does not return the main id
        :return: set
        """
        if self.uniprot_id and 'secondary' in self.uniprot_id:
            return set(self.uniprot_id['secondary'])
        return set()

    def length(self, seq='aa'):
        """
        :param seq: aa | na
        :return: return length of amino acid or nucleic acid sequence
        """
        return len(self.aa_seq) if seq == 'aa' else len(self.na_seq)

    def all_snvs(self, outpath='', index=False):
        """
        Creates a DataFrame of all nonsynonymous single nucleotide variants (SNVs) in the gene.

        Notes:
            - The last codon (3 nucleotides) is skipped (assumed to be the stop codon).
            - If len(self.na_seq) % 3 != 0, we return empty DataFrame.

        :param index: bool, include index column in DataFrame if True
        :param outpath: str, if provided, saves the DataFrame as a CSV to this path
        :return: pandas DataFrame with SNV data
        """
        # SKIP if the nucleic acid sequence is not a multiple of 3
        # if len(self.na_seq) % 3 != 0:
        #     print(f"Gene: {self.kegg_id} has {len(self.na_seq)} nucleotides, <!%3==0>!")
        #     return pd.DataFrame(columns=FAMANALYSIS_COLUMNS)

        # SKIP the sequence if it is not a CDS
        if self.coding_type != "CDS":
            print(f"Gene: {self.kegg_id} is a {self.coding_type}, skip.")
            return

        def read_in_chunks(seq, chunk_size=3):
            """Yield only full codons (3 bases)."""
            for i in range(0, len(seq) - (len(seq) % chunk_size), chunk_size):
                yield seq[i:i + chunk_size]

        row_data = lambda idx, ref_na, alt_na, ref_aa, alt_aa: \
            ['-', idx, idx, ref_na, alt_na, self.uid, f'{ref_aa}{idx}{alt_aa}']

        mutate_codon = lambda codon, idx, alt: codon[:idx] + alt + codon[idx + 1:]

        df = pd.DataFrame(columns=FAMANALYSIS_COLUMNS)

        for chunk, codon in enumerate(read_in_chunks(self.na_seq[:-3], chunk_size=CODON_LENGTH)):  # skip stop codon
            codon = codon.lower()  # ensure lowercase for consistency with CODON_TRANSLATOR
            if codon not in CODON_TRANSLATOR:
                continue  # skip unknown or invalid codons
            ref_aa = CODON_TRANSLATOR[codon]

            for idx, ref_na in enumerate(codon):
                if ref_na not in NA_CHANGE:
                    continue  # skip invalid nucleotide
                index = (CODON_LENGTH * chunk) + idx

                for alt_na in NA_CHANGE[ref_na]:
                    alt_codon = mutate_codon(codon, idx, alt_na)
                    alt_codon = alt_codon.lower()
                    if alt_codon not in CODON_TRANSLATOR:
                        continue  # skip invalid mutated codons
                    alt_aa = CODON_TRANSLATOR[alt_codon]
                    if alt_aa == ref_aa or alt_aa == STOP_AA:
                        continue  # ignore synonymous and nonsense variants

                    df.loc[len(df)] = row_data(index, ref_na, alt_na, ref_aa, alt_aa)

        if outpath:
            df.to_csv(outpath, index=index)

        return df