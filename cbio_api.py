# The cBioPortal API class

from definitions import *
import pandas as pd
from bravado.client import SwaggerClient
import os

snake_format = lambda s: s.replace(' ', '_').replace('-', '_').lower()

class CbioApi:
    """api for cbio portal"""
    # cbio is portal for cancer genomics data

    def __init__(self, ):
        """Constructor for Cbio"""
        cbioportal = SwaggerClient.from_url(CBIO_API_URL,
                                            config={"validate_requests": False, "validate_responses": False,
                                                    "validate_swagger_spec": False})
        for a in dir(cbioportal):
            cbioportal.__setattr__(snake_format(a), cbioportal.__getattr__(a))
        self.api = cbioportal

    @staticmethod
    def set_cbio_api_call():
        cbioportal = SwaggerClient.from_url(CBIO_API_URL,
                                            config={"validate_requests": False, "validate_responses": False,
                                                    "validate_swagger_spec": False})
        for a in dir(cbioportal):
            cbioportal.__setattr__(snake_format(a), cbioportal.__getattr__(a))
        return cbioportal

    def get_sample_cancer_type(self, study_id, sample_id):
        data = self.api.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId=study_id,
                                                                                sampleId=sample_id,
                                                                                attributeId='CANCER_TYPE')
        return snake_format(data.result()[0].value)

    def download_all_tcga_study_mutations(self):
        """
        Downloads mutations for all TCGA studies and saves them as CSV files in path defined by CBIO_MUTATION_STUDIES.
        """
        all_studies = self.api.Studies.getAllStudiesUsingGET().result()
        tcga_studies = [study for study in all_studies if 'tcga' in study.studyId.lower()]
        num_studies = len(tcga_studies)

        for i, study in enumerate(tcga_studies):
            study_id = study.studyId
            print(f"    Downloading mutations for study: {study_id}")

            mutations = self.api.mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
                molecularProfileId=f"{study_id}_mutations",
                # {study_id}_mutations gives default mutations profile for study
                sampleListId=f"{study_id}_all",  # {study_id}_all includes all samples
                projection="DETAILED")  # include gene info
            output_path = pjoin(CBIO_MUTATION_STUDIES, f"{study_id}.csv")

            if os.path.exists(output_path):
                print(f"{study_id} mutations already downloaded ({i+1}/{num_studies}).")
            else:
                try:
                    self.study_to_csv(mutations, str(output_path))
                    print(f"    Saved {study_id}.csv ({i+1}/{num_studies}).")
                except Exception as e:
                    print(f"    Skipping {study_id} due to error: {e}")


    @staticmethod
    def study_to_csv(results, output_path='', remove_duplicates=True):
        """
        :param remove_duplicates:
        :param output_path:
        :param results: bravado.http_future.HttpFuture object
        """
        mutations = results.result()
        data = [(m.chr, m.startPosition, m.endPosition, m.referenceAllele, m.variantAllele, m.gene.hugoGeneSymbol,
                 m.proteinChange, m.patientId, m.uniquePatientKey, m.sampleId, m.studyId, m.ncbiBuild)
                for m in mutations if m.mutationType == MISSENSE_MUTATION]
        df = pd.DataFrame.from_records(data, columns=STUDY_COLUMNS)
        #  drop duplicate mutations of the same patient
        #  mutations can repeat in the same patient in the same study if there are multiple samples per patient
        if remove_duplicates:
            df.drop_duplicates(keep='first', inplace=True, ignore_index=True, subset=DUPLICATE_EXCLUSION_COLUMNS)
        if output_path:
            df.to_csv(output_path)

    def cancer_types_dict(self):
        """
        :return: dict {cancer_type : cbio_short_name}
        """
        all_types = self.api.Cancer_Types.getAllCancerTypesUsingGET().result()
        return {snake_format(cancer_type.name): cancer_type.shortName for cancer_type in all_types}

    def all_studies_by_keyword(self, keyword, outpath=''):
        """
        :param keyword: abbreviated cancer type
        :param outpath: path to save results
        :return: all studies with samples of cancer_type == keyword
        """
        studies = self.api.Studies.getAllStudiesUsingGET(keyword=keyword).result()
        #  make sure cancer type is correct
        studies = [study for study in studies if study.cancerTypeId == keyword]
        study_ids = [study.studyId for study in studies]
        study_names = [study.name for study in studies]
        if outpath:
            with open(outpath, 'w+') as file:
                for id, name in zip(study_ids, study_names):
                    file.write(f"{name} \t {id}\n")
        return study_ids, study_names

    def get_cancer_type(self, cancer_short_name: str) -> str:
        """
        :param cancer_short_name: str abbreviated cancer type
        :return: full cancer type name
        """
        all_types = self.api.Cancer_Types.getAllCancerTypesUsingGET().result()
        for cancer_type in all_types:
            if cancer_type.shortName.lower() == cancer_short_name:
                return cancer_type.name
        return ''

    def get_cancer_short_name(self, cancer_type: str) -> str:
        """
        :param cancer_type: str full cancer type name
        :return: abbreviated cancer type name
        """
        all_types = self.api.Cancer_Types.getAllCancerTypesUsingGET().result()
        for cancer_t in all_types:
            if cancer_t.name.lower() == cancer_type.lower():
                return cancer_t.shortName
        return ''

