# Downloading clinical data on patients form all cBioPortal TCGA mutations studies

from cbio_api import CbioApi

if __name__ == '__main__':
    print("----- Downloading clinical data on patients from all cBioPortal TCGA mutations studies... -----")
    # Initialize cBioPortal API
    cbio = CbioApi()

    # Download all TCGA study mutations
    cbio.download_all_tcga_patient_clinical()