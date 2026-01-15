# Downloading all cBioPortal TCGA mutations studies

from cbio_api import CbioApi

if __name__ == '__main__':
    print("----- Downloading all cBioPortal TCGA mutations studies... -----")
    # Initialize cBioPortal API
    cbio = CbioApi()

    # Download all TCGA study mutations
    cbio.download_all_tcga_study_mutations()