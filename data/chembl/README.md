# data.chembl

## about

Scripts to handle the ChEMBL access to data: molecules, targets, assays, and activities through the API

Reference: https://www.ebi.ac.uk/chembl/api/data/docs

## testing

This toolkit uses pytest and associated tests and data within the tests/ directory as part of its development.

These test should have close to full coverage across all the various objects, methods, and scripts.

    python3 -m pytest tests/data/pharos/test_pharos_utils.py

## dependencies

### Python Modules
  requests
  requests_cache

## TODO
* implement requests_cache
