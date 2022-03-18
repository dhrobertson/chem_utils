"""
test_probeminer_tsv.py
------------
Author: Daniel H Robertson

Part of the IBRI cheminformatics system

test the utility functions for accessing ProbeMiner information through downloaded tsv files
"""
from data import probeminer_tsv
import time
import json
import logging

# create logger
logger = logging.getLogger()
logger.level = logging.DEBUG

# tests

# PAK1 UNIPROT_ACCESSION Q13153
# PAK2 UNIPROT_ACCESSION Q13177
# BRD4 UNIPROT ACCESSION O60885
_pak1_uniprot_ = 'Q13153'
_pak2_uniprot_ = 'Q13177'
_brd4_uniprot_ = 'O60885'


def test_core_functions():
    """ test sending the core functions to read the downloaded probeminer tsv file """

    t0 = time.perf_counter()
    # no cache
    data = probeminer_tsv._read_tsv_file()
    assert 'targets' in data
    assert 'compounds' in data
    assert len(data['targets']) == 3038
    assert len(data['compounds']) == 554989
    t1 = time.perf_counter()
    logger.debug("elapsed time: {}".format(t1-t0))
    t0 = t1

    # should be cached
    data2 = probeminer_tsv._read_tsv_file()
    assert 'targets' in data2
    assert 'compounds' in data2
    assert len(data2['targets']) == 3038
    assert len(data2['compounds']) == 554989
    t1 = time.perf_counter()
    logger.debug("elapsed time: {}".format(t1-t0))
    t0 = t1
    t1 = time.perf_counter()
    logger.debug("elapsed time: {}".format(t1-t0))
    t0 = t1

    # force fail to see timings/debug messages
    #assert True == False

def test_wrapper_functions():
    """ test the wrapper functions to return results in json for one or
        more targets """

    # single value -- returns list of possible probes for single target
    results = probeminer_tsv.get_target_info(_pak1_uniprot_) # as single value
    assert type(results) == type(list())
    assert len(results) == 1
    assert len(results[0]) == 254
    #print(results)

    # list value -- returns list of lists of possible probes
    results2 = probeminer_tsv.get_target_info([_pak1_uniprot_, _pak2_uniprot_]) # as list
    assert type(results2) == type(list())
    assert len(results2) == 2
    assert len(results2[0]) == 254
    assert len(results2[1]) == 53
    #print(results2)

    # list value -- returns list of lists of possible probes
    results3 = probeminer_tsv.get_target_info([_brd4_uniprot_]) # as list
    assert type(results3) == type(list())
    assert len(results3) == 1
    assert len(results3[0]) == 1806
    #print(results3)

    # force fail for debugging
    #assert True == False
