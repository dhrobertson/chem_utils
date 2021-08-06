"""
subselect_data.py
------------
Author: Daniel H Robertson

Part of the IBRI informatics system

compile the data needed for the ade tool
"""
import argparse
import os
import sys
import csv
import logging

def read_ids_file(ids_file):

    id_list = list()
    if not os.path.isfile(ids_file):
        logger.critical("Unable to locate defined input file \"{}}\" for ids ... skipping".format(ids_file))
        return []

    # read file
    nread = 0
    with open(ids_file) as file_ref:
        for line in file_ref.readlines():
            id = line.rstrip()
            # hack
            id = id.split(" ")[1]
            id_list.append(id)
            nread += 1
    file_ref.close()

    logger.info("\nRead {} ids from \"{}\"\n\n".format(len(id_list), ids_file))

    return id_list


def subselect_data(in_file, ids_file, id_column, out_file):

    ids_list = read_ids_file(ids_file)

    if not len(ids_list):
        logger.critical("Did not read any ids from input file {} ... exiting".format(ids_file))
        return

    if not os.path.isfile(in_file):
        logger.critical("Unable to locate defined input file \"{}}\" ... exiting".format(in_file))
        exit()

    #
    # either use out_file if defined or stdout otherwise
    #
    out_stream = sys.stdout
    if (out_file):
        out_stream = open(out_file, 'w', newline = '')
    csvwriter = csv.writer(out_stream, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    nrow = 0
    nselected = 0
    with open(in_file) as file_ref:
        csv_reader = csv.reader(file_ref, delimiter=',')
        for row in csv_reader:
            logger.debug("{}".format(row))
            if nrow == 0:
                headers = row
                csvwriter.writerow(row)
                col_no = -1
                for i in range(0, len(row)):
                    if row[i] == id_column:
                        col_no = i
                if col_no == -1:
                    logger.critical("Cannot find id_column {} in headers {} ... exiting".format(id_column, row))
                    file_ref.close()
                    exit()
            else:
                id = row[col_no]
                if id in ids_list:
                    csvwriter.writerow(row)
                    nselected += 1
            nrow += 1
    file_ref.close()
    logger.info("Keep {} of {} rows".format(nselected, nrow-1))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='subselect_data.py',
        description='subselect records to only those within the list.')
    parser.add_argument('-f', '--file', dest='in_file', help='csv file containing full data', required=True)
    parser.add_argument('-i', '--ids_file', dest='ids_file', help='txt file of ids to keep', required=True)
    parser.add_argument('-c', '--column', dest='id_column', help='column with id for subselection', required=True)
    parser.add_argument('-o', '--output', dest='out_file', help='optional output file -- otherwise write to stdout')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", default=False, help='turn on addition information that is sent to stderr')
    parser.add_argument('-d', '--debug', dest='debug', action="store_true", default=False, help='turn on very verbose debug information to stderr')

    args = parser.parse_args()
    # create logger
    logger = logging.getLogger()
    logger.level = logging.WARN
    if args.verbose:
        logger.level = logging.INFO
    if args.debug:
        logger.level = logging.DEBUG
    stream_handler = logging.StreamHandler(sys.stderr)
    logger.addHandler(stream_handler)

    # process the out_file argument
    out_file = args.out_file if args.out_file else ""

    subselect_data(args.in_file, args.ids_file, args.id_column, out_file)
