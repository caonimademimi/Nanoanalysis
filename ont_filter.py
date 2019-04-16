#!/usr/bin/env python

import re
import os
import sys
import json
import gzip
import logging
import os.path
import argparse

from multiprocessing import Pool


def checkDir(out_filter, out_clean, overwrite):
#Check if the folder exists, if it exists, reminder, if it is not pure, create a new folder
    if os.path.exists(out_filter):
        if not overwrite:
            logging.error('This out filter [' + out_filter + '] is existed, please check.')
            sys.exit(2)
        else:
            logging.warning('This out filter [' + out_filter + '] is existed, will be overwrited.')
            os.system("rm -rf " + out_filter)
    if os.path.exists(out_clean):
        if not overwrite:
            logging.error('This out clean [' + out_clean + '] is existed, please check.')
            sys.exit(2)
        else:
            logging.warning('This out clean [' + out_clean + '] is existed, will be overwrited.')
            os.system("rm -rf " + out_clean)
    try :
        os.makedirs(out_filter)
        os.makedirs(out_clean)
    except OSError, err:
        logging.error('Error: %s' %(str(err)))
        sys.exit(2)


def jusge_fastq_format(yield_string):
#Determine if the format of fastq is correct
    for lines in yield_string:

        lines = lines.strip()
        lines = lines.splitlines()

        if len(lines)==4 and len(lines)>=1:
            if len(lines[1].split('\x00'))==1 and len(lines[1])==len(lines[3]):
                if lines[2].startswith("+") and lines[0].startswith("@"):
                    yield lines[0],lines[1],lines[2],lines[3]
                else:
                    print lines[0]
            else:
                print lines[0]
        else:
            print lines[0]



def choose_fastq_quality(yield_string, min_length, mean_q_weight):
#Choose a sequence that fits the quality and length

    for lines in yield_string:  

        lines = lines.strip()
        lines = lines.splitlines()
        mean_quality = 0

        if len(lines)==4:

            quality = [ord(i) for i in lines[3]]
            max_quality = max(quality)
            mean_quality = sum(quality)*1.0/len(quality)

            if max_quality>=0 and max_quality<=33:

                mean_quality = mean_quality

            elif max_quality>=34 and max_quality<=75:

                mean_quality = mean_quality-33

            if mean_quality>=mean_q_weight and len(lines[1])>=min_length:
                yield lines[0],lines[1],lines[2],lines[3]



def yield_read_fastq(stream):
    """
    yield fastq records from stream
    :param stream: a stream object
    :return:
    """

    n = 0
    string = ""

    for line in stream:

        line = line.strip()
        lines = line.split()

        if not line:
            continue
        if lines[0].startswith("@") and len(lines)>=2:
            if len(string)>=1:
               yield string

            string = ""

            string += "%s\n" % line
        else:
            string += "%s\n" % line
    if len(string)>=1:
        yield string


def output_fastq(yield_fastq, filter_fastq):

    out_fastq = open(filter_fastq, 'w') 
    
    for line in yield_fastq:
        out_fastq.write(line[0]+"\n"+line[1]+"\n"+line[2]+"\n"+line[3]+"\n")
    out_fastq.close()


def open_fastq(input_fastq):

    input_name = os.path.abspath(input_fastq)
    mode = 'r'

    if input_name.endswith(".gz"):

        stream = gzip.open(input_name, mode)

    else:
        stream = open(input_name,mode)
    yield_fastq = yield_read_fastq(stream)

    return yield_fastq


def run_raw_to_clean(input_fastq, filter_fastq, clean_fastq, min_length, mean_q_weight, control):

    yield_fastq = open_fastq(input_fastq)
    yield_fastq = jusge_fastq_format(yield_fastq)

    output_fastq(yield_fastq,filter_fastq)

    if control.upper() == "YES":
        yield_fastq = open_fastq(filter_fastq)
        yield_fastq = choose_fastq_quality(yield_fastq,min_length,mean_q_weight)

        output_fastq(yield_fastq,clean_fastq)
    elif control.upper() == "NO":
        print("Does not quality control the removal of malformed sequences")
    else:
        print("Your output may be incorrect, we will not quality control the sequence that removes the formatting error.")


def read_files(fastqs):
#Read the feature type file under the file

    fastq_list = []

    for files in fastqs:
        if  os.path.isfile(files):
            fastq_list.append(files)
    return fastq_list


def run_pool(fastq_list, out_filter, concurrent, out_clean, min_length, mean_q_weight, control):

    pool = Pool(processes=concurrent)
    results=[]
    i = 0

    for line in fastq_list:

        line = line.strip().split()
        i = i+1

        if len(line)==1:

            input_fastq = line[0]
            filter_fastq = out_filter+'/'+str(i)+'.fastq'
            clean_fastq = out_clean+'/'+str(i)+'_clean.fastq'

            results.append(pool.apply_async(run_raw_to_clean,(input_fastq,filter_fastq,clean_fastq,min_length,mean_q_weight, control)))
    pool.close()
    pool.join()


def stat_filter(args):

    fastq_list = read_files(args.fastq)
    

    out_filter, out_clean, concurrent, min_length, mean_q_weight, overwrite = [os.path.abspath(args.out_filter), os.path.abspath(args.out_clean), args.concurrent, args.min_length, args.mean_q_weight, args.overwrite]
    checkDir(out_filter, out_clean, overwrite)
    run_pool(fastq_list, out_filter, concurrent, out_clean, min_length, mean_q_weight, args.control)

    
def filter_args(parser):
    parser.add_argument('-fq', '--fastq', metavar='FILE', nargs='+', type=str, required=True,
        help='set the input fastq.')
    parser.add_argument('-qc', '--control', metavar='INT', type=str,  default='NO',
        help='Whether to quality control the sequence that removes the wrong format,default=NO.')
    parser.add_argument("-ml", "--min_length", metavar='INT', type=int, default=500,
        help="min length,default=500.")
    parser.add_argument("-mqw", "--mean_q_weight", metavar='INT', type=int, default=7,
        help="mean qscore template,default=7.")
    parser.add_argument("-c", "--concurrent", metavar='INT', type=int, default=1,
        help="number of concurrent process,default=1.")
    parser.add_argument('-ow', '--overwrite', action='store_true',
        help = 'set overwrite the workdir/outdir if existed.')
    parser.add_argument('-of', '--out_filter', metavar='DIRECTORY', default='01_raw_data',
        help='set the work directory,default=01_raw_data.')
    parser.add_argument('-oc', '--out_clean', metavar='DIRECTORY', default='02_clean_data',
        help='set the out directory,default=02_clean_data.')
    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:

    ont_filter.py -- Quality control of ont data

attention:
    ont_filter.py -i input.info -w work -o results
Input file format:
    You need to enter the absolute path of the ont fastq file and the summary.txt file.
    example:
    /export/sda1/Nanopore_raw_data/1.fastq
    /export/sda1/Nanopore_raw_data/2.fastq
''')
    args = filter_args(parser).parse_args()

    stat_filter(args)


if __name__ == "__main__":
    main()

