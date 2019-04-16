# Nanoanalysis
ont_filter.py mainly corrects the format of ont data and quality control of data.
Instructions:
./ont_filter.py -h
usage: ont_filter.py [-h] -fq FILE [FILE ...] [-qc INT] [-ml INT] [-mqw INT]
                     [-c INT] [-ow] [-of DIRECTORY] [-oc DIRECTORY]

name:

    ont_filter.py -- Quality control of ont data

attention:
    ont_filter.py -i *.fastq -w work -o results
Input file format:
    You need to enter the absolute path of the ont fastq file and the summary.txt file.
    example:
    /export/sda1/Nanopore_raw_data/1.fastq
    /export/sda1/Nanopore_raw_data/2.fastq

optional arguments:
  -h, --help            show this help message and exit
  -fq FILE [FILE ...], --fastq FILE [FILE ...]
                        set the input fastq.
  -qc INT, --control INT
                        Whether to quality control the sequence that removes
                        the wrong format,default=NO.
  -ml INT, --min_length INT
                        min length,default=500.
  -mqw INT, --mean_q_weight INT
                        mean qscore template,default=7.
  -c INT, --concurrent INT
                        number of concurrent process,default=1.
  -ow, --overwrite      set overwrite the workdir/outdir if existed.
  -of DIRECTORY, --out_filter DIRECTORY
                        set the work directory,default=01_raw_data.
  -oc DIRECTORY, --out_clean DIRECTORY
                        set the out directory,default=02_clean_data.


Building Nanopore's assembly notes and other analysis software
