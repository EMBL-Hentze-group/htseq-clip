import argparse
import logging
import os
import re
import sys
import traceback
from datetime import datetime

from .bamCLIP import bamCLIP
from .countCLIP import countCLIP
from .createMatrix import MatrixConverter
from .gffCLIP import FeatureOrderException, gffCLIP


'''
--------------------------------------------------
htseq-clip main
Authors: Marko Fritz, marko.fritz@embl.de
         Thomas Schwarzl, schwarzl@embl.de
         Nadia Ashraf, nadia.ashraf@embl.de
Modified by: Sudeep Sahadevan, sudeep.sahadevan@embl.de
Institution: EMBL Heidelberg
Date: October 2015
--------------------------------------------------
'''

def _annotation(args):
    '''
    Parse annotations from given GFF file
    @TODO use logging module
    '''
    logging.info('Paring annotations') 
    logging.info('GFF file {}, output file {}'.format(args.gff,args.output))
    gffc = gffCLIP(args)
    try:
        gffc.process(args.unsorted)
    except FeatureOrderException as se:
        if args.unsorted:
            raise(se)
        else:
            logging.warning(str(se))
            logging.warning('Trying to parse {} with "--unsorted" option.'.format(args.gff))
            logging.warning('This step is memory hungry')
            gffc.process(True)

def _createSlidingWindows(args):
    '''
    Create sliding windows from the given annotation file
    '''
    logging.info('Create sliding windows')
    logging.info('input file {}, output file {}'.format(args.input,args.output))
    logging.info('Window size {} step size {}'.format(args.windowSize,args.windowStep))
    gffc = gffCLIP(args)
    gffc.slidingWindow(args.input)

def _mapToId(args):
    logging.info('Creating mapping file from annotations')
    logging.info('Input file {} output file {}'.format(args.annotation,args.output))
    mapC = countCLIP(args)
    mapC.annotationToIDs()

def _extract(args):
    '''
    Extract cross-link sites
    '''
    if args.choice == 's':
        logging.info('Extracting start sites')
        logging.info('Bam file : {}, output file: {}, offset: {}'.format(args.input,args.output,args.offset))
        with bamCLIP(args) as bh:
            bh.extract_start_sites(offset = args.offset)
    elif args.choice == 'i':
        logging.info('Extracting insertion sites')
        logging.info('Bam file : {}, output file: {}'.format(args.input,args.output))
        with bamCLIP(args) as bh:
            bh.extract_insertion_sites()
    elif args.choice == 'd':
        logging.info('Extracting deletion sites')
        logging.info('Bam file : {}, output file: {}'.format(args.input,args.output))
        with bamCLIP(args) as bh:
            bh.extract_deletion_sites()
    elif args.choice == 'm':
        logging.info('Extracting middle sites')
        logging.info('Bam file : {}, output file: {}'.format(args.input,args.output))
        with bamCLIP(args) as bh:
            bh.extract_middle_sites()
    elif args.choice == 'e':
        logging.info('Extracting end sites')
        logging.info('Bam file : {}, output file: {}, offset: {}'.format(args.input,args.output,args.offset))
        with bamCLIP(args) as bh:
            bh.extract_end_sites(offset = args.offset)

def _count(args):
    '''
    Count crosslink sites per sliding window
    '''
    logging.info('Count crosslink sites')
    logging.info('Annotation file {} crosslink sites file {} output file {}'.format(args.annotation,args.input,args.output))
    countC = countCLIP(args)
    stranded = True
    if args.unstranded:
        stranded = False
    countC.count(stranded)

def _countMatrix(args):
    logging.info('Generate matrix from files')
    logging.info('Input folder {}, output file {}'.format(args.input,args.output))
    mC = MatrixConverter(args.input,args.prefix,args.postfix,args.output)
    mC.read_samples()
    mC.write_matrix()

logger = logging.getLogger()

def main():
    prog = 'htseq-clip'
    description = '''
    {0}:  A flexible toolset for the analysis of iCLIP and eCLIP sequencing data

    The function (as a positional argument) should be one of:

    [Annotation]
        annotation              flattens a gff formatted annotation file
        createSlidingWindows    creates sliding windows based on given annotation file
        mapToId                 map entries in "name" column to unique ids and write in tab separated format
    
    [Extraction]
        extract                 extracts crosslink sites, insertions or deletions
    
    [Counting]
        count                   count sites in annotation
    
    [Helpers]
        createMatrix            create R friendly matrix from count function output files

    '''.format(prog)
    epilog = "For command line options of each argument, use: {} <positional argument> -h".format(prog)
    parser = argparse.ArgumentParser(prog=prog, description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)
    # log levels
    loglevels = ['debug','info','warn','quiet']
    # subparsers
    subps = parser.add_subparsers(help='Need positional arguments',dest='subparser')

    ''' ____________________ [Annotation] ___________________  '''
    # annotation
    ahelp = 'annotation: flattens (to BED format) the given annotation file (in GFF format)'
    annotation = subps.add_parser('annotation',description=ahelp, formatter_class=argparse.RawTextHelpFormatter) # help='flatten annotation',
    annotation.add_argument('-g','--gff',metavar='annotation',dest='gff',help='GFF formatted annotation file, supports gzipped (.gz) files',required=True)
    annotation.add_argument('-o','--output',metavar = 'output file',dest='output',help='output file (.bed[.gz], default: print to console)',default=None,type=str)
    annotation.add_argument('-u','--geneid',metavar='gene id',dest='id',help='Gene id attribute in GFF file (default: gene_id for gencode gff files)',default='gene_id',type=str)
    annotation.add_argument('-n','--genename',metavar='gene name',dest='name',help='Gene name attribute in GFF file (default: gene_name for gencode gff files)',default='gene_name',type=str)
    annotation.add_argument('-t','--genetype',metavar='gene type',dest='type',help='Gene type attribute in GFF file (default: gene_type for gencode gff files)',default='gene_type',type=str)
    annotation.add_argument('--splitExons',dest='splitExons',help='use this flag to split exons into exonic features such as 5\'UTR, CDS and 3\' UTR',action='store_true')
    annotation.add_argument('--unsorted',dest='unsorted',help='use this flag if the GFF file is unsorted',action='store_true')
    annotation.add_argument('-v','--verbose',metavar='Verbose level',dest='log',help='Allowed choices: '+', '.join(loglevels)+' (default: info)',choices=loglevels,default='info')
    # createSlidingWindows
    cshelp = 'createSlidingWindows: creates sliding windows out of the flattened annotation file'
    createSlidingWindows = subps.add_parser('createSlidingWindows',description=cshelp, formatter_class=argparse.RawTextHelpFormatter) # help='create sliding windows',
    createSlidingWindows.add_argument('-i','--input',metavar='input file',dest='input',help='flattend annotation file, see "{} annotation -h"'.format(prog),required=True)
    createSlidingWindows.add_argument('-o','--output',metavar = 'output file',dest='output',help='annotation sliding windows file (.bed[.gz], default: print to console)',default=None,type=str)
    createSlidingWindows.add_argument('-w','--windowSize',metavar = 'window size',dest='windowSize',help='window size (in number of base pairs) for sliding window (default: 50)',default=50,type=int)
    createSlidingWindows.add_argument('-s','--windowStep',metavar = 'step size',dest='windowStep',help='window step size for sliding window (default: 20)',default=20,type=int)
    createSlidingWindows.add_argument('-v','--verbose',metavar='Verbose level',dest='log',help='Allowed choices: '+', '.join(loglevels)+' (default: info)',choices=loglevels,default='info')
    # mapToIds
    maphelp = 'mapToId: extract "name" column from the annotation file and map the entries to unique id and print out in tab separated format'
    mapToId = subps.add_parser('mapToId',description=maphelp, formatter_class = argparse.RawTextHelpFormatter)
    mapToId.add_argument('-a','--annotation',metavar= 'annotation file', help = 'flattened annotation file from "{0} annotation -h" or sliding window file from "{0} createSlidingWindows -h"'.format(prog),required=True)
    mapToId.add_argument('-o','--output',metavar = 'output file',dest='output',help='region/window annotation mapped to a unique id (.txt[.gz], default: print to console)',default=None,type=str)
    mapToId.add_argument('-v','--verbose',metavar='Verbose level',dest='log',help='Allowed choices: '+', '.join(loglevels)+' (default: info)',choices=loglevels,default='info')

    ''' ____________________ [Extraction] ___________________ '''
    # extract
    ehelp = 'extract:  extracts crosslink sites, insertions or deletions'
    echoices = ['s','i','d','m','e']
    mates = [1,2]
    extract = subps.add_parser('extract',description=ehelp,formatter_class=argparse.RawTextHelpFormatter) #,help='extract crosslinks'
    extract.add_argument('-i','--input', metavar='input file',dest='input',help='input file (.bam, MUST be co-ordinate sorted and indexed)',required=True)
    extract.add_argument('-o','--output', metavar = 'output file',dest='output',help='output file (.bed, default: print to console)',default=None,type=str)
    extract.add_argument('-e','--mate', dest='mate',help='for paired end sequencing, select the read/mate to extract the crosslink sites from.\n Must be one of: {}'.format(', '.join([str(i) for i in mates])),type=int,choices=mates,required=True) # make it required ?
    extract.add_argument('-s','--site',dest='choice',
        help='Crosslink site choices, must be one of: {0}\n s: start site \n i: insertion site \n d: deletion site \n m: middle site \n e: end site (default: e).'.format(', '.join(echoices)),choices=echoices,default='e')
    extract.add_argument('-g','--offset',metavar='offset length',dest='offset',help='Number of nucleotides to offset for crosslink sites (default: 0)',type=int,default=0)
    extract.add_argument('--ignore',dest='ignore',help='flag to ignore crosslink sites outside of genome',action='store_true')
    extract.add_argument('-q','--minAlignmentQuality',metavar = 'min. alignment quality',dest='minAlignmentQuality',help='minimum alignment quality (default: 10)',type=int,default=10)
    extract.add_argument('-m','--minReadLength',metavar='min. read length',dest='minReadLength',help='minimum read length (default: 0)',type=int,default=0)
    extract.add_argument('-x','--maxReadLength',metavar='max. read length',dest='maxReadLength',help='maximum read length (default: 500)',type=int,default=500)
    extract.add_argument('-l','--maxReadInterval',metavar='max. read interval',dest='maxReadIntervalLength',help='maximum read interval length (default: 10000)',type=int,default=10000)
    extract.add_argument('--primary',dest='primary',help='flag to use only primary positions of multimapping reads',action='store_true')
    extract.add_argument('-c','--cores',dest='cores',metavar='cpus',help='Number of cores to use for alignment parsing (default: 5)',default=5,type=int)
    extract.add_argument('-t','--tmp',dest='tmp',metavar='tmp',help='Path to create and store temp files (default behavior: use folder from "--output" parameter)',default=None,type=str)
    extract.add_argument('-v','--verbose',metavar='Verbose level',dest='log',help='Allowed choices: '+', '.join(loglevels)+' (default: info)',choices=loglevels,default='info')

    ''' ____________________ [Counting] ___________________ '''
    # count
    chelp = 'count: counts the number of crosslink/deletion/insertion sites'
    count = subps.add_parser('count',description=chelp,formatter_class=argparse.RawTextHelpFormatter) #help='count crosslinks',
    count.add_argument('-i','--input',metavar='input bed',dest='input',help='extracted crosslink, insertion or deletion sites (.bed[.gz]), see "{} extract -h"'.format(prog),required=True)
    count.add_argument('-o','--output',metavar = 'output file',dest='output',help='output count file (.txt[.gz], default: print to console)',default=None,type=str)
    count.add_argument('-a','--ann',metavar = 'annotation',dest='annotation',help='flattened annotation file (.bed[.gz]), see "{0} annotation -h" OR sliding window annotation file (.bed[.gz]), see "{0} createSlidingWindows -h"'.format(prog),required=True)
    count.add_argument('--unstranded',dest='unstranded', help='crosslink site counting is strand specific by default. Use this flag for non strand specific crosslink site counting',action='store_true')
    count.add_argument('-v','--verbose',metavar='Verbose level',dest='log',help='Allowed choices: '+', '.join(loglevels)+' (default: info)',choices=loglevels,default='info')
    
    ''' ____________________ [Helpers] ___________________ '''
    # createMatrix
    cmhelp = 'createMatrix: create R friendly output matrix file from count function output files'
    createMatrix = subps.add_parser('createMatrix',description=cmhelp,formatter_class = argparse.RawTextHelpFormatter)
    createMatrix.add_argument('-i','--inputFolder', dest='input', metavar = 'input folder', help='Folder name with output files from count function, see "{} count -h ", supports .gz (gzipped files)'.format(prog), required = True)
    createMatrix.add_argument('-b','--prefix', dest='prefix', metavar = 'file name prefix', help='Use files only with this given file name prefix (default: None)', default="", type=str)
    createMatrix.add_argument('-e','--postfix', dest='postfix', metavar = 'file name postfix', help='Use files only with this given file name postfix (default: None). WARNING! either "--prefix" or "--postfix" argument must be given!', default="", type=str)
    createMatrix.add_argument('-o','--output',metavar = 'output file',dest='output',help='output junction file (.txt[.gz], default: print to console)',default=None,type=str)
    createMatrix.add_argument('-v','--verbose',metavar='Verbose level',dest='log',help='Allowed choices: '+', '.join(loglevels)+' (default: info)',choices=loglevels,default='info')

    # Now read in arguments and process
    try:
        args = parser.parse_args()
        if args.subparser is None:
            parser.print_help(sys.stderr)
            sys.exit(1)
        # set logging level and handler
        if args.log== 'quiet':
            logger.addHandler(logging.NullHandler())
        else:
            logger.setLevel(logging.getLevelName(args.log.upper()))
            if len(logger.handlers)>=1:
                # ugly fix for multiple logging handlers
                logger.handlers = []
            consHandle = logging.StreamHandler(sys.stderr)
            consHandle.setLevel(logging.getLevelName(args.log.upper()))
            consHandle.setFormatter(logging.Formatter(' [%(levelname)s]  %(message)s'))
            logger.addHandler(consHandle)
        logging.info('run started at {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M')))
        # check subparsers
        if args.subparser == 'annotation':
            # parse annotations
            _annotation(args)
        elif args.subparser == 'createSlidingWindows':
            # create sliding windows
            _createSlidingWindows(args)
        elif args.subparser == 'mapToId':
            _mapToId(args)
        elif args.subparser == 'extract':
            # extract crosslink sites based on annotations
            _extract(args)
        elif args.subparser == 'count':
            # count extracted crosslink sites
            _count(args)
        elif args.subparser == 'createMatrix':
            # collect output files from count function and generate an R friendly matrix
            if args.prefix == '' and args.postfix == '':
                createMatrix.print_help()
                raise argparse.ArgumentTypeError('Input values for both arguments "--prefix" and "--postfix" cannot be empty! Either one of the values MUST be given')
            _countMatrix(args)
        logging.info('run completed at {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M')))
    except KeyboardInterrupt:
        sys.stderr.write('Keyboard interrupt... good bye\n')
        sys.exit(1)
    except Exception:
        traceback.print_exc(file=sys.stdout)
        sys.exit(1)
    sys.exit(0)

if __name__=='__main__':
    main()
