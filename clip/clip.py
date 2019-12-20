import argparse
import os
import re
import sys
import traceback

from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from bokehCLIP import bokehCLIP
from countCLIP import countCLIP
from createMatrix import MatrixConverter
from featureCLIP import feature
from gffCLIP import gffCLIP
from heatmap import HeatMap


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
    gffc = gffCLIP(args)
    try:
        gffc.process(args.unsorted)
    except SyntaxError as se:
        if args.unsorted:
            raise(se)
        else:
            sys.stderr.write(str(se)+'\n')
            sys.stderr.write('Trying to parse {} with "--unsorted" option.\n'.format(args.gff))
            sys.stderr.write('Warning! this step is memory hungry\n')
            gffc.process(True)

def _createSlidingWindows(args):
    '''
    Create sliding windows from the given annotation file
    '''
    gffc = gffCLIP(args)
    gffc.slidingWindow(args.input)

def _mapToId(args):
    mapC = countCLIP(args)
    mapC.annotationToIDs()

def _extract(args):
    '''
    Extract cross-link sites
    '''
    bamC = bamCLIP(args)
    if args.choice == 's':
        bamC.extract_StartSites(offset=args.offset,ignore=args.ignore)
    elif args.choice == 'i':
        bamC.extract_InsertionSites()
    elif args.choice == 'd':      
        bamC.extract_DeletionSites()
    elif args.choice == 'm': 
        bamC.extract_MiddleSites()
    elif args.choice == 'e':
        bamC.extract_EndSites(offset=args.offset,ignore=args.ignore)

def _count(args):
    '''
    Count crosslink sites per sliding window
    '''
    countC = countCLIP(args)
    stranded = True
    if args.unstranded:
        stranded = False
    countC.count(stranded)

def _junction(args):
    raise NotImplementedError('This function is not yet implemented yet')
    bedC = bedCLIP(args)
    bedC.junction()

def _countMatrix(args):
    mC = MatrixConverter(args.input,args.prefix,args.postfix,args.output)
    mC.read_samples()
    mC.write_matrix()

'''
Given below is should be the completed help file, functions with a '*' needs to be reworked extensively
    {0}:  A flexible toolset for the analysis of iCLIP and eCLIP sequencing data

    The function (as a positional argument) should be one of:

    [Annotation]
        annotation              flattens an annotation gtf file
        createSlidingWindows    creates sliding windows based on given annotation file
        mapToId               map entries in "name" column to unique ids and write in tab separated format 
    
    [Extraction]
        extract                 extracts crosslink sites, insertions or deletions
    
    [Counting]
        count                   count sites in annotation
        feature*                count sites in repeated regions
            
    [Distances]
        junction                calculates distances to junctions
        dist*                    calculates nearest cross link site to a feature

    [Visualisation] 
        plot*                   visualisation 

    [Helpers]
        createMatrix            create R friendly matrix from count function output files

    [General help]
        -h, --help               help
    '''

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
        
    [Distances]
        junction                calculates distances to junctions
    
    [Helpers]
        createMatrix            create R friendly matrix from count function output files

    '''.format(prog)
    epilog = "For command line options of each argument, use: {} <positional argument> -h".format(prog)
    # @TODO: some of the argument names are confusing, needs fixing
    parser = argparse.ArgumentParser(prog=prog, description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)
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
    # createSlidingWindows
    cshelp = 'createSlidingWindows: creates sliding windows out of the flattened annotation file'
    createSlidingWindows = subps.add_parser('createSlidingWindows',description=cshelp, formatter_class=argparse.RawTextHelpFormatter) # help='create sliding windows',
    createSlidingWindows.add_argument('-i','--input',metavar='input file',dest='input',help='flattend annotation file, see "{} annotation -h"'.format(prog),required=True)
    createSlidingWindows.add_argument('-o','--output',metavar = 'output file',dest='output',help='annotation sliding windows file (.bed[.gz], default: print to console)',default=None,type=str)
    createSlidingWindows.add_argument('-w','--windowSize',metavar = 'window size',dest='windowSize',help='window size (in number of base pairs) for sliding window (default: 50)',default=50,type=int)
    createSlidingWindows.add_argument('-s','--windowStep',metavar = 'step size',dest='windowStep',help='window step size for sliding window (default: 20)',default=20,type=int)
    # mapToIds
    maphelp = 'mapToId: extract "name" column from the annotation file and map the entries to unique id and print out in tab separated format'
    mapToId = subps.add_parser('mapToId',description=maphelp, formatter_class = argparse.RawTextHelpFormatter)
    mapToId.add_argument('-a','--annotation',metavar= 'annotation file', help = 'flattened annotation file from "{0} annotation -h" or sliding window file from "{0} createSlidingWindows -h"'.format(prog),required=True)
    mapToId.add_argument('-o','--output',metavar = 'output file',dest='output',help='region/window annotation mapped to a unique id (.txt[.gz], default: print to console)',default=None,type=str)

    ''' ____________________ [Extraction] ___________________ '''
    # extract
    ehelp = 'extract:  extracts crosslink sites, insertions or deletions'
    echoices = ['s','i','d','m','e']
    mates = [1,2]
    extract = subps.add_parser('extract',description=ehelp,formatter_class=argparse.RawTextHelpFormatter) #,help='extract crosslinks'
    extract.add_argument('-i','--input', metavar='input file',dest='input',help='input file (.bam)',required=True)
    extract.add_argument('-o','--output', metavar = 'output file',dest='output',help='output file (.bed, default: print to console)',default=None,type=str)
    extract.add_argument('-e','--mate', dest='mate',help='for paired end sequencing, select the read/mate to extract the crosslink sites from.\n Must be one of: {}'.format(', '.join([str(i) for i in mates])),type=int,choices=mates,required=True) # make it required ?
    extract.add_argument('-s','--site',dest='choice',
        help='Crosslink site choices, must be one of: {0}\n s: start site \n i: insertion site \n d: deletion site \n m: middle site \n e: end site (default: e).'.format(', '.join(echoices)),choices=echoices,default='e')
    extract.add_argument('-g','--offset',metavar='offset length',dest='offset',help='Number of nucleotides to offset for crosslink sites (default: 0)',type=int,default=0)
    extract.add_argument('--ignore',dest='ignore',help='flag to ignore crosslink sites outside of genome',action='store_true')
    extract.add_argument('-q','--minAlignmentQuality',metavar = 'min. alignment quality',dest='minAlignmentQuality',help='minimum alignment quality (default: 10)',type=int,default=10)
    extract.add_argument('-m','--minReadLength',metavar='min. read length',dest='minReadLength',help='minimum read length (default: 0)',type=int,default=0)
    extract.add_argument('-x','--maxReadLength',metavar='max. read length',dest='maxReadLength',help='maximum read length (default: 0)',type=int,default=0)
    extract.add_argument('-l','--maxReadInterval',metavar='max. read interval',dest='maxReadIntervalLength',help='maximum read interval length (default: 10000)',type=int,default=10000)
    extract.add_argument('--primary',dest='primary',help='flag to use only primary positions of multimapping reads',action='store_true')
    ''' ____________________ [Counting] ___________________ '''
    # count
    chelp = 'count: counts the number of crosslink/deletion/insertion sites'
    count = subps.add_parser('count',description=chelp,formatter_class=argparse.RawTextHelpFormatter) #help='count crosslinks',
    count.add_argument('-i','--input',metavar='input bed',help='extracted crosslink, insertion or deletion sites (.bed[.gz]), see "{} extract -h"'.format(prog),required=True)
    count.add_argument('-o','--output',metavar = 'output file',dest='output',help='output count file (.txt[.gz], default: print to console)',default=None,type=str)
    count.add_argument('-a','--ann',metavar = 'annotation',dest='annotation',help='flattened annotation file (.bed[.gz]), see "{0} annotation -h" OR sliding window annotation file (.bed[.gz]), see "{0} createSlidingWindows -h"'.format(prog),required=True)
    count.add_argument('--unstranded',dest='unstranded', help='by default, crosslink site counting is strand specific. Use this flag for non strand specific crosslink site counting',action='store_true')
    ''' ____________________ [Distances] ___________________ '''
    # junction 
    jhelp = 'junction: calculates the distance from crosslink/deletion/insertion sites to the junction'
    junction = subps.add_parser('junction',description=jhelp, formatter_class=argparse.RawTextHelpFormatter) # help='crosslink junctions',
    junction.add_argument('-i','--input',metavar='input bed',help='extracted crosslink, insertion or deletion sites (.bed[.gz]), see "{} extract -h"'.format(prog),required=True)
    junction.add_argument('-a','--ann',metavar = 'annotation',dest='compare',help='flattened annotation file (.bed[.gz]), see "{} annotation -h"'.format(prog),required=True)
    junction.add_argument('-o','--output',metavar = 'output file',dest='output',help='output junction file (.txt[.gz], default: print to console)',default=None,type=str)
    ''' ____________________ [Helpers] ___________________ '''
    # createMatrix
    cmhelp = 'createMatrix: create R friendly output matrix file from count function output files'
    createMatrix = subps.add_parser('createMatrix',description=cmhelp,formatter_class = argparse.RawTextHelpFormatter)
    createMatrix.add_argument('-i','--inputFolder', dest='input', metavar = 'input folder', help='Folder name with output files from count function, see "{} count -h ", supports .gz (gzipped files)'.format(prog), required = True)
    createMatrix.add_argument('-b','--prefix', dest='prefix', metavar = 'file name prefix', help='Use files only with this given file name prefix (default: None)', default="", type=str)
    createMatrix.add_argument('-e','--postfix', dest='postfix', metavar = 'file name postfix', help='Use files only with this given file name postfix (default: None). WARNING! either "--prefix" or "--postfix" argument must be given!', default="", type=str)
    createMatrix.add_argument('-o','--output',metavar = 'output file',dest='output',help='output junction file (.txt[.gz], default: print to console)',default=None,type=str)

    # # distance
    # dhelp = 'distance: calculates the nearest crosslink site to a region/feature'
    # distance = subps.add_parser('distance',description=dhelp) # help='nearest crosslink site',
    # # feature ?
    # fhelp = 'feature:  counts the number of crosslink/deletion/insertion sites'
    # feature = subps.add_parser('feature',description=fhelp) # help='count crosslinks',
    
    # # plot
    # phelp  = 'plot: plots htseq-clip analysis results'
    # plot = subps.add_parser('plot',description=phelp) # help = 'plot results',
    # # heatmap
    # hhelp = 'heatmap: plots the heatmap results'
    # heatmap = subps.add_parser('heatmap',description=hhelp) # help='plot heatmap',
    # Now read in arguments and process
    try:
        args = parser.parse_args()
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
        elif args.subparser == 'junction':
            # generate junction info from extracted crosslink sites
            _junction(args)
        elif args.subparser == 'createMatrix':
            # collect output files from count function and generate an R friendly matrix
            if args.prefix == '' and args.postfix == '':
                createMatrix.print_help()
                raise argparse.ArgumentTypeError('Input values for both arguments "--prefix" and "--postfix" cannot be empty! Either one of the values MUST be given')
            _countMatrix(args)
    except KeyboardInterrupt:
        sys.stderr.write('Keyboard interrupt... good bye\n')
        sys.exit(1)
    except Exception:
        traceback.print_exc(file=sys.stdout)
        sys.exit(1)
    sys.exit(0)

if __name__=='__main__':
    main()
