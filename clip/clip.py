import argparse
import os
import re
import sys
import traceback

from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from bokehCLIP import bokehCLIP
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
    @TODO catch the exception thrown with unsorted files and re call gffc.process(unsoted=True)
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
    Count crosslink sites per feature
    '''
    bedC = bedCLIP(args)
    if args.choice == 'a':
        bedC.count_all()
    elif args.choice == 'o':
        bedC.count_only()

def _countSlidingWindows(args):
    '''
    Count crosslink sites per sliding window
    '''
    bedC = bedCLIP(args)
    bedC.countSlidingWindow()

def _slidingWindowsToDEXSeq(args):
    '''
    Convert sliding window crosslink sites to DEXSeq format
    '''
    bedC = bedCLIP(args)
    bedC.toDEXSeq()

def _junction(args):
    bedC = bedCLIP(args)
    bedC.junction()

'''
Given below is should be the completed help file, functions with a '*' needs to be reworked extensively
    {0}:  A flexible toolset for the analysis of iCLIP and eCLIP sequencing data

    The function (as a positional argument) should be one of:

    [Annotation]
        annotation              flattens an annotation gtf file
        createSlidingWindows    creates sliding windows based on given annotation file
    
    [Extraction]
        extract                 extracts crosslink sites, insertions or deletions
    
    [Counting]
        count                   count sites in annotation
        countSlidingWindows     count sites in sliding windows
        feature*                 count sites in repeated regions
    
    [Transformation]
        slidingWindowsToDEXSeq  transform sliding window counts to DEXSeq format
        
    [Distances]
        junction                calculates distances to junctions
        dist*                    calculates nearest cross link site to a feature

    [Visualisation] 
        plot*                   visualisation 
           
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
    
    [Extraction]
        extract                 extracts crosslink sites, insertions or deletions
    
    [Counting]
        count                   count sites in annotation
        countSlidingWindows     count sites in sliding windows
    
    [Transformation]
        slidingWindowsToDEXSeq  transform sliding window counts to DEXSeq format
        
    [Distances]
        junction                calculates distances to junctions
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
    annotation.add_argument('--unsorted',dest='unsorted',help='use this flag if the GFF file is unsorted',action='store_true')
    # createSlidingWindows
    cshelp = 'createSlidingWindows: creates sliding windows out of the flattened annotation file'
    createSlidingWindows = subps.add_parser('createSlidingWindows',description=cshelp, formatter_class=argparse.RawTextHelpFormatter) # help='create sliding windows',
    createSlidingWindows.add_argument('-i','--input',metavar='input file',dest='input',help='flattend annotation file, see "{} annotation -h"'.format(prog),required=True)
    createSlidingWindows.add_argument('-o','--output',metavar = 'output file',dest='output',help='annotation sliding windows file (.bed[.gz], default: print to console)',default=None,type=str)
    createSlidingWindows.add_argument('-w','--windowSize',metavar = 'window size',dest='windowSize',help='window size (in number of base pairs) for sliding window (default: 50)',default=50,type=int)
    createSlidingWindows.add_argument('-s','--windowStep',metavar = 'step size',dest='windowStep',help='window step size for sliding window (default: 20)',default=20,type=int)
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
    cchoices = ['o','a']
    chelp = 'count: counts the number of crosslink/deletion/insertion sites'
    count = subps.add_parser('count',description=chelp,formatter_class=argparse.RawTextHelpFormatter) #help='count crosslinks',
    count.add_argument('-i','--input',metavar='input bed',help='extracted crosslink, insertion or deletion sites (.bed[.gz]), see "{} extract -h"'.format(prog),required=True)
    count.add_argument('-o','--output',metavar = 'output file',dest='output',help='output count file (.txt[.gz], default: print to console)',default=None,type=str)
    count.add_argument('-a','--ann',metavar = 'annotation',dest='compare',help='flattened annotation file (.bed[.gz]), see "{} annotation -h"'.format(prog),required=True)
    count.add_argument('-c','--count',dest='choice',
        help='Crosslink site count choices, must be one of: {0} \n o: only include features (exons/introns) with crosslink site \n a: include all features, even with 0 crosslink sites. (default: o).'.format(', '.join(cchoices)),choices=cchoices,default='o')
    # countSlidingWindows
    cswhelp = 'countSlidingWindows: counts the number of crosslink/deletion/insertion sites in a certain sliding window'
    countSlidingWindows = subps.add_parser('countSlidingWindows', description=cswhelp, formatter_class=argparse.RawTextHelpFormatter) # help='count sliding windows',
    countSlidingWindows.add_argument('-i','--input',metavar='input bed',help='extracted crosslink, insertion or deletion sites (.bed[.gz]), see "{} extract -h"'.format(prog),required=True)
    countSlidingWindows.add_argument('-o','--output',metavar = 'output file',dest='output',help='output count file (.txt[.gz], default: print to console)',default=None,type=str)
    countSlidingWindows.add_argument('-s','--sw',metavar = 'sliding window',dest='compare',help='sliding window annotation file (.bed[.gz]), see "{} createSlidingWindows -h"'.format(prog),required=True)
    # @TODO: finish counting part
    ''' ____________________ [Transformation] ___________________ '''
    # slidingWindowsToDEXSeq
    # @TODO: create a wrapper for countSlidingWindows and slidingWindowsToDEXSeq
    swhelp = 'slidingWindowsToDEXSeq: transforms the sliding window counts into DEXSeq format'
    slidingWindowsToDEXSeq = subps.add_parser('slidingWindowsToDEXSeq',description=swhelp, formatter_class=argparse.RawTextHelpFormatter) # help='to DEXSeq format',
    slidingWindowsToDEXSeq.add_argument('-i','--input',metavar='input counts',help='sliding window count data file [.gz], see "{} countSlidingWindows -h"'.format(prog),required=True)
    slidingWindowsToDEXSeq.add_argument('-o','--output',metavar = 'output file',dest='output',help='output DEXSeq formatted count (.txt[.gz], default: print to console)',default=None,type=str)
    ''' ____________________ [Distances] ___________________ '''
    # junction 
    jhelp = 'junction: calculates the distance from crosslink/deletion/insertion sites to the junction'
    junction = subps.add_parser('junction',description=jhelp, formatter_class=argparse.RawTextHelpFormatter) # help='crosslink junctions',
    junction.add_argument('-i','--input',metavar='input bed',help='extracted crosslink, insertion or deletion sites (.bed[.gz]), see "{} extract -h"'.format(prog),required=True)
    junction.add_argument('-a','--ann',metavar = 'annotation',dest='compare',help='flattened annotation file (.bed[.gz]), see "{} annotation -h"'.format(prog),required=True)
    junction.add_argument('-o','--output',metavar = 'output file',dest='output',help='output junction file (.txt[.gz], default: print to console)',default=None,type=str)
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
        elif args.subparser == 'extract':
            # extract crosslink sites based on annotations
            _extract(args)
        elif args.subparser == 'count':
            # count extracted crosslink sites
            _count(args)
        elif args.subparser == 'countSlidingWindows':
            # count sliding window from extracted crosslink sites
            _countSlidingWindows(args)
        elif args.subparser == 'slidingWindowsToDEXSeq':
            # convert sliding window counts to DEXSeq format
            _slidingWindowsToDEXSeq(args)
        elif args.subparser == 'junction':
            # generate junction info from extracted crosslink sites
            _junction(args)
    except KeyboardInterrupt:
        sys.stderr.write('Keyboard interrupt... good bye\n')
    except Exception:
        traceback.print_exc(file=sys.stdout)
    sys.exit(0)

if __name__=='__main__':
    main()