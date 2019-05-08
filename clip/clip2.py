import argparse
import os
import re
import sys
import traceback

from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from bokehCLIP import bokehCLIP
#from fastaCLIP import fastaCLIP
from featureCLIP import feature
from gffCLIP import gffClip
from gtf import gtfClip
from gtfCLIP import gtfCLIP
from heatmap import HeatMap

'''
--------------------------------------------------
htseq-clip main
Authors: Marko Fritz, marko.fritz@embl.de
         Thomas Schwarzl, schwarzl@embl.de
         Nadia Ashraf, nadia.ashraf@embl.de
Institution: EMBL Heidelberg
Date: October 2015
Modified by: Sudeep Sahadevan, sudeep.sahadevan@embl.de
--------------------------------------------------
'''


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


def main(argv):
    prog = 'htseq-clip'
    description = '''
    {0}:  A flexible toolset for the analysis of iCLIP and eCLIP sequencing data

    The functions (as positional arguments) should be one of:

    [Annotation]
        annotation              flattens an annotation gtf file
        createSlidingWindows    creates sliding windows based on given annotation file
    
    [Extraction]
        extract                 extracts crosslink sites, insertions or deletions
    
    [Counting]
        count                   count sites in annotation
        countSlidingWindows     count sites in sliding windows
        feature                 count sites in repeated regions
        
    [Distances]
        junction                calculates distances to junctions
        dist                    calculates nearest cross link site to a feature

    [Visualisation] 
        plot                    visualisation 
        
    [Transformation]
        slidingWindowsToDEXSeq  transform sliding window counts to DEXSeq format
        
    [General help]
        -h, --help               help
    '''.format(prog)
    epilog = "For command line options of each argument, use: {} <positional argument> -h".format(prog)
    # @TODO: some of the argument names are confusing, needs fixing
    parser = argparse.ArgumentParser(prog=prog, description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)
    subps = parser.add_subparsers(help='Need positional arguments',dest='subparser')
    # extract
    ehelp = 'extract:  extracts crosslink sites, insertions or deletions'
    echoices = ['s','i','d','m','e']
    mates = [1,2]
    extract = subps.add_parser('extract',description=ehelp,formatter_class=argparse.RawTextHelpFormatter) #,help='extract crosslinks'
    extract.add_argument('-i','--input',metavar='input file',dest='input',help='input file (.bam)',required=True)
    extract.add_argument('-o','--output',metavar = 'output file',dest='output',help='output file (.bed, default: print to console)',default=None,type=str)
    extract.add_argument('-e','--mate',metavar='mate info',dest='mate',help='select which read to extract the crosslink sites from. Allowed choices: {} (default: 1)'.format(', '.join([str(i) for i in mates])),choices=mates,default=1,type=int) # make it required ?
    extract.add_argument('-s','--site',dest='choice',
        help='Crosslink site choices: Allowed choices: {0}\n s: start site \n i: insertion site \n d: deletion site \n m: middle site \n e: end site (default: e).'.format(', '.join(echoices)),choices=echoices,default='e')
    extract.add_argument('-g','--offset',metavar='offset length',dest='offset',help='Number of nucleotides to offset for crosslink sites (default: 0)',type=int,default=0)
    extract.add_argument('--ignore',dest='ignore',help='flag to ignore crosslink sites outside of genome',action='store_true')
    extract.add_argument('-q','--minAlignmentQuality',metavar = 'min. alignment quality',dest='minAlignmentQuality',help='minimum alignment quality (default: 10)',type=int,default=10)
    extract.add_argument('-m','--minReadLength',metavar='min. read length',dest='minReadLength',help='minimum read length (default: 0)',type=int,default=0)
    extract.add_argument('-x','--maxReadLength',metavar='max. read length',dest='maxReadLength',help='maximum read length (default: 0)',type=int,default=0)
    extract.add_argument('-l','--maxReadInterval',metavar='max. read interval',dest='maxReadIntervalLength',help='maximum read length (default: 0)',type=int,default=0)
    extract.add_argument('--primary',dest='primary',help='flag to use only primary positions of multimapping reads',action='store_true')
    # annotation
    ahelp = 'annotation: flattens the given annotation file'
    annotation = subps.add_parser('annotation',description=ahelp) # help='flatten annotation',
    annotation.add_argument('-g','--gff',metavar='annotation',dest='gff',help='GFF3 file for annotation processing supports gzipped (.gz) files',required=True)
    annotation.add_argument('-o','--output',metavar = 'output file',dest='output',help='output file (.bed/.bed.gz, default: print to console)',default='',type=str)
    # count
    cchoices = ['o','a']
    chelp = 'count: counts the number of crosslink/deletion/insertion sites'
    count = subps.add_parser('count',description=chelp,formatter_class=argparse.RawTextHelpFormatter) #help='count crosslinks',
    count.add_argument('-i','--input',metavar='input bed',help='extracted crosslink, insertion or deletion sites (.bed[.gz])',required=True)
    count.add_argument('-o','--output',metavar = 'output file',dest='output',help='output count file (.txt[.gz], default: print to console)',default=None,type=str)
    count.add_argument('-f','--compare',metavar = 'annotation',dest='compare',help='flattened annotation file (.bed[.gz])',required=True)
    count.add_argument('-c','--count',dest='choice',
        help='Crosslink site count choices: Allowed choices: {0} \n o: only include features (exons/introns) with crosslink site \n a: include all features, even with 0 crosslink sites. (default: o).'.format(', '.join(cchoices)),choices=cchoices,default='o')
    # distance
    dhelp = 'distance: calculates the nearest crosslink site to a region/feature'
    distance = subps.add_parser('distance',description=dhelp) # help='nearest crosslink site',
    # feature ?
    fhelp = 'feature:  counts the number of crosslink/deletion/insertion sites'
    feature = subps.add_parser('feature',description=fhelp) # help='count crosslinks',
    # junction 
    jhelp = 'junction: calculates the distance from crosslink/deletion/insertion sites to the junction'
    junction = subps.add_parser('junction',description=jhelp) # help='crosslink junctions',
    # createSlidingWindows
    cshelp = 'createSlidingWindows: creates sliding windows out of the flattened annotation file'
    createSlidingWindows = subps.add_parser('createSlidingWindows',description=cshelp) # help='create sliding windows',
    # countSlidingWindows
    cswhelp = 'countSlidingWindows: counts the number of crosslink/deletion/insertion sites in a certain sliding window'
    createSlidingWindows = subps.add_parser('countSlidingWindows', description=cswhelp) # help='count sliding windows',
    # slidingWindowsToDEXSeq
    swhelp = 'slidingWindowsToDEXSeq: transforms the sliding window counts into DEXSeq format'
    slidingWindowsToDEXSeq = subps.add_parser('slidingWindowsToDEXSeq',description=swhelp) # help='to DEXSeq format',
    # plot
    phelp  = 'plot: plots htseq-clip analysis results'
    plot = subps.add_parser('plot',description=phelp) # help = 'plot results',
    # heatmap
    hhelp = 'heatmap: plots the heatmap results'
    heatmap = subps.add_parser('heatmap',description=hhelp) # help='plot heatmap',
    

    # Now read in arguments and process
    try:
        args = parser.parse_args()
        if args.subparser=='extract':
            # extract crosslink sites based on annotations
            print(args)
            print(hasattr(args,'output'))
            # print(len(args.output))
            _extract(args)
        elif args.subparser=='count':
            _count(args)
    except KeyboardInterrupt:
        sys.stderr.write('Keyboard interrupt... good bye\n')
    except Exception:
        traceback.print_exc(file=sys.stdout)
    sys.exit(0)

if __name__=='__main__':
    main(sys.argv)