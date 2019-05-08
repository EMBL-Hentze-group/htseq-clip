# --------------------------------------------------
# htseq-clip main
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# --------------------------------------------------
    
import argparse, traceback, os, sys
from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from gtfCLIP import gtfCLIP
from gtf import gtfClip
from bokehCLIP import bokehCLIP
#from fastaCLIP import fastaCLIP
from featureCLIP import feature
from gffCLIP import gffClip
from heatmap import HeatMap

VERSION = "0.1.1"

#======================================================================================
# changed function description from slidingWindows to createSlidingWindows
def usage():
    print '''
htseq-clip:  A flexible toolset for the analysis of iCLIP and eCLIP sequencing data
usage:       htseq-clip <command> [options]
    
The functions include:
    
[Annotation]
    annotation              flattens an annotation gtf file
    createSlidingWindows           creates sliding windows based on given annotation file
    
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
    
[In development]
    genomeToReads           splits up an genome fasta to reads for mappability tests 
    
[General help]
    -h, --help               help
    --version                version
    
'''
    
def usage_extract():
    print '''
htseq-clip extract:  extracts crosslink sites, insertions or deletions
usage:               htseq-clip extract [options]

Options:

 -i, --input                      input file (.bam)
 
 -o, --output                     output file (.bed)
 
 -c, --choice                     parameter for the choice of the cross-link sites:

                                  s for start sites. 

                                  s<int>[i] position + offset e.g.
                                     s-1 start site minus one nucleotide
                                     s-1i ignore crosslink sites outside of genome
                                  
                                  m for middle sites (center of reads)
                                  
                                  e<int>[i] for end sites
                                     e+1 end site plus one nucleotide 
                                     e+1i ignores crosslink sites outside of genome  
                                                                 
                                  i for insertion sites
                                  
                                  d for deletion sites

 -e, --mate                       select which read to extract the crosslink sites from
                                  
                                  f for first read
                                  s for second read
                                  
 -q, --minAlignmentQuality        minimum alignment quality (default: 10)
 
 -m, --maxReadLength              minimum read length (default: 0)
 
 -x, --maxReadLength              maximum read length (default: 0)
 
 -l, --maxReadIntervalLength      maximum read interval length (default: 100000)
 
 -p, --primary                    set if you want only primary positions of 
                                  multimapping reads (default: False)
  
 -h, --help                       help
 --version                        version
'''    
    
def usage_annotation():
    print '''
htseq-clip annotation:  flattens the given annotation file
usage:                  htseq-clip annotation [options]

Options:

 -g, --gtf         GTF or GFF3 file for annotation processing (.gtf[.gz] / .gff[.gz] / .gff3[.gz])
 
 -o, --output      output file (.bed[.gz])
 
 -t, --type        Gene type identifier in gtf / gff file

 -r, --region      True if you want exons to be split into cds and utr regions. Default is False.
                               
 -h, --help        help
 --version         version
'''
    
def usage_count():
    print '''
htseq-clip count:  counts the number of crosslink/deletion/insertion sites 
usage:             htseq-clip count [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed[.gz])
 
 -o, --output        output count file (.txt[.gz])
 
 -f, --compare       flattened annotation file (.bed[.gz])
 
 -c, --choice        parameter for the choice of included counts:
 
                     a if you want to have all counts included
                     even if a chromosome or strand does not
                     contain any crosslink sites
                     
                     o if you want have only the counts included
                     of exons/introns which possess crosslink
                     sites

                              
 -h, --help          help
 --version           version
'''  

def usage_dist():
    print '''
htseq-clip distance:  calculates the nearest cross link site to a region/feature
usage:                htseq-clip dist [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed)

 -o, --output        output file (.txt)

 -f, --compare       file containing regions/features for example rmsk file for repeat elements (.bed)

 -sc, --score        If both files contain alignment score as 5th column. Default is 'y'. If no score information is provided then 0 is
                     written in score column of output file.

 -st, --strand       If both files contain strand information as 6th column. Default is 'y'. If no information about strand is provided
                     then '.' is written in strand column of output file representing either of the two strands (+/-).

 -h, --help          help
 --version           version
'''

def usage_feature():
    print '''
htseq-clip feature:  counts the number of crosslink/deletion/insertion sites
usage:             htseq-clip count [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed)

 -o, --output        output file (.txt)

 -f, --compare       features file e.g. rmsk file for repeat elements (.bed)

 -c, --choice        parameter for the choice of included counts:

                     a if you want to have all counts included
                     even if a chromosome or strand does not
                     contain any crosslink sites

                     o if you want to have only the counts included
                     of exons/introns which possess crosslink
                     sites

 -sc, --score        If both files contain alignment score as 5th column. Default is 'y'. If no score information is provided then 0 is
                     written in score column of output file.

 -st, --strand       If both files contain strand information as 6th column. Default is 'y'. If no information about strand is provided
                     then '.' is written in strand column of output file representing either of the two strands (+/-).

 -h, --help          help
 --version           version
'''

def usage_junction():
    print '''
htseq-clip junction:  calculates the distance from crosslink/deletion/insertion sites to the junction
usage:                htseq-clip junction [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed)
 
 -o, --output        output file (.txt)
 
 -f, --compare       flattened annotation file (.bed)

                             
 -h, --help          help
 --version           version
'''

def usage_dist():
    print '''
htseq-clip distance:  calculates the nearest cross link site to a region/feature
usage:                htseq-clip dist [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed)

 -o, --output        output file (.txt)

 -f, --compare       file containing regions/features for example rmsk file for repeat elements (.bed)

 -sc, --score        If both files contain alignment score as 5th column. Default is 'y'. If no score information is provided then 0 is
                     written in score column of output file.

 -st, --strand       If both files contain strand information as 6th column. Default is 'y'. If no information about strand is provided
                     then '.' is written in strand column of output file representing either of the two strands (+/-).

 -h, --help          help
 --version           version
'''
    
def usage_createSlidingWindows():
    print '''
htseq-clip createSlidingWindows:     creates sliding windows out of the flattened annotation file 
usage:                               htseq-clip createSlidingWindows [options]

Options:

 -i, --input           flattened annotation file (.bed)
 
 -o, --output          output file (.bed)
 
 -w, --windowSize      window size for sliding window
 
 -s, --windowStep      window step for sliding window

                             
 -h, --help            help
 --version             version
''' 
    
def usage_countSlidingWindows():
    print '''
htseq-clip countSlidingWindows:    counts the number of crosslink/deletion/insertion sites in a certain sliding window
usage:                             htseq-clip countSlidingWindows [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed)
 
 -o, --output        output file (.txt)
 
 -f, --compare       created sliding windows out of the flattened annotation file (.bed)

                              
 -h, --help          help
 --version           version
'''  
    
def usage_slidingWindowsToDEXSeq():
    print '''
htseq-clip slidingWindowsToDEXSeq:    transforms the sliding window counts into DEXSeq format
usage:                                htseq-clip slidingWindowsToDEXSeq [options]

Options:

 -i, --input         counted number of crosslink/deletion/insertion sites in a certain sliding window (.txt)
 
 -o, --output        output file (.txt)

                             
 -h, --help          help
 --version           version
'''      
    
    
def usage_plot():
    print '''
htseq-clip plot:  plots the results
usage:            htseq-clip plot [options]

Options:

 -i, --input         count or junction input file (.bed[.gz], .txt[.gz])
 
 -o, --output        output file (.html)
 
 -c, --choice        parameter for the choice of the desired plot:
 
                     c for count information plots
                     
                     j for junction information plots
                     
                     r read information plots

                              
 -h, --help          help
 --version           version
'''  

def usage_heat():
    print '''
htseq-clip heatmap:  plots the heatmap
usage:            htseq-clip heatmap [options]

Options:

 -i, --input         input files. 2 .txt files; 1st input is cross link count inside the feature.
                     2nd input is cross link sites outside (upstream and downstream) the feature.

 -o, --output        output files. 2 output files. 1st is .txt for storing the densities of cross link
                     in each feature. 2nd is .pdf for saving the heatmap.

 -h, --height        height of heatmap in inches. Default is 5 inches.

 -wd, --width        width of heatmap in inches. Default is 5 inches.

 -e, --element       Feature for which heatmap is to be generated. This is case sensitive, so make sure
                     to use same case for the feature as in the input file.


 -h, --help          help
 --version           version
'''

#======================================================================================
#======================================================================================
'''
Counting the cross-link sites per feature
'''    
def count(program, parser, args):
    
    bedC = bedCLIP(args)
    
    if hasattr(args, 'choice') and args.choice != None:
        if args.choice == 'a':
            bedC.count_all() 
        elif args.choice == 'o':      
            bedC.count_only()
        else:
            parser.error('Invalid option for count')
    else:
        parser.error('You need -c option for the correct counting of your data!')
        
'''
Sliding window counts
updated function name to countSlidingWindows to match with the rest of the main module
'''            
def countSlidingWindows(args):
    
    bedC = bedCLIP(args)
    bedC.countSlidingWindow()

'''
Extract cross-link sites
'''
def extract(parser,args):
    
    bamC = bamCLIP(args)
    

    if hasattr(args, 'choice') and args.choice != None:
        if args.choice.startswith("s") :
            if hasattr(args,'mate') and args.mate != None:
                (ignore, offset) = bamC.extractOptions(args.choice.split("s"))
                bamC.extract_StartSites(offset, ignore)
            else:
                parser.error('You need --mate option. Indicate which read to extract cross linking sites from') 
        elif args.choice == 'm': 
            bamC.extract_MiddleSites()
        elif args.choice.startswith("e"):   
            (ignore, offset) = bamC.extractOptions(args.choice.split("e"))
            bamC.extract_EndSites(offset, ignore)
        elif args.choice == 'd':      
            bamC.extract_DeletionSites()
        elif args.choice == 'i':      
            bamC.extract_InsertionSites()
        else:
            parser.error("""Invalid option for parameter -c (can be s, m, e, d or i).
Please find a detailed description in the usage (htseq-clip extract --help)""")
    else:
        parser.error('You need to specify the -c option for extraction.')
        
'''
Processing genome.fa file into fasta format 
@TODO: this function needs biopython!
'''    
def genomeToReads(args):
    
    fastaC = fastaCLIP(args)   
    fastaC.genomeToReads()

'''
Calculate distance from cross-link site to exon/intron junction site
'''        
def junction(program, args):
    
    bedC = bedCLIP(args)
    bedC.junction()

'''
Plotting function
'''

def plot(parser, args):

    bokehC = bokehCLIP(args)

    if hasattr(args, 'choice') and args.choice != None:
        if args.choice == 'c':
            bokehC.plot_count()
        elif args.choice == 'd':
            bokehC.plot_comparison()
        elif args.choice == 'j':
            bokehC.plot_junction()
        elif args.choice == 'r':
            bokehC.plot_read()

        else:
            parser.error('Invalid option for plotting')
    else:
        parser.error('You need -c option for the correct plotting of your data!')

'''
Processing the annotation file
'''       
def process(args):
    if args.gtf.endswith('.gtf'):
        gtf = gtfClip(args)
        gtfC = gtfCLIP(args)
        if args.region:
            gtf.processGTF()
        else:
            gtfC.processGTF()
    elif args.gtf.endswith('.gff') or args.gtf.endswith('.gff3'):
        gff = gffClip(args)
        gtfC = gtfCLIP(args)
        if args.region:
            gff.processGFF()
        else:
            gtfC.processGTF()

'''
Processing the annotation file into sliding windows
'''    
def slidingWindow(args):
    
    gtfC = gtfCLIP(args)
    gtfC.slidingWindow()

'''
Converting the sliding window counts into DEXSeq format
'''    
def toDEXSeq(args):
    
    bedC = bedCLIP(args)
    bedC.toDEXSeq()
"""
Calculate distance to nearest cross link site
"""

def cl_dist(parser, args):
    feat = feature(args)
    feat.dist_cl()

"""
Function for counting cross link sites inside repeat regions
"""
def features(parser, args):
    feat = feature(args)
    if hasattr(args, 'choice') and args.choice != None:
        if args.choice == 'a':
            feat.count_all()
        elif args.choice == 'o':
            feat.count_only()
        else:
            parser.error('Invalid option for count')
    else:
        parser.error('You need -c option for the correct counting of your data!')

"""
Function for creating heatmap for visualizing the cross link trends.
"""
def heatmap(parser,args):
    heat = HeatMap(args)
    heat.heatmap()
#======================================================================================
#-------------------------------------------------------------
def checkFileExists(filename, parser):
    if (notEmpty(filename)):
        parser.error ('%s is not a valid filename' % filename)
    if (os.path.isfile(filename) == False):
        parser.error ('%s does not exist' % filename)

#-------------------------------------------------------------
def notEmpty(filename):
    return filename == False or filename == ''

'''
        test_parser = argparse.ArgumentParser(version="htseq-clip", conflict_handler='resolve')
        test_parser.add_argument('-v', '--verbose', action='store_true', default=argparse.SUPPRESS, help='verbose output')
        test_parser.add_argument('-o', '--output', action='store',nargs = '+', type= str,  default='/tmp/test_count.csv', dest='output', help='output file name')
        test_parser.add_argument('-f', '--compare', action='store', type= str,  default=argparse.SUPPRESS, dest='compare', help='file which you want to compare with your input file')
        test_parser.add_argument('-c', '--choice', action='store', type=str,  default='', dest='choice', help='option')
        test_parser.add_argument('-a', '--mate', action='store', type=int, default=1, dest='mate', help='select first or second read')
        test_parser.add_argument('-q', '--minAlignmentQuality',   action='store', type=int, default=argparse.SUPPRESS,   dest='minAlignmentQuality', help='minimum alignment quality')
        test_parser.add_argument('-m', '--minReadLength', action='store', type=int, default=argparse.SUPPRESS, dest='minReadLength', help='minimum read length')
        test_parser.add_argument('-x', '--maxReadLength', action='store', type=int, default=argparse.SUPPRESS, dest='maxReadLength', help='maximum read length')
        test_parser.add_argument('-h', '--height', action='store', type=int, default=argparse.SUPPRESS, dest='height', help='height of heatmap')
        test_parser.add_argument('-wd', '--width', action='store', type=int, default=argparse.SUPPRESS, dest='width', help='width of heatmap')
        test_parser.add_argument('-l', '--maxReadIntervalLength', action='store', type=int, default=argparse.SUPPRESS, dest='maxReadIntervalLength', help='maximum read interval length')
        test_parser.add_argument('-p', '--primary', action='store_true', default=argparse.SUPPRESS,  dest='primary', help='set if you want to use only primary positions of multimapping reads')
        test_parser.add_argument('-d', '--dist', action='store', type=int, default=argparse.SUPPRESS, dest='dist', help='Maximum distance between two sites of a read if its higher than this value it will be count as this value')
        test_parser.add_argument('-g', '--gtf', action='store', type=str, default=argparse.SUPPRESS, dest='gtf', help='gtf file for annotation')
        test_parser.add_argument('-t', '--type', action='store', type=str, default=argparse.SUPPRESS, dest='type', help='gene type for annotation')
        test_parser.add_argument('-n', '--name', action='store', type=str, default=argparse.SUPPRESS, dest='name', help='gene name for annotation')
        test_parser.add_argument('-sc', '--score', action='store', type=str, default=argparse.SUPPRESS, dest='score', help='alignment score')
        test_parser.add_argument('-st', '--strand', action='store', type=str, default=argparse.SUPPRESS, dest='strand', help='region strand')
        test_parser.add_argument('-w', '--windowSize', action='store', type=int, default=argparse.SUPPRESS, dest='windowSize', help='window size for sliding window')
        test_parser.add_argument('-s', '--windowStep', action='store', type=int, default=argparse.SUPPRESS, dest='windowStep', help='window step for sliding window')
        test_parser.add_argument('-r', '--region', action='store_true', dest='region', help='set if you want exons to be split into cds and utr regions. Default is False.')
        test_parser.add_argument('-sc', '--score', action='store', type=str, default=argparse.SUPPRESS, dest='score', help='alignment score')
        test_parser.add_argument('-st', '--strand', action='store', type=str, default=argparse.SUPPRESS, dest='strand', help='region strand')
        test_parser.add_argument('-e', '--element', action='store_true', dest='element', help='Region for heatmap e.g. Alu.')
        test_parser.add_argument('-h', '--height', action='store', type=int, default=argparse.SUPPRESS, dest='height', help='height of heatmap')
        test_parser.add_argument('-wd', '--width', action='store', type=int, default=argparse.SUPPRESS, dest='width', help='width of heatmap')
        test_parser.add_argument('command', nargs = '?', help='name of program to run ')
        test_parser.add_argument('-i', '--input', action='store',nargs = '+', type= str,  default='/home/sahadeva/hentze/projects/Meta-Analysis/CLIPRNAlength/clip_encode_bam/input/ENCFF027WDY.bam', dest='input', help='input file')
'''
def main():
    
    try:
        # get parser and add arguments
        parser = argparse.ArgumentParser(version="htseq-clip {}".format(VERSION), conflict_handler='resolve')

        parser.add_argument('-v', '--verbose', action='store_true', default=argparse.SUPPRESS, help='verbose output')
        parser.add_argument('-i', '--input', action='store',nargs = '+', type= str,  default=argparse.SUPPRESS, dest='input', help='input file')
        parser.add_argument('-o', '--output', action='store',nargs = '+', type= str,  default=argparse.SUPPRESS, dest='output', help='output file name')
        parser.add_argument('-f', '--compare', action='store', type= str,  default=argparse.SUPPRESS, dest='compare', help='file which you want to compare with your input file')
        parser.add_argument('-c', '--choice', action='store', type=str,  default=argparse.SUPPRESS, dest='choice', help='option')
        parser.add_argument('-a', '--mate', action='store', type=int, default=1, dest='mate', help='select first or second read')
        parser.add_argument('-q', '--minAlignmentQuality',   action='store', type=int, default=argparse.SUPPRESS,   dest='minAlignmentQuality', help='minimum alignment quality')
        parser.add_argument('-m', '--minReadLength', action='store', type=int, default=argparse.SUPPRESS, dest='minReadLength', help='minimum read length')
        parser.add_argument('-x', '--maxReadLength', action='store', type=int, default=argparse.SUPPRESS, dest='maxReadLength', help='maximum read length')
        parser.add_argument('-h', '--height', action='store', type=int, default=argparse.SUPPRESS, dest='height', help='height of heatmap')
        parser.add_argument('-wd', '--width', action='store', type=int, default=argparse.SUPPRESS, dest='width', help='width of heatmap')
        parser.add_argument('-l', '--maxReadIntervalLength', action='store', type=int, default=argparse.SUPPRESS, dest='maxReadIntervalLength', help='maximum read interval length')
        parser.add_argument('-p', '--primary', action='store_true', default=argparse.SUPPRESS,  dest='primary', help='set if you want to use only primary positions of multimapping reads')
        parser.add_argument('-d', '--dist', action='store', type=int, default=argparse.SUPPRESS, dest='dist', help='Maximum distance between two sites of a read if its higher than this value it will be count as this value')
        parser.add_argument('-g', '--gtf', action='store', type=str, default=argparse.SUPPRESS, dest='gtf', help='gtf file for annotation')
        parser.add_argument('-t', '--type', action='store', type=str, default=argparse.SUPPRESS, dest='type', help='gene type for annotation')
        parser.add_argument('-n', '--name', action='store', type=str, default=argparse.SUPPRESS, dest='name', help='gene name for annotation')
        parser.add_argument('-sc', '--score', action='store', type=str, default=argparse.SUPPRESS, dest='score', help='alignment score')
        parser.add_argument('-st', '--strand', action='store', type=str, default=argparse.SUPPRESS, dest='strand', help='region strand')
        parser.add_argument('-w', '--windowSize', action='store', type=int, default=argparse.SUPPRESS, dest='windowSize', help='window size for sliding window')
        parser.add_argument('-s', '--windowStep', action='store', type=int, default=argparse.SUPPRESS, dest='windowStep', help='window step for sliding window')
        parser.add_argument('-r', '--region', action='store_true', dest='region', help='set if you want exons to be split into cds and utr regions. Default is False.')
        parser.add_argument('-sc', '--score', action='store', type=str, default=argparse.SUPPRESS, dest='score', help='alignment score')
        parser.add_argument('-st', '--strand', action='store', type=str, default=argparse.SUPPRESS, dest='strand', help='region strand')
        parser.add_argument('-e', '--element', action='store_true', dest='element', help='Region for heatmap e.g. Alu.')
        parser.add_argument('-h', '--height', action='store', type=int, default=argparse.SUPPRESS, dest='height', help='height of heatmap')
        parser.add_argument('-wd', '--width', action='store', type=int, default=argparse.SUPPRESS, dest='width', help='width of heatmap')
        parser.add_argument('command', nargs = '?', help='name of program to run ')
        args = parser.parse_args()
        d = vars(args)
        # Print the default usage when there is no program specified
        if args.command == None:
            usage()
            os._exit(1)
        # Print program specific messages  
        else:
            program = args.command
            # changed all occurences of len(d) <3 to len(d) <=4, by default d = {'region': False, 'mate': 1, 'command': 'slidingWindowsToDEXSeq', 'element': False}
            # so, len(d)<3 never statisfies for missing inputs of any sub parsers
            subpars = ['extract','countSlidingWindows','junction','count','annotation','plot','createSlidingWindows','slidingWindowsToDEXSeq','genomeToReads','feature','dist','heatmap']
            if program == 'extract':
                if len(d) <= 4:
                    usage_extract()
                    os._exit(1)
                else:   
                    checkFileExists(args.input[0], parser)
                    extract(parser, args)
            elif program == 'countSlidingWindows':
                if len(d) <= 4:
                    usage_countSlidingWindows()
                    os._exit(1)
                else:
                    checkFileExists(args.input[0], parser)
                    checkFileExists(args.compare, parser)
                    countSlidingWindows(args)
            elif program == 'junction':
                if len(d) <= 4:
                    usage_junction()
                    os._exit(1)
                else: 
                    checkFileExists(args.input[0], parser)
                    checkFileExists(args.compare, parser)
                    junction(program, args)
            elif program == 'count':
                if len(d) <= 4:
                    usage_count()
                    os._exit(1)
                else: 
                    checkFileExists(args.input[0],parser)
                    checkFileExists(args.compare, parser)
                    count(program, parser, args)
            elif program == 'annotation':
                if len(d) <= 4:
                    usage_annotation()
                    os._exit(1)
                else:
                    checkFileExists(args.gtf, parser)
                    process(args)
            elif program == 'plot':
                if len(d) <= 4:
                    usage_plot()
                    os._exit(1)
                else:
                    checkFileExists(args.input[0], parser)
                    plot(parser, args)
            elif program =='createSlidingWindows':
                if len(d) <= 4:
                    usage_createSlidingWindows()
                    os._exit(1)
                else: 
                    checkFileExists(args.input[0], parser)
                    slidingWindow(args)
            elif program =='slidingWindowsToDEXSeq':
                if len(d) <= 4:
                    usage_slidingWindowsToDEXSeq()
                    os._exit(1)
                else:
                    checkFileExists(args.input[0], parser)
                    toDEXSeq(args)
            elif program =='genomeToReads':
                checkFileExists(args.input[0], parser)
                genomeToReads(args)
            elif program == 'feature':
                if len(d) <= 4:
                    usage_feature()
                    os._exit(1)
                else:
                    checkFileExists(args.input[0], parser)
                    checkFileExists(args.compare, parser)
                    features(parser,args)
            elif program == 'dist':
                if len(d) <= 4:
                    usage_dist()
                    os._exit(1)
                else:
                    checkFileExists(args.input[0], parser)
                    checkFileExists(args.compare, parser)
                    cl_dist(parser,args)

            elif program == 'heatmap':
                if len(d) <3:
                    usage_heat()
                    os.exit(1)
                else:
                    checkFileExists(args.input[0],parser)
                    checkFileExists(args.input[1],parser)
                    heatmap(parser,args)

            else:
                parser.error ('Incorrect function argument "{0}". Function must be one of: {1}'.format(program,", ".join(subpars)))

    # Exception Handling: Interruption / Errors
    except KeyboardInterrupt, exception: # Control - C
        raise exception
    except SystemExit, exception:
        raise exception
    except Exception, exception:
        print '______________________________[ Error ]______________________________ '
        print str(exception)
        print '____________________________[ Traceback ]____________________________  '
        traceback.print_exc()

        os._exit(1)

#=====================================================================================

# program can be used as stand-alone or as module
if __name__ == '__main__':
    main()
    

