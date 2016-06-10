# --------------------------------------------------
# htseq-clip main
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# --------------------------------------------------
    
import argparse, traceback, os
from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from gtfCLIP import gtfCLIP
from gtf import gtfClip
from bokehCLIP import bokehCLIP
from fastaCLIP import fastaCLIP
from featureCLIP import feature

VERSION = "0.1.1"

#======================================================================================
def usage():
    print '''
htseq-clip:  A flexible toolset for the analysis of iCLIP sequencing data
usage:       htseq-clip <function> [options]
    
The functions include:
    
[annotation]
    annotation              flattens an annotation gtf file
    slidingWindow           creates sliding windows based on given annotation file
    
[iCLIP]
    extract                 extracts crosslink, insertion or deletion sites
    
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
    slidingWindowToDEXSeq  transform sliding window counts to DEXSeq format
    
[In development]
    genomeToReads           splits up an genome fasta to reads for mappability tests 
    
[General help]
    -h, --help               help
    --version                version
    
'''
    
def usage_extract():
    print '''
htseq-clip extract:  extracts crosslink, insertion or deletion sites
usage:               htseq-clip extract [options]

Options:

 -i, --input                      input file (.bam)
 
 -o, --output                     output file (.bed)
 
 -c, --choice                     parameter for the choice of the cross-link sites:

                                  s for start sites. 
                                  In addition you can write e.g 
                                  s-1 then you will get the position 1 before the 
                                  start site, this is only available for the extraction 
                                  of the start sites. You can also write s-1i that means 
                                  if a postion should be below 0 it ignores it to write 
                                  out otherwise you will get an exception.
                                  
                                  m for middle sites
                                  
                                  e for end sites
                                  
                                  i for insertion sites
                                  
                                  d for deletion sites

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

 -g, --gtf         GTF or GFF file for annotation processing (.gtf)
 
 -o, --output      output file (.bed)
 
 -t, --type        Gene type for annotation

 -r, --region      True if you want exons to be split into cds and utr regions. Default is False.
  
                               
 -h, --help        help
 --version         version
'''
    
def usage_count():
    print '''
htseq-clip count:  counts the number of crosslink/deletion/insertion sites 
usage:             htseq-clip count [options]

Options:

 -i, --input         extracted crosslink, insertion or deletion sites (.bed)
 
 -o, --output        output file (.txt)
 
 -f, --compare       flattened annotation file (.bed)
 
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

 -i, --input         input file (.bed, .txt)
 
 -o, --output        output file (.html)
 
 -c, --choice        parameter for the choice of the desired plot:
 
                     c for count information plots
                     
                     j for junction information plots
                     
                     r read information plots

                              
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
'''            
def countSlidingWindow(args):
    
    bedC = bedCLIP(args)
    bedC.countSlidingWindow()

'''
Extract cross-link sites
'''
def extract(parser,args):
    
    bamC = bamCLIP(args)
    
    if hasattr(args, 'choice') and args.choice != None:
        if args.choice.startswith("s"):
            bamC.extract_StartSites() 
        elif args.choice == 'm':      
            bamC.extract_MiddleSites()
        elif args.choice == 'e':      
            bamC.extract_EndSites()
        elif args.choice == 'd':      
            bamC.extract_DeletionSites()
        elif args.choice == 'i':      
            bamC.extract_InsertionSites()
        else:
            parser.error('Invalid option for extraction')
    else:
        parser.error('You need -c option for the correct extraction of your data!')
        
'''
Processing genome.fa file into fasta format 
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

    gtf = gtfClip(args)
    gtfC = gtfCLIP(args)
    if args.region:
        gtf.processGTF()
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

def main():
    
    try:
        # get parser and add arguments
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage=globals()['__doc__'], version="%prog 1.0", conflict_handler='resolve', argument_default=argparse.SUPPRESS)

        parser.add_argument('-v', '--verbose', action='store_true', default=argparse.SUPPRESS, help='verbose output')
        parser.add_argument('-i', '--input', action='store', type= str,  default=argparse.SUPPRESS, dest='input', help='input file')
        parser.add_argument('-o', '--output', action='store', type= str,  default=argparse.SUPPRESS, dest='output', help='output file name')
        parser.add_argument('-f', '--compare', action='store', type= str,  default=argparse.SUPPRESS, dest='compare', help='file which you want to compare with your input file')
        parser.add_argument('-c', '--choice', action='store', type=str,  default=argparse.SUPPRESS, dest='choice', help='option')
        parser.add_argument('-q', '--minAlignmentQuality',   action='store', type=int, default=argparse.SUPPRESS,   dest='minAlignmentQuality', help='minimum alignment quality')
        parser.add_argument('-m', '--minReadLength', action='store', type=int, default=argparse.SUPPRESS, dest='minReadLength', help='minimum read length')
        parser.add_argument('-x', '--maxReadLength', action='store', type=int, default=argparse.SUPPRESS, dest='maxReadLength', help='maximum read length')
        parser.add_argument('-l', '--maxReadIntervalLength', action='store', type=int, default=argparse.SUPPRESS, dest='maxReadIntervalLength', help='maximum read interval length')
        parser.add_argument('-p', '--primary', action='store_true', default=argparse.SUPPRESS,  dest='primary', help='set if you want to use only primary positions of multimapping reads')
        parser.add_argument('-d', '--dist', action='store', type=int, default=argparse.SUPPRESS, dest='dist', help='Maximum distance between two sites of a read if its higher than this value it will be count as this value')
        parser.add_argument('-g', '--gtf', action='store', type=str, default=argparse.SUPPRESS, dest='gtf', help='gtf file for annotation')
        parser.add_argument('-t', '--type', action='store', type=str, default=argparse.SUPPRESS, dest='type', help='gene type for annotation')
        parser.add_argument('-n', '--name', action='store', type=str, default=argparse.SUPPRESS, dest='name', help='gene name for annotation')
        parser.add_argument('-w', '--windowSize', action='store', type=int, default=argparse.SUPPRESS, dest='windowSize', help='window size for sliding window')
        parser.add_argument('-s', '--windowStep', action='store', type=int, default=argparse.SUPPRESS, dest='windowStep', help='window step for sliding window')
        parser.add_argument('-r', '--region', action='store_true', default= argparse.SUPPRESS, dest='region', help='set if you want exons to be split into cds and utr regions. Default is False.')
        parser.add_argument('command', help='name of program to run ')
        args= parser.parse_args()
        d = vars(args)

        # check if there are all arguments
        if len(d) < 1:
            usage()
            os._exit(1)
        else:
            print len(vars(args))
            print args
            program = args.command

            if program == "extract":
                if len(d) < 3:
                    usage_extract()
                    os._exit(1)
                else:   
                    checkFileExists(args.input, parser)
                    extract(parser, args)
            elif program == 'countSlidingWindows':
                if len(d) < 3:
                    usage_countSlidingWindows()
                    os._exit(1)
                else:
                    checkFileExists(args.input, parser)
                    checkFileExists(args.compare, parser)
                    countSlidingWindow(args)
            elif program == 'junction':
                if len(d) < 3:
                    usage_junction()
                    os._exit(1)
                else: 
                    checkFileExists(args.input, parser)
                    checkFileExists(args.compare, parser)
                    junction(program, args)
            elif program == 'count':
                if len(d) < 3:
                    usage_count()
                    os._exit(1)
                else: 
                    checkFileExists(args.input, parser)
                    checkFileExists(args.compare, parser)
                    count(program, parser, args)
            elif program == 'annotation':
                if len(d) < 3:
                    usage_annotation()
                    os._exit(1)
                else:
                    checkFileExists(args.gtf, parser)
                    process(args)
            elif program == 'plot':
                if len(d) < 3:
                    usage_plot()
                    os._exit(1)
                else: 
                    checkFileExists(args.input, parser)
                    plot(parser, args)
            elif program =='createSlidingWindows':
                if len(d) < 3:
                    usage_createSlidingWindows()
                    os._exit(1)
                else: 
                    checkFileExists(args.input, parser)
                    slidingWindow(args)
            elif program =='slidingWindowsToDEXSeq':
                if len(d) < 3:
                    usage_slidingWindowsToDEXSeq()
                    os._exit(1)
                else:
                    checkFileExists(args.input, parser)
                    toDEXSeq(args)
            elif program =='genomeToReads':
                checkFileExists(args.input, parser)
                genomeToReads(args)
            else:
                parser.error ('Incorrect argument for execution') 

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
    


