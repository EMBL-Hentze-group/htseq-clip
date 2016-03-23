# --------------------------------------------------
# htseq-clip main
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# --------------------------------------------------
    
import optparse, traceback, os
from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from gtfCLIP import gtfCLIP
from bokehCLIP import bokehCLIP
from fastaCLIP import fastaCLIP

VERSION = "0.1.0"

'''
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

[Distances]
  junction                calculates distances to junctions

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
    
[Distances]
    junction                calculates distances to junctions
    
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

<<<<<<< HEAD
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
def count(program, parser, options, args):
    
    bedC = bedCLIP(options)
    
    if hasattr(options, 'choice') and options.choice != None:
        if options.choice == 'a':
            bedC.count_all() 
        elif options.choice == 'o':      
            bedC.count_only()
        else:
            parser.error('Invalid option for count')
    else:
        parser.error('You need -c option for the correct counting of your data!')
        
'''
Sliding window counts
'''            
def countSlidingWindow(options, args):
    
    bedC = bedCLIP(options)
    bedC.countSlidingWindow()

'''
Extract cross-link sites
'''
def extract(parser, options ,args):
    
    bamC = bamCLIP(options)
    
    if hasattr(options, 'choice') and options.choice != None:
        if options.choice.startswith("s"):
            bamC.extract_StartSites() 
        elif options.choice == 'm':      
            bamC.extract_MiddleSites()
        elif options.choice == 'e':      
            bamC.extract_EndSites()
        elif options.choice == 'd':      
            bamC.extract_DeletionSites()
        elif options.choice == 'i':      
            bamC.extract_InsertionSites()
        else:
            parser.error('Invalid option for extraction')
    else:
        parser.error('You need -c option for the correct extraction of your data!')
        
'''
Processing genome.fa file into fasta format 
'''    
def genomeToReads(options, args):
    
    fastaC = fastaCLIP(options)   
    fastaC.genomeToReads()

'''
Calculate distance from cross-link site to exon/intron junction site
'''        
def junction(program, options, args):
    
    bedC = bedCLIP(options)
    bedC.junction()

'''
Plotting function
'''        
def plot(parser, options, args):
    
    bokehC = bokehCLIP(options)
    
    if hasattr(options, 'choice') and options.choice != None:
        if options.choice == 'c':      
            bokehC.plot_count()
        elif options.choice == 'd':
            bokehC.plot_comparison() 
        elif options.choice == 'j':      
            bokehC.plot_junction()
        elif options.choice == 'r':      
            bokehC.plot_read()
        
        else:
            parser.error('Invalid option for plotting')
    else:
        parser.error('You need -c option for the correct plotting of your data!')
        
'''
Processing the annotation file
'''       
def process(options, args):
    
    gtfC = gtfCLIP(options)
    gtfC.processGTF()

'''
Processing the annotation file into sliding windows
'''    
def slidingWindow(options, args):
    
    gtfC = gtfCLIP(options)
    gtfC.slidingWindow()

'''
Converting the sliding window counts into DEXSeq format
'''    
def toDEXSeq(options, args):
    
    bedC = bedCLIP(options)
    bedC.toDEXSeq()

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
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version="%prog 1.0")

        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-i', '--input', action='store', type='string',  default="../test/dummy_test/test.bam", dest='input', help='input file ')
        parser.add_option('-o', '--output', action='store', type='string',  default="../test/dummy_test/output.bed", dest='output', help='output file name')
        parser.add_option('-f', '--compare', action='store', type='string',  default="../test/dummy_test/output_130000SS.bed", dest='compare', help='file which you want to compare with your input file')
        parser.add_option('-c', '--choice', action='store', type='string',  default=None, dest='choice', help='option')
        parser.add_option('-q', '--minAlignmentQuality',   action='store', type='int', default=10,   dest='minAlignmentQuality', help='minimum alignment quality')
        parser.add_option('-m', '--minReadLength', action='store', type='int', default=0, dest='minReadLength', help='minimum read length')
        parser.add_option('-x', '--maxReadLength', action='store', type='int', default=0, dest='maxReadLength', help='maximum read length')
        parser.add_option('-l', '--maxReadIntervalLength', action='store', type='int', default=10000, dest='maxReadIntervalLength', help='maximum read interval length')
        parser.add_option('-p', '--primary',   action='store_true', default=False,  dest='primary', help='set if you want to use only primary positions of multimapping reads')
        parser.add_option('-d', '--dist', action='store', type='int', default=4000, dest='dist', help='Maximum distance between two sites of a read if its higher than this value it will be count as this value')
        parser.add_option('-g', '--gtf', action='store', type='string', default='../test/gff_gtf/human.gtf', dest='gtf', help='gtf file for annotation')
        parser.add_option('-t', '--type', action='store', type='string', default='gene_biotype', dest='type', help='gene type for annotation')
        parser.add_option('-w', '--windowSize', action='store', type='int', default=50, dest='windowSize', help='window size for sliding window')
        parser.add_option('-s', '--windowStep', action='store', type='int', default=20, dest='windowStep', help='window step for sliding window')
            
        (options, args) = parser.parse_args()

        # check if there are all arguments
        if len(args) != 1:
            usage()
            os._exit(1)
        else:

            program = args[0]

            if program == "extract":
                if len(args) < 3:
                    usage_extract()
                    os._exit(1)
                else:   
                    checkFileExists(options.input, parser)
                    extract(parser, options, args)
            elif program == 'countSlidingWindows':
                if len(args) < 3:
                    usage_countSlidingWindows()
                    os._exit(1)
                else:
                    checkFileExists(options.input, parser)
                    checkFileExists(options.compare, parser)
                    countSlidingWindow(options, args)
            elif program == 'junction':
                if len(args) < 3:
                    usage_junction()
                    os._exit(1)
                else: 
                    checkFileExists(options.input, parser)
                    checkFileExists(options.compare, parser)
                    junction(program, options, args)
            elif program == 'count':
                if len(args) < 3:
                    usage_count()
                    os._exit(1)
                else: 
                    checkFileExists(options.input, parser)
                    checkFileExists(options.compare, parser)
                    count(program, parser, options, args)
            elif program == 'annotation':
                if len(args) < 3:
                    usage_annotation()
                    os._exit(1)
                else: 
                    checkFileExists(options.gtf, parser)
                    process(options, args)
            elif program == 'plot':
                if len(args) < 2:
                    usage_plot()
                    os._exit(1)
                else: 
                    checkFileExists(options.input, parser)
                    plot(parser, options, args)
            elif program =='createSlidingWindows':
                if len(args) < 3:
                    usage_createSlidingWindows()
                    os._exit(1)
                else: 
                    checkFileExists(options.input, parser)
                    slidingWindow(options, args)
            elif program =='slidingWindowsToDEXSeq':
                if len(args) < 3:
                    usage_slidingWindowsToDEXSeq()
                    os._exit(1)
                else:
                    checkFileExists(options.input, parser)
                    toDEXSeq(options, args)
            elif program =='genomeToReads':
                checkFileExists(options.input, parser)
                genomeToReads(options, args)
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
    


