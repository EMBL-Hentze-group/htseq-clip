# --------------------------------------------------
# htseq-clip main
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# --------------------------------------------------

'''
    SYNOPSIS

        clip.py execution [options]

    DESCRIPTION
    
        --execution
                    Functionality you want to execute; Available functions:
                          bokeh
                          count
                          countSlidingWindow
                          extract
                          genomeToReads
                          junction    
                          process
                          slidingWindow
                          toDexSeq      
                    
        [-i, --input=FILE]
                          Input file

        [-o, --output=FILE]
                          Output file
                          
        [-f, --compare=FILE]
                          File which you want to compare with your input file
        
        [-c, --choice]
                          Option for different functions
                          default: s for start sites extraction (can vary in function and option of course)
        
        [-q, --minAlignmentQuality]
                          Minimum alignment quality
                          default: 10

        [-m, --minReadLength]
                          Minimum read length
                          default: 0
        
        [-x, --maxReadLength]
                          Maximum read length
                          default: 0
                          
        [-l, --maxReadIntervalLength]
                          Maximum read interval length
                          default: 100000

        [-p, --primary]
                          Set if you want to use only primary positions of multimapping reads
                          default: False
                          
        [-d, --dist]
                          Distance between two sites of a read e.g if its higher than this value
                          it will be count as this value; you can use this variable for different
                          distances matters
                          default: 4000
                          
        [-g, --gtf]
                          GTF or GFF file for annotation processing
                          
        [-t, --type]
                          Gene type for annotation
                          
        [-w, --windowSize]
                          window size for sliding window
                          
        [-s, --windowStep]
                          window step for sliding window
        
        [-h, --help]      display help message with all parameters and options


        [--verbose]       display verbose output for debugging


    EXAMPLES

        python clip.py extract -i input.bam -o output.bed -c s

        python clip.py --help

'''
    
import optparse, traceback, os
from bamCLIP import bamCLIP
from bedCLIP import bedCLIP
from gtfCLIP import gtfCLIP
from bokehCLIP import bokehCLIP
from fastaCLIP import fastaCLIP

VERSION = "0.1.0"

#=====================================================================================
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
            parser.error ('Incorrect number of user-specified arguments') 
        else:

            program = args[0]

            if program == "extract":
                checkFileExists(options.input, parser)
                extract(parser, options, args)
            elif program == 'countSlidingWindow':
                checkFileExists(options.input, parser)
                checkFileExists(options.compare, parser)
                countSlidingWindow(options, args)
            elif program == 'junction':
                checkFileExists(options.input, parser)
                checkFileExists(options.compare, parser)
                junction(program, options, args)
            elif program == 'count':
                checkFileExists(options.input, parser)
                checkFileExists(options.compare, parser)
                count(program, parser, options, args)
            elif program == 'process':
                checkFileExists(options.gtf, parser)
                process(options, args)
            elif program == 'bokeh':
                checkFileExists(options.input, parser)
                plot(parser, options, args)
            elif program =='slidingWindow':
                checkFileExists(options.input, parser)
                slidingWindow(options, args)
            elif program =='toDexSeq':
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
    


