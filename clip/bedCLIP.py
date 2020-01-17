# --------------------------------------------------
# bedCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf
# Institution: EMBL Heidelberg
# Modified by Sudeep Sahadevan, sahadeva@embl.de
# Date: October 2015
# --------------------------------------------------

from builtins import str
from builtins import range
from builtins import object
import gzip
import sys
import numpy as np
from collections import Counter, defaultdict

import HTSeq
from output import Output

class bedCLIP(object):
    
    data = {}
    fInput = ""
    fOutput = ""
    fCompare = ""
    choice = ""
    dist = 4000
    
    def __init__(self, options):
        
        self.fInput = options.input
        self.fOutput = options.output
        self.output = Output(self.fOutput)
        if hasattr(options,'compare'):
            self.fCompare = options.compare
        if hasattr(options,'choice'):
            self.choice = options.choice
        # self.dist = options.dist                           
        # self.data = {'dist': self.dist}
        # region types and encoding letters
        self.rtypes = {'exon':'E','intron':'I','CDS':'CDS','3UTR':'3U','5UTR':'5U'}
         
    def buildDictForComparison(self, almnt_file):
        '''
        This method builds up a dictionary for comparison analysis
        The Dictionary looks like: { chromosome : { strand : [(Start postion, end postion, name, alignment score), (...), ...] }}
        On Assumption that the read name is unique
        @TODO: needs revision, most likely a bottleneck
        '''
        d = defaultdict(dict)
        for almnt in almnt_file:
            if almnt.iv.strand == '+':
                almntInfo = [almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)]
            elif almnt.iv.strand == '-':
                almntInfo = [almnt.iv.end_d, almnt.iv.start_d, almnt.name, int(almnt.score)]
            else:
                raise ValueError('Missing strand information, strand column must be "+" or "-", found: {}'.format(almnt.iv.strand))
            try:
                d[almnt.iv.chrom][almnt.iv.strand].append(almntInfo)
            except KeyError:
                d[almnt.iv.chrom][almnt.iv.strand] = [almntInfo]
        return d
    #===================================================================================
            
    #===================================================================================
    '''
    This method calculates all the counts of cross-link sites in the 
    given reference
    ''' 
    def count_all(self):
        #Get the information for normalisation of the plots  
        if self.fCompare.endswith(".gz"):
            f = gzip.open(self.fCompare, 'r') 
        else:        
            f = open(self.fCompare, 'r')
        fn = f.readlines()
        seq = ('Chromosome','Region start pos','Region end pos','Gene ID','Gene name','Flag','Strand','Type of region','Number of exon or intron','Total exons or introns',
                'Functional type','Length in nt', 'Total cross-link sites in region','Positions where crosslinks are located', 'Max height',
                'Density','Total before duplication removal','Max height before duplication removal ')
        self.output.write("\t".join(seq)+'\n')
        for line in fn:
            if line.startswith("track"):
                self.output.write('#'+line)
        self.output.write('\n')
        f.close()
        
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
        
        d1 = self.buildDictForComparison(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)      
        
        for chrom in d2:
            if chrom not in d1:
                if self.choice is None:
                    continue
                for strand, B in d2[chrom].items():
                    # if the input file does not contain reads on the current chromosome, then write out all positions with zero
                    self._writeZeroCount(chrom,strand,B)
            for strand, B in d2[chrom].items():
                if strand not in d1[chrom]:
                    #if the input file contains reads on the current chromosome but not on the same strand, the write out all positons with zero
                    self._writeZeroCount(chrom,strand,B)
                A = d1[chrom][strand]
                self.calculateCount(A, B, chrom, strand)                                   
        self.output.close()

    def _writeZeroCount(self,chrom,strand,iv):
        '''
        Helper function: for intervals in annotion file, write those with zero reads in the input count file
        Arguments:
         chrom: chromosome name
         strand: strand info
         iv: list of intervals
        '''
        for b in iv:
            length = b[1] - b[0]
            name = b[2].split("@")
            posi = name[4].split("/")
            seq = (chrom, str(b[0]+1), str(b[1]+1), name[0], name[1], str(1), strand, name[3], posi[0], posi[1], name[2], str(length), str(0), str(0), str(0), str(0), str(0), str(0))
            self.output.write("\t".join(seq) + "\n")

    #===================================================================================
    #===================================================================================
    '''
    This method calculates only the counts of cross-link sites
    if there are counts in a region of the reference annotation
    ''' 
    def count_only(self):

        #Get the information for normalisation of the plots  
        if self.fCompare.endswith(".gz"):
            f = gzip.open(self.fCompare, 'r') 
        else:        
            f = open(self.fCompare, 'r')
        fn = f.readlines()
        seq = ('Chromosome','Region start pos','Region end pos','Gene ID','Gene name','Flag','Strand','Type of region','Number of exon or intron','Total exons or introns',
                'Functional type','Length in nt', 'Total cross-link sites in region','Positions where crosslinks are located', 'Max height',
                'Density','Total before duplication removal','Max height before duplication removal ')
        self.output.write("\t".join(seq)+'\n')
        for line in fn:
            if line.startswith("track"):
                self.output.write('#'+line)
        self.output.write('\n')
        f.close()
        
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
             
        d1 = self.buildDictForComparison(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)     
        # count reads in the commmon chromosome, common strand
        for chrom in d1:  
            if chrom not in d2:
                continue
            for strand in d1[chrom]:          
                if strand not in d2[chrom]:       
                    continue
                A = d1[chrom][strand]
                B = d2[chrom][strand]
                self.calculateCount(A, B, chrom, strand)
        self.output.close()     
    #===================================================================================
    #===================================================================================
    '''
    This method calculates the counts of cross-link sites
    '''    
    def calculateCount(self, A, B, chrom, strand):
        
        if len(B) > 0:
            
            bi = 0
            
            #boolean to stop writing out if there is no other cross-linking site
            finished = False
            #first region
            b_First = B[0]
            #last region
            b_Last  = B[-1]
            #current region
            b_Curr = B[bi]
            
            intergenicCounts = 0
            
            #length of feature
            length = 0
            #data structure for analysis
            d_count = {}
            d_dup = {}
                 
            for a in A:
                
                check = True
                
                while check:
                    #if smaller then intergenic region before first region
                    #if bigger intergenic region after last region
                    #else its in the regions         
                    if a[1] < b_First[0]:
                        intergenicCounts += 1
                        check = False
                        continue
                    elif a[0] > b_Last[1]:
                        intergenicCounts += 1
                        check = False
                    else:
                        
                        b_Curr = B[bi]
                        
                        #if bigger search go for the next region and write out
                        #the last region if zero nothing is found if not zero there are positions found
                        #else count in the current region where the cross-link site is in the positions
                        if a[0] > b_Curr[1] or finished == True:
                            
                            length = b_Curr[1] - b_Curr[0] 
                            
                            if not self.choice == "o":
                                self.writeOut(chrom, strand, b_Curr, d_count, d_dup, length)
                            elif self.choice == "o" and len(d_count) != 0:
                                self.writeOut(chrom, strand, b_Curr, d_count, d_dup, length)

                            
                            if not b_Curr == b_Last:
                                bi = bi + 1
                            else:
                                check = False
                            d_count = {}
                            d_dup = {}             
                        else:
                                     
                            if a[0] >= b_Curr[0] and a[1] <= b_Curr[1]:
                                #if nothing is in the feature generate the data structure
                                if len(d_count) == 0:
                                    length = b_Curr[1] - b_Curr[0] 
             
                                    for i in range(length+1):
                                        d_count[i+b_Curr[0]] = 0
                                        d_dup[i+b_Curr[0]] = 0
                                                                                                
                                d_count[a[0]] += 1
                                d_dup[a[0]] += a[3] 
                                
                                #if last cross-link site found write out
                                if a == A[-1]:
                                    finished = True
                                    
                                    length = b_Curr[1] - b_Curr[0] 
                                    
                                    self.writeOut(chrom, strand, b_Curr, d_count, d_dup, length)
                            
                                    if not b_Curr == b_Last:
                                        bi = bi + 1
                                    else:
                                        check = False

                                    d_count = {}
                                    d_dup = {}               
                                else:
                                    check = False
                            else:
                                intergenicCounts += 1
                                check = False
                                
            if not intergenicCounts == 0:
                seq = (chrom, '~', '~', '~', '~','~', strand, 'intergenic', '~', '~', 'intergenic', '~', str(intergenicCounts), '~', '~', '~', '~', '~')
                self.output.write("\t".join(seq) + "\n")
                                
                    
    #===================================================================================
     
    #===================================================================================

    def junction(self):
        '''
        Method that calculates the distances from cross-link sites to exon/intron regions
        '''
        d1 = self.buildDictForComparison(HTSeq.BED_Reader(self.fInput))
        d2 = self.buildDictForComparison(HTSeq.BED_Reader(self.fCompare))
        for chrom in d1:
            if chrom not in d2:
                continue
            for strand in d1[chrom]:
                if strand not in d2[chrom]:
                    continue
                A = d1[chrom][strand]
                B = d2[chrom][strand]
                if (len(A)==0) or (len(B)==0):
                    continue
                self.calculateJunction(A, B, chrom, strand)
        self.output.close()
    #===================================================================================
    #=================================================================================== 
    '''
    Calculate the distances to the junction
    '''
    def calculateJunction(self, A, B, chrom, strand):
        #first exon position in chromosome
        b_First = B[0]
        #last exon position in chromosome
        b_Last = B[-1]
        bi = 0
        #for each cl
        for a in A:
            check = True
            while check:
                #if smaller than first its intergenic region before first gene
                #else if bigger than its intergenic region after last gene
                #else its in a region or in an intergenic region between two genes
                if a[1] < b_First[0]:
                    seq = (chrom, str(a[0]), str(a[1]), '~', '~', '~', '~', strand, '~', 'intergenic', 'intergenic', str(1))
                    self.output.write("\t".join(seq) + "\n")
                    check = False
                elif a[0] > b_Last[1]:
                    seq = (chrom, str(a[0]), str(a[1]), '~', '~', '~', '~', strand, '~', 'intergenic', 'intergenic', str(3))
                    self.output.write("\t".join(seq) + "\n")
                    check = False
                else:
                        
                    #current exon/intron position
                    b_Curr = B[bi]
                    #name of current position
                    bn = b_Curr[2].split('@')
                    
                    flag = b_Curr[3]
    
                    #if bigger than count until the region where it is is found
                    if a[0] > b_Curr[1]:
                        bi = bi + 1
                    else:
                    #if between calculate the distance to the current feature
                    #else its in an intergenic region between 2 genes
                        if a[0] >= b_Curr[0] and a[1] <= b_Curr[1]:
                            
                            d1 = a[0] - b_Curr[0]
                            d2 = a[1] - b_Curr[1]
                            if strand == '-':
                                d1 = d1 * -1
                                d2 = d2 * -1
                            # modified from the original module to account for 'gene_name' in the annotation tags
                            seq = (chrom, str(a[0]), str(a[1]), bn[0], str(d2), str(d1), str(flag), strand, bn[1], bn[2], bn[3], bn[4])
                            self.output.write("\t".join(seq) + "\n")
                            check = False
                            # else:         
                            #     seq = (chrom, str(a[0]), str(a[1]), bn[0], str(d1) , str(d2), str(flag), strand, bn[1], bn[2], bn[3], bn[4])
                            #     output.write("\t".join(seq) + "\n")
                            #     check = False
                        else:
                            seq = (chrom, str(a[0]), str(a[1]), '~', '~', '~', '~', strand, '~', 'intergenic', 'intergenic', str(2))
                            self.output.write("\t".join(seq) + "\n")
                            check = False                              
    #===================================================================================
    
    def _countCrosslinks(self,almnt_file):
        '''
        test function for counting crosslink sites
        '''
        clMap = defaultdict(dict)
        d = defaultdict(dict)
        for almnt in almnt_file:
            if almnt.iv.strand == '+':
                almntInfo = almnt.iv.end_d
            elif almnt.iv.strand == '-':
                almntInfo = almnt.iv.start_d
            else:
                raise ValueError('Missing strand information, strand column must be "+" or "-", found: {}'.format(almnt.iv.strand))
            try:
                d[almnt.iv.chrom][almnt.iv.strand].append(almntInfo)
            except KeyError:
                d[almnt.iv.chrom][almnt.iv.strand] = [almntInfo]
        for chrom, strandDict in d.items():
            for strand, pos in strandDict.items():
                countPos = Counter(pos)
                countarray = np.zeros(shape=max(countPos.keys()))
                for pos, clCount in countPos.items():
                    countarray[pos] = clCount
                    clMap[chrom][strand] = countarray
        return clMap
    
    #===================================================================================
    
    def countSlidingWindow(self):
        '''
        This method is used counting the sliding window counts
        @TODO: needs improvement
        '''
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
        d1 = self._countCrosslinks(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)      
        for chrom in d1:
            if chrom not in d2:
                continue
            for strand in d1[chrom]:
                if strand not in d2[chrom]:
                    continue
                
                A = d1[chrom][strand]
                B = d2[chrom][strand]
                if len(A)==0 or len(B)==0:
                    continue
                # self.countSW(A, B, chrom, strand)
                # testing sliding window count function      
                # print("{} {} found".format(chrom,strand))
                self.countSW(A,B,chrom,strand)  
        self.output.close()
    #===================================================================================
    #===================================================================================
    '''
    This method calculates the counts and density of cross-link sites
    in a given window of an exon/intron
    ''' 
    def countSW(self, A, B, chrom, strand):
               
        #first window
        b_First = B[0]
        #last window
        b_Last  = B[-1]
        
        #data structure for analysis
        d_count = {}
        d_dup = {}
                
        ai = 0
        stepCount = 0         
                
        for b in B:  
                
            check = True
            
            #length of sliding window
            length = b[1]-b[0]
                        
            while check:
                
                #current CL
                if not ai > len(A)-1:
                    a_curr = A[ai]
                
                #if smaller then intergenic region before first region
                #if bigger intergenic region after last region
                #else its in the regions         
                if a_curr[1] < b_First[0]:            
                    check = False
                elif a_curr[0] > b_Last[1]:
                    check = False
                else:              
                    
                    #if bigger search go for the next region and write out
                    #the last region if zero nothing is found if not zero there are positions found
                    #else count in the current region where the cross-link site lies
                    if a_curr[0] > b[1] or a_curr == A[-1]:
                        
                        self.writeOut(chrom, strand, b, d_count, d_dup, length)
                        
                        ai = ai - stepCount
                        stepCount = 0
                        d_count = {}
                        d_dup = {}
                        check = False
                                
                    else:
                                
                        if a_curr[0] >= b[0] and a_curr[1] <= b[1]:
                            #if nothing is counted in the feature yet, generate the data structures
                            if len(d_count) == 0:
            
                                for i in range(length+1):
                                    d_count[i+b[0]] = 0
                                    d_dup[i+b[0]] = 0
                                                                                            
                            d_count[a_curr[0]] += 1
                            d_dup[a_curr[0]] += a_curr[3]
                            
                            ai = ai + 1 
                            stepCount = stepCount + 1                                       
                        else:
                            ai = ai + 1
                                    
    #=================================================================================== 
    #===================================================================================
    # def _countSW_numpy(self,A,B,chrom,strand):
    #     for b in B:
    #         if b[0] > A.shape[0]:
    #             # window start pos bigger than max pos in the crosslink site array
    #             continue
    #         if b[1]+1 > A.shape[0]:
    #             crosslinks = B[b[0]+1:len(A)]
    #         else:
    #             pass


    def _countSW(self, A, B, chrom, strand):
        '''
        testing sliding window count function
        Arguments:
         A: list of input positions
         B: list of window positions
         chrom: chromosome name
         strand: strand name
        '''
        npA = np.core.records.fromarrays(np.array(A).transpose(),names='begin, end, name, score',formats='i8, i8, S70, i8')
        multiMappers = set()
        for b in B:
            countInd = np.intersect1d( np.where(npA['end']>b[0]),np.where(npA['end']<=b[1]) )
            if countInd.shape[0]==0:
                continue
            readPosMap = {} # read to position map
            posReadMap = {} # position to read map
            for ci in countInd:
                if npA[ci]['name'] in multiMappers:
                    continue
                try:
                    readPosMap[npA[ci]['name']].add(npA[ci]['end'])
                except KeyError:
                    readPosMap[npA[ci]['name']] = set([npA[ci]['end']])
                try:
                    posReadMap[npA[ci]['end']].add(npA[ci]['name'])
                except KeyError:
                    posReadMap[npA[ci]['end']] = set([npA[ci]['name']])
            for readId, pos in readPosMap.items():
                # collect read id of all multimappers in this strand
                if len(pos)==1:
                    continue
                # multiMappers.add(readId)
            del readPosMap # save some space
            crosslinkCount = 0 # crosslink count
            for pos in posReadMap.keys():
                # iterate through all positions and remove multimappers
                readIds = posReadMap[pos] - multiMappers
                crosslinkCount += len(readIds)
                posReadMap[pos] = readIds
            density = float(crosslinkCount)/float(b[1]-b[0])
            clMax = max([ len(reads) for reads in list(posReadMap.values())]) # max number of crosslink sites in one pos
            dupCount = npA[countInd]['name'].shape[0] - crosslinkCount
            self._outWriter(chrom=chrom,strand=strand,windowData=b,crosslinkCount=crosslinkCount,crosslinkPosCount=len(posReadMap),maxPosCount=clMax,density=density,
                dupCount=dupCount,maxDupCountPos=0)
            # self.writeOut(chrom, strand, b, countInd.shape[0], len(readIdCount), b[1]-b[0])
    
    def _outWriter(self,chrom,strand,windowData,crosslinkCount,crosslinkPosCount,maxPosCount,density,dupCount,maxDupCountPos):
        '''
        Helper function, write outputs
        '''
        names = windowData[2].split('@')
        featInd,featCount = names[4].split('/')
        wlen = windowData[1]-windowData[0]
        outDat = [chrom,str(windowData[0]+1),str(windowData[1]+1),names[0],names[1],str(windowData[3]),strand,names[3],featInd,featCount,names[2],str(wlen),
            str(crosslinkCount),str(crosslinkPosCount),str(maxPosCount),str(density),str(dupCount),str(maxDupCountPos)]
        self.output.write("\t".join(outDat)+"\n")
    '''
    This functions converts the sliding window counts into DEXSeq format
    '''
    def toDEXSeq(self):
        
        if self.fInput.endswith(".gz"):
            almnt_file = gzip.open(self.fInput, 'r') 
        else:        
            almnt_file = open(self.fInput, 'r')
        
        for line in almnt_file:
            line = line.split('\t')
            idx = line[3]
            feature = line[7] # changed from line[6] to line[7]
            featureNr = line[8].zfill(3) # changed from line[7] to line[8]
            windowNr = line[4].zfill(4) # changed from line[9] to line[4]
            counts = line[12] # changed from line[11] to line[12]            
            if feature in self.rtypes:
                letter = self.rtypes[feature]
            else:
                # adhoc warning, should be changed to proper logging module
                sys.stderr.write('WARNING! Skipping {}, found uknown region type: {}\n'.format(line.strip('\n'),feature))
                continue
            seq = (idx+":"+letter+featureNr+"W"+windowNr, counts)
            self.output.write("\t".join(seq) + "\n")
        self.output.close()
    #===================================================================================
    
    #===================================================================================
    '''
    Write in output file
    '''
    def writeOut(self, chrom, strand, b, d_count, d_dup, length):

        name = b[2].split("@")


        if len(d_count) > 0:

            counts = 0
            for key in d_count:
                if not d_count[key] == 0:
                    counts += 1
            m_count = max(d_count.keys(), key=(lambda k: d_count[k]))

            dup_counts = 0
            for dup in d_dup:
                if not d_dup[dup] == 0:
                    dup_counts += d_dup[dup]
            m_dup = max(d_dup.keys(), key=(lambda k: d_dup[k]))

            density = float(counts) / float(length)
            
            posi = name[4].split("/")
            seq = (chrom, str(b[0]+1), str(b[1]+1), name[0],name[1], str(b[3]), strand, name[3], posi[0], posi[1], name[2], str(length), str(sum(d_count.values())), str(counts), str(d_count[m_count]), str(density), str(dup_counts), str(d_dup[m_dup]))
            self.output.write("\t".join(seq) + "\n")
        else:
            posi = name[4].split("/")
            seq = (chrom, str(b[0]+1), str(b[1]+1), name[0],name[1], str(b[3]), strand, name[3], posi[0], posi[1], name[2], str(length), str(0), str(0), str(0), str(0), str(0), str(0))
            self.output.write("\t".join(seq) + "\n")

    #===================================================================================
    
    
    #===================================================================================
    '''
    Function that calculates the minimum distance between two sites
    '''
    def compare(self, A, B, chrom, strand, out):
          
        b = 0
        bi = 0
        
        if len(B) > 0:
            
            dTMP = None
            
            #foreach site in a
            for a in A:
                
                check = True
                            
                while check:
                    
                    #if bigger than no bigger position than a was found
                    #and the distance between the last one will be used
                    if bi > len(B)-1:
                        
                        b = B[bi-1]
                        if a[0] > b[1]:
                            d = a[0] - b[1] + 1
                            if strand == '-':
                                d = d * (-1)
                            out.write(str(d) + '\n') 
                        else:
                            if a[0] == b[1]:
                                d = 1
                                if strand == '-':
                                    d = d * (-1)
                                out.write(str(d) + '\n')     
                            else: 
                                out.write(str(0) + '\n')                              
                        check = False
                    else:
                        b = B[bi]
                        
                        #If the two sites descend from the same read skip
                        if a[2] == b[2]:
                            bi += 1
                            continue
                        
                        #if smaller than site, memorize the distance and go for the next one
                        if b[0] < a[0] and b[1] < a[0]:
                            dTMP = a[0] - b[1] + 1
                            bi += 1
                        else:
                            #if 2 reads overlap distance is zero
                            if b[0] < a[0] and b[1] > a[0] and a[0] != b[1]:
                                out.write(str(0) + '\n')
                            #if 2 reads overlap only(!) in a[start] and b[end] position distance is 1
                            elif b[0] < a[0] and b[1] == a[0]:
                                d = 1
                                if strand == '-':
                                    d = d * (-1)
                                out.write(str(d) + '\n') 
                            else:
                                #if bigger compare the distances between the smaller one and bigger one
                                #and the minimum distance will be written out
                                #else if the read lies in the other read or the overlap so distance is zero
                                #else if 2 reads overlap only(!) in a[end] and b[start] position distance is 1
                                if b[0] > a[1]:
                                    d = b[0] - a[1] + 1
                                    if d < dTMP or dTMP == None:
                                        if strand == '+':
                                            d = d * (-1)
                                        out.write(str(d) + '\n') 
                                    else:
                                        if strand == '-':
                                            dTMP = dTMP * (-1)
                                        out.write(str(dTMP) + '\n')     
                                else:
                                    if a[1] == b[0]:
                                        d = 1
                                        if strand == '+':
                                            d = d * (-1)
                                        out.write(str(d) + '\n') 
                                    else:
                                        out.write(str(0) + '\n') 
                                    
                                if not bi == 0:
                                    bi = bi - 1
                            
                            check = False                                     
    #===================================================================================
    #================================================================================= 
    '''
    This method calculates the minimum distances between the sites in the same file
    '''
    def calcDistofSite(self, almnt_file):
        
        d = {}
        
        for almnt in almnt_file:
            if almnt.iv.chrom not in d:
                d[almnt.iv.chrom] = {almnt.iv.strand : [almnt.iv.start_d]}
            else:
                if almnt.iv.strand not in d[almnt.iv.chrom]:
                    d[almnt.iv.chrom][almnt.iv.strand] = [almnt.iv.start_d]
                else:
                    if almnt.iv.strand == '+':
                        if  almnt.iv.start_d - d[almnt.iv.chrom][almnt.iv.strand][-1] > self.data['Dist']:
                            self.output.write(str(self.data['Dist']) + "\n")
                        else:
                            self.output.write(str(almnt.iv.start_d - d[almnt.iv.chrom][almnt.iv.strand][-1]) + "\n")
                        d[almnt.iv.chrom][almnt.iv.strand].append(almnt.iv.start_d)
                    elif almnt.iv.strand == '-':
                        if d[almnt.iv.chrom][almnt.iv.strand][-1] - almnt.iv.start_d < self.data['Dist']*(-1):
                            self.output.write(str(self.data['Dist']*(-1)) + "\n")
                        else:
                            self.output.write(str(d[almnt.iv.chrom][almnt.iv.strand][-1] - almnt.iv.start_d) + "\n")
                        d[almnt.iv.chrom][almnt.iv.strand].append(almnt.iv.start_d)
                        
    #===================================================================================         
    #===================================================================================
    '''
    This method calculates the distances between each cross-link site on the same strand on
    the same chromosome
    '''
    def calcDistancesFromSite(self):
        
        almnt_file = HTSeq.BED_Reader(self.fInput)   
        self.calcDistofSite(almnt_file)

        self.output.close()
      
    #===================================================================================
