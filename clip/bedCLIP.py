# --------------------------------------------------
# bedCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# --------------------------------------------------

import os, gzip

try:
    import HTSeq
except Exception:
    print "Please install the HTSeq framework e.g. like this"
    print "pip install HTSeq"
    print "pip install HTSeq --user"
    os._exit(1)
    
class bedCLIP:
    
    data = {}
    fInput = ""
    fOutput = ""
    fCompare = ""
    choice = ""
    dist = 4000
    
    def __init__(self, options):
        
        if hasattr(options, 'input'):
            self.fInput = options.input
        
        if hasattr(options, 'output'):
            self.fOutput = options.output
            
        if hasattr(options, 'compare'):
            self.fCompare = options.compare
            
        if hasattr(options, 'choice'):
            self.choice = options.choice
            
        if hasattr(options, 'dist'):
            self.dist = options.dist
                           
        self.data = {'dist': self.dist}
         
    #=================================================================================
    '''
    This method builds up a dictionary for comparison analysis
    The Dictionary looks like: { chromosome : { strand : [(Start postion, end postion, name, alignment score), (...), ...] }}
    On Assumption that the read name is unique
    '''
    def buildDictForComparison(self, almnt_file):
        
        d = {}
        
        for almnt in almnt_file:
            
            if almnt.iv.strand == '+':
                if not d.has_key(almnt.iv.chrom):
                    d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)]]}
                else:
                    if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                        d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)]]
                    else:
                        d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)]) 
            else:
                if not d.has_key(almnt.iv.chrom):
                    d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.end_d, almnt.iv.start_d, almnt.name, int(almnt.score)]]}
                else:
                    if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                        d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.end_d, almnt.iv.start_d, almnt.name, int(almnt.score)]]
                    else:
                        d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.end_d, almnt.iv.start_d, almnt.name, int(almnt.score)])        
                        
        return d       
    #===================================================================================
            
    #===================================================================================
    '''
    This method calculates all the counts of cross-link sites in the 
    given reference
    ''' 
    def count_all(self):
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w') 
        
        #Get the information for normalisation of the plots  
        if self.fCompare.endswith(".gz"):
            f = gzip.open(self.fCompare, 'r') 
        else:        
            f = open(self.fCompare, 'r') 
            
        for line in f:
            if line.startswith("track"):
                output.write(line)
                
        f.close()
        
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
        
        d1 = self.buildDictForComparison(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)      
        
        for chrom in d2:
            
            if not d1.has_key(chrom):
                if self.choice == None:
                    continue
            for strand in d2[chrom]:
                
                B = d2[chrom][strand] 
                
                #if the input file does not contain reads 
                #on the current chromosome, then write out all positions with zero
                if not d1.has_key(chrom):
                    for b in B:
                        
                        length = b[1] - b[0]
                        
                        name = b[2].split("@")
                        posi = name[3].split("/")  
                        
                        seq = (chrom, str(b[0]+1), str(b[1]+1), name[0], str(1), strand, name[2], posi[0], posi[1], name[1], str(length), str(0), str(0), str(0), str(0), str(0), str(0))
                        output.write(str("\t").join(seq) + "\n")
                    
                    continue
                
                #if the input file contains reads on the current chromosome but not on the same
                #strand, the write out all positons with zero
                elif not d1[chrom].has_key(strand):          
                    for b in B:
                        
                        length = b[1] - b[0]
                        
                        name = b[2].split("@")
                        posi = name[3].split("/")  
                        
                        seq = (chrom, str(b[0]+1), str(b[1]+1), name[0], str(1), strand, name[2], posi[0], posi[1], name[1], str(length), str(0), str(0), str(0), str(0), str(0), str(0))
                        output.write(str("\t").join(seq) + "\n")
                    continue
                         
                A = d1[chrom][strand]
      
                self.calculateCount(A, B, chrom, strand, output)
                                   
        output.close()     
    #===================================================================================
    #===================================================================================
    '''
    This method calculates only the counts of cross-link sites
    if there are counts in a region of the reference annotation
    ''' 
    def count_only(self):
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w') 
        
        #Get the information for normalisation of the plots  
        if self.fCompare.endswith(".gz"):
            f = gzip.open(self.fCompare, 'r') 
        else:        
            f = open(self.fCompare, 'r') 
            
        for line in f:
            if line.startswith("track"):
                output.write(line)
                
        f.close()
        
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
             
        d1 = self.buildDictForComparison(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)     
        
        #only if there are reads on the current chromsome
        #and on the same strand, the counting is performed
        for chrom in d1:  
            if not d2.has_key(chrom):
                continue
            for strand in d1[chrom]:          
                if not d2[chrom].has_key(strand):       
                    continue
                         
                A = d1[chrom][strand]
                B = d2[chrom][strand]
      
                self.calculateCount(A, B, chrom, strand, output)
                                   
        output.close()     
    #===================================================================================
    #===================================================================================
    '''
    This method calculates the counts of cross-link sites
    '''    
    def calculateCount(self, A, B, chrom, strand, output):
        
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
                                self.writeOut(chrom, strand, b_Curr, d_count, d_dup, length, output)
                            elif self.choice == "o" and len(d_count) != 0:
                                self.writeOut(chrom, strand, b_Curr, d_count, d_dup, length, output)
                            
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
                                    
                                    self.writeOut(chrom, strand, b_Curr, d_count, d_dup, length, output)
                            
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
                seq = (chrom, '~', '~', '~', '~', strand, 'intergenic', '~', '~', 'intergenic', '~', str(intergenicCounts), '~', '~', '~', '~', '~')
                output.write(str("\t").join(seq) + "\n")
                                
                    
    #===================================================================================
     
    #===================================================================================
    '''
    Method that calculates the distances from cross-link sites to exon/intron regions
    '''
    def junction(self):
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w')
                
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
        
        d1 = self.buildDictForComparison(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)
         
        for chrom in d1:
            if not d2.has_key(chrom):
                continue
            for strand in d1[chrom]:
                if not d2[chrom].has_key(strand):
                    continue
                
                A = d1[chrom][strand]
                B = d2[chrom][strand]
            
                self.calculateJunction(A, B, chrom, strand, output)
          
        output.close()
    #===================================================================================
    #=================================================================================== 
    '''
    Calculate the distances to the junction
    '''
    def calculateJunction(self, A, B, chrom, strand, output):
         
        if len(B) > 0:
              
            #first exon position in chromosom   
            b_First = B[0]
            #last exon postion in chromosom
            b_Last  = B[-1]
              
            bi = 0
                
            #foreach cl       
            for a in A:
                  
                check = True
                  
                while check:
                              
                    #if smaller than first its intergenic region before first gene
                    #else if bigger than its intergenic region after last gene
                    #else its in a region or in an intergenic region between two genes
                    if a[1] < b_First[0]:
                        seq = (chrom, str(a[0]), str(a[1]), '~', '~', '~', '~', strand, 'intergenic', 'intergenic', str(1))
                        output.write(str("\t").join(seq) + "\n")
                        check = False
                    elif a[0] > b_Last[1]:
                        seq = (chrom, str(a[0]), str(a[1]), '~', '~', '~', '~', strand, 'intergenic', 'intergenic', str(3))
                        output.write(str("\t").join(seq) + "\n")
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
                                    seq = (chrom, str(a[0]), str(a[1]), bn[0], str(d2), str(d1), str(flag), strand, bn[1], bn[2], bn[3])
                                    output.write(str("\t").join(seq) + "\n")
                                    check = False
                                else:         
                                    seq = (chrom, str(a[0]), str(a[1]), bn[0], str(d1) , str(d2), str(flag), strand, bn[1], bn[2], bn[3])
                                    output.write(str("\t").join(seq) + "\n")
                                    check = False
                            else:
                                seq = (chrom, str(a[0]), str(a[1]), '~', '~', '~', '~', strand, 'intergenic', 'intergenic', str(2))
                                output.write(str("\t").join(seq) + "\n")
                                check = False                              
    #===================================================================================
    
    
    #===================================================================================
    '''
    This method is used counting the sliding window counts
    '''
    def countSlidingWindow(self):
        
        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)
        
        d1 = self.buildDictForComparison(almnt_file1)
        d2 = self.buildDictForComparison(almnt_file2)      
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w')  
        
        for chrom in d1:
            if not d2.has_key(chrom):
                continue
            for strand in d1[chrom]:
                if not d2[chrom].has_key(strand):
                    continue
                
                A = d1[chrom][strand]
                B = d2[chrom][strand]  
                
                self.countSW(A, B, chrom, strand, output)
                                   
        output.close() 
    #===================================================================================
    #===================================================================================
    '''
    This method calculates the counts and density of cross-link sites
    in a given window of an exon/intron
    ''' 
    def countSW(self, A, B, chrom, strand, output):
        
        if len(A) > 0 and len(B) > 0:
               
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
                            
                            self.writeOut(chrom, strand, b, d_count, d_dup, length, output)
                            
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
    '''
    This functions converts the sliding window counts into DEXSeq format
    '''
    def toDEXSeq(self):
        
        if self.fInput.endswith(".gz"):
            almnt_file = gzip.open(self.fInput, 'r') 
        else:        
            almnt_file = open(self.fInput, 'r')
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w')
        
        for line in almnt_file:
            line = line.split('\t')
            
            idx = line[3]
            feature = line[6]
            featureNr = line[7].zfill(3)
            windowNr = line[9].zfill(4)
            counts = line[11]
            
            letter = ""
            
            if feature == "exon":
                letter = "E"
            elif feature == "intron":
                letter = "I"
            else:
                raise ValueError("Wrong feature detected! Check your data!")
                   
            seq = (idx+":"+letter+featureNr+"W"+windowNr, counts)
            output.write(str("\t").join(seq) + "\n")
            
        
        output.close()
    #===================================================================================
    
    #===================================================================================
    '''
    Write in output file
    '''
    def writeOut(self, chrom, strand, b, d_count, d_dup, length, output):
        
        name = b[2].split("@")
        posi = name[3].split("/")
                            
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
            
            print b, name
                                
            seq = (chrom, str(b[0]+1), str(b[1]+1), name[0], str(b[3]), strand, name[2], posi[0], posi[1], name[1], str(length), str(sum(d_count.values())), str(counts), str(d_count[m_count]), str(density), str(dup_counts), str(d_dup[m_dup]))
            output.write(str("\t").join(seq) + "\n")
        else:
            seq = (chrom, str(b[0]+1), str(b[1]+1), name[0], str(b[3]), strand, name[2], posi[0], posi[1], name[1], str(length), str(0), str(0), str(0), str(0), str(0), str(0))
            output.write(str("\t").join(seq) + "\n")
     
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
    def calcDistofSite(self, almnt_file, output):
        
        d = {}
        
        for almnt in almnt_file:
            if not d.has_key(almnt.iv.chrom):
                    d[almnt.iv.chrom] = {almnt.iv.strand : [almnt.iv.start_d]}
            else:
                if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                    d[almnt.iv.chrom][almnt.iv.strand] = [almnt.iv.start_d]
                else:
                    if almnt.iv.strand == '+':
                        if  almnt.iv.start_d - d[almnt.iv.chrom][almnt.iv.strand][-1] > self.data['Dist']:
                            output.write(str(self.data['Dist']) + "\n")
                        else:
                            output.write(str(almnt.iv.start_d - d[almnt.iv.chrom][almnt.iv.strand][-1]) + "\n")
                        d[almnt.iv.chrom][almnt.iv.strand].append(almnt.iv.start_d)
                    elif almnt.iv.strand == '-':
                        if d[almnt.iv.chrom][almnt.iv.strand][-1] - almnt.iv.start_d < self.data['Dist']*(-1):
                            output.write(str(self.data['Dist']*(-1)) + "\n")
                        else:
                            output.write(str(d[almnt.iv.chrom][almnt.iv.strand][-1] - almnt.iv.start_d) + "\n")
                        d[almnt.iv.chrom][almnt.iv.strand].append(almnt.iv.start_d)
                        
    #===================================================================================         
    #===================================================================================
    '''
    This method calculates the distances between each cross-link site on the same strand on
    the same chromosome
    '''
    def calcDistancesFromSite(self):
        
        almnt_file = HTSeq.BED_Reader(self.fInput)
          
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w')     
        
        self.calcDistofSite(almnt_file, output)

        output.close()
      
    #===================================================================================
      
    
    
     
    
    
    
    
    
    
    