# --------------------------------------------------
# gtfCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: November 2015
# --------------------------------------------------


import gzip, HTSeq
from collections import OrderedDict

    
class gtfCLIP:
    
    features = {}
    fInput = ''
    gtfFile = ''
    fOutput = ''
    geneType = ''
    windowSize = 50
    windowStep = 20
    
    def __init__(self, options):
        
        if hasattr(options, 'gtf'):
            self.gtfFile = options.gtf
            
        if hasattr(options, 'input'):
            self.fInput = options.input
        
        if hasattr(options, 'output'):
<<<<<<< HEAD
            self.fOutput = options.output  
=======
            self.fOutput = options.output   
>>>>>>> b1fbc1dde60f3c475e7c8d5afcb0a3c6d16eb6df
            
        if hasattr(options, 'windowSize'):
            self.windowSize = options.windowSize  
            
        if hasattr(options, 'windowStep'):
            self.windowStep = options.windowStep 
            
        if hasattr(options, 'type'):
            self.geneType = options.type 
            
            
    '''
    This method processes through the gtf file and determinates the positions of exons and introns
    '''        
    def processGTF(self):
        
        gtf = HTSeq.GFF_Reader(self.gtfFile)
                
        gas = None
        
        name = ''
        chrom = ''
        strand = ''
        t = ''
        start = 0
        end = 0
        ec = 1
        ic = 1
        ps = []
        eFlag = None
        fs = []
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w')
            
        chromFooter = {}
        typeFooter = {}
        
        #for each feature in gtf file                    
        for feature in gtf:
            
            #if gene a new region is found
            #else a new exon is found
            if feature.type == "gene" or feature.type == "tRNAscan":
                
                #if genomic array of sets (gas) is none generate a new one
                #else calculate positions of exon or intron (always starts and ends with exon
                #between are introns
                if gas == None:
                    gas = HTSeq.GenomicArrayOfSets('auto', stranded = True)
                else:
                    #Generate a list for all exons in gene
                    l = list(gas[HTSeq.GenomicInterval(chrom, start, end, strand)].steps())
                    
                    #Sort the list where all exons are in
                    ps.sort(key=lambda tup: tup[0])
                    
                    #calculate number of introns and exons
                    iLength = int(len(l)/2)
                    eLength = int(len(l)/2) + 1
                    
                    exon = True
                         
                    #Depending on if current feature is exon or intron calculate the correct position
                    #and flag the exons
                    for e in l:
                        if exon:
                            
                            while len(ps) > 0 and ps[0][0] < e[0].end:
                                if ps[0][0] == e[0].start and ps[0][1] == e[0].end:
                                    eFlag = 3
                                elif ps[0][0] == e[0].start and ps[0][1] < e[0].end and eFlag != 0:
                                    eFlag = 2
                                elif ps[0][0] > e[0].start and ps[0][1] == e[0].end and eFlag != 0:
                                    eFlag = 1
                                elif  ps[0][0] > e[0].start and ps[0][1] < e[0].end and eFlag != 0:
                                    eFlag = 0
                                del ps[0]
                            
                            #Calculate the values which are need for the normalization of the plots    
                            if t == "protein_coding":
                                if not typeFooter.has_key("protein_coding_exon"):
                                    typeFooter["protein_coding_exon"] = (e[0].end - e[0].start)
                                else:
                                    typeFooter["protein_coding_exon"] += (e[0].end - e[0].start)
                            else:    
                                if not typeFooter.has_key(t):
                                    typeFooter[t] = (e[0].end - e[0].start)
                                else:
                                    typeFooter[t] += (e[0].end - e[0].start)
                                                   
                            seq = [e[0].chrom, str(e[0].start), str(e[0].end), name+'@exon@'+str(ec)+'/'+str(eLength), eFlag, e[0].strand, "e"]
                            fs.append(seq)
                            ec += 1
                            exon = False
    
                        else:
                            
                            #Calculate the values which are need for the normalization of the plots   
                            if t == "protein_coding":
                                if not typeFooter.has_key("protein_coding_intron"):
                                    typeFooter["protein_coding_intron"] = (e[0].end - e[0].start)
                                else:
                                    typeFooter["protein_coding_intron"] += (e[0].end - e[0].start)
                            else:    
                                if not typeFooter.has_key(t):
                                    typeFooter[t] = (e[0].end - e[0].start)
                                else:
                                    typeFooter[t] += (e[0].end - e[0].start)
                                
                            seq = [e[0].chrom, str(e[0].start), str(e[0].end), name+'@intron@'+str(ic)+'/'+str(iLength), 1, e[0].strand, "i"]
                            fs.append(seq)
                            ic += 1
                            exon = True
                            
                    #Flagging of introns and writing out in output file
                    #else there is only one exon
                    if not len(fs) == 1:
                            
                        #Initialisation of Flags
                        p = 0
                        eFlag1 = fs[p][4]
                        p = p + 2
                        eFlag2 = fs[p][4]
                             
                        for line in fs:
                            if line[6] == "e":
                                seq = (line[0], line[1], line[2], line[3], str(line[4]), line[5], "\n")
                                output.write(str('\t').join(seq))
                            elif line[6] == "i":
                                
                                #Simple binary operations to calculate the flag for the introns
                                #The intron flag is calculated by using the left exon flag and the 
                                #right exon flag by using the bits of the binary coded flag
                                #for further information how they are calculated look up the
                                #documentation
                                line[4] = ((eFlag1 & 1) << 1) | (eFlag2 >> 1)    
                                
                                seq = (line[0], line[1], line[2], line[3], str(line[4]), line[5], "\n")
                                output.write(str('\t').join(seq))
                 
                                eFlag1 = eFlag2
                                p = p + 2
                                if not p > len(fs):
                                    eFlag2 = fs[p][4]
                    else:
                        #if gene type is is something else, you have to flag it with 3 because there only exists one isoform
                        if fs[0][4] == None:
                            fs[0][4] = 3
                        seq = (fs[0][0], fs[0][1], fs[0][2], fs[0][3], str(fs[0][4]), fs[0][5], "\n")
                        output.write(str('\t').join(seq))                                 
                                
                            
                    gas = HTSeq.GenomicArrayOfSets('auto', stranded = True)
                
                chrom = feature.iv.chrom
                
                #Calculate the values which are need for the normalization of the plots   
                if not chromFooter.has_key(chrom):
                    chromFooter[chrom] = (feature.iv.end - feature.iv.start)
                else:
                    chromFooter[chrom] += (feature.iv.end - feature.iv.start)
                    
                strand = feature.iv.strand
                name = feature.name
                if feature.attr.has_key(self.geneType):
                    attribute = str(feature.attr[str(self.geneType)])
                else:
                    error = "Wrong gene type: "+self.geneType+". Check your annotation file!!"
                    raise KeyError(error)
                         
                name = name + '@' + attribute
                t = attribute
                start = feature.iv.start
                end = feature.iv.end
                ec = 1
                ic = 1
                ps = []
                fs = []
                                
            elif feature.type == "exon" or feature.type == "exonic":
                giv = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.strand)
                ps.append((feature.iv.start, feature.iv.end))
                s = 'exon'
                gas[giv] += s    
        
        #Write out information for normalization of plots
        chromFooter = OrderedDict(sorted(chromFooter.items(), key=lambda x: (-x[1], x[0])))  
        
        for k in chromFooter:
            output.write("track chr "+str(k)+" "+str(chromFooter[k])+"\n")
        
        typeFooter = OrderedDict(sorted(typeFooter.items(), key=lambda x: (-x[1], x[0])))  
        
        for k in typeFooter:
            output.write("track type "+str(k)+" "+str(typeFooter[k])+"\n")
        
        output.close()                 
    #================================================================================= 
       
    #=================================================================================
    '''
    This functions calculates the sliding window positions
    '''
    def slidingWindow(self):
          
        if self.fInput.endswith(".gz"):
            almnt_file = gzip.open(self.fInput, 'r') 
        else:        
            almnt_file = open(self.fInput, 'r')
        
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w') 
        else:        
            output = open(self.fOutput, 'w')
        
        currentName = None
        
        for line in almnt_file:
            if line.startswith("track"):
                continue
            
            line = line.split('\n')
            line = line[0].split('\t')
            
            name = line[3].split('@')
            
            if currentName == None or currentName != name[0]:
                currentName = name[0]
                windowCount = 1
            
            strand = line[5]
            
            start = int(line[1])
            end = int(line[2])
            
            pos1 = start
            pos2 = start + self.windowSize
                  
            #if length shorter than given windowsize then the whole feature is one window
            #else split the current feature up by the given options
            if pos2 >= end:
                seq = (line[0], str(start), str(end), name[0]+"@"+str(windowCount)+"@"+name[2]+"@"+name[3], line[4], strand)     
                output.write(str('\t').join(seq) + "\n")
                
                windowCount = windowCount + 1
            else:
                while pos2 < end:

                    seq = (line[0], str(pos1), str(pos2), name[0]+"@"+str(windowCount)+"@"+name[2]+"@"+name[3], line[4], strand)
                    output.write(str('\t').join(seq) + "\n")
                    
                    pos1 = pos1 + self.windowStep
                    pos2 = pos2 + self.windowStep
                    
                    windowCount = windowCount + 1
                          
                    if pos2 > end:
                        pos2 = end
                        seq = (line[0], str(pos1), str(pos2), name[0]+"@"+str(windowCount)+"@"+name[2]+"@"+name[3], line[4], strand)
                        output.write(str('\t').join(seq) + "\n")
                        
                        windowCount = windowCount + 1
                
        output.close()          
    #==================================================================================
