# --------------------------------------------------
# bamCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# --------------------------------------------------

import gzip, decimal, HTSeq
from output import Output

class bamCLIP:
    
    # default parameters
    data = {}
    fInput = ""
    fOutput = ""
    maxReadLength = 0
    minReadLength = 0
    maxReadIntervalLength = 10000
    minAlignmentQuality = 10
    primary = False
    choice = ''
    mate = 1
    
    count = 0
    
    def __init__(self, options):
        
        if hasattr(options, 'input'):
            # Only first element of list will be used as input
            if type(options.input) is list:
                self.fInput = options.input[0]
            else:
                self.fInput = options.input
        
        if hasattr(options, 'output'):
            # Write to file
            self.writeFile = True

            # Only first element of list will be used as output
            if type(options.output) is list:
                self.fOutput = options.output[0]
            else:
                self.fOutput = options.output
        # Write to stdout
        else:
            self.fOutput = ""
        
        if hasattr(options, 'minAlignmentQuality'):
            self.minAlignmentQuality = options.minAlignmentQuality
        
        if hasattr(options, 'minReadLength'):
            self.minReadLength = options.minReadLength
            
        if hasattr(options, 'maxReadLength'):
            self.maxReadLength = options.maxReadLength
        
        if hasattr(options, 'maxReadIntervalLength'):
            self.maxReadIntervalLength = options.maxReadIntervalLength
            
        if hasattr(options, 'primary'):
            self.primary = options.primary
            
        if hasattr(options, 'choice'):
            self.choice = options.choice
		
        if hasattr(options, 'mate'):
			self.mate = options.mate
           
        self.data = {'maxReadLength' : self.maxReadLength,
                     'minReadLength' : self.minReadLength,
                     'primary': self.primary,
                     'maxReadIntervalLength': self.maxReadIntervalLength,
                     'minAlignmentQuality': self.minAlignmentQuality}

    #================================================================================ 
    '''
    Printing function
    '''
    def write(self, s):
        if self.writeFile:
            self.fOutput.write(s)
        else:
            print(s)
              
    #================================================================================ 
    '''
    This method determines if a read fullfills the critera to be included in the analysis.
    '''
    def readFullfillsQualityCriteria(self, almnt):
        if self.mate == 2: #select the second read of a pair to extract
            if almnt.paired_end and almnt.pe_which == "second":
                return(almnt.aligned and almnt.iv.length <= self.data['maxReadIntervalLength'] and
                       almnt.aQual >= self.data['minAlignmentQuality'] and not almnt.failed_platform_qc
                       and self.primaryFilter(almnt))
            elif almnt.paired_end and almnt.pe_which == "first":
                self.count += 1
                return False
            elif not almnt.paired_end:
                raise("Data are not paired end %s" % almnt.paired_end)
                
                
        elif self.mate == 1: #select the first read of a pair to extract
            if ( not almnt.paired_end ) or (almnt.paired_end and almnt.pe_which =="first"):
                return(almnt.aligned and
                       almnt.iv.length <= self.data['maxReadIntervalLength'] and
                       almnt.aQual >= self.data['minAlignmentQuality'] and
                       not almnt.failed_platform_qc and #SAM flag 0x0200
                       self.primaryFilter(almnt))
            elif almnt.paired_end and almnt.pe_which == "second":
                self.count += 1
                return False
        else:
            raise ValueError("mate argument can only be 1 or 2")

    '''
    If the primary filter is activated (option primary is set),
    for multimapping reads it will filter out the best location
    using the SAM flag 0x0100
    '''
    def primaryFilter(self, almnt):
        if self.data['primary']:
            return(almnt.not_primary_alignment)
        else:
            return(True) 

   '''
    Calculates the read length and stores the min and max length
    '''
    def getSequenceLength(self, almnt):
        return(len(almnt.read.seq))   
       
    def calcMinMax(self, almnt):
        length = self.getSequenceLength(almnt)
      
        if not (self.data['maxReadLength'] == 0 and self.data['minReadLength'] == 0):
            if length > self.data['maxReadLength']:
                self.data['maxReadLength'] = length
            if length < self.data['minReadLength']:
                self.data['minReadLength'] = length
        return(length)

   '''
    Returns a list of Genomic intervals for the positions of deletions/insertions/mutations
    '''
    def parseCigar(self, almnt):
        variations = {'deletions': list(),
                     'insertions': list(),
                     'hits': list() }
      
        for i in almnt.cigar:
            if i.type == 'D':
                variations['deletions'].append(i)
            elif i.type == 'I':
                variations['insertions'].append(i)
            elif i.type == 'M':
                variations['hits'].append(i)
         
        return(variations)

   def posCalcStartSite(self, pos_d, strand, x):
        if strand == "+":
            return(pos_d + x)
        elif strand == "-":
            return(pos_d + x)
        else:
            raise("Strand not known %s" % strand)

    def posCalcMiddleSite(self, pos_d, strand, x):
        if strand == "+":
            return(pos_d + x)
        elif strand == "-":
            return(pos_d - x)
        else:
            raise("Strand not known %s" % strand)

    def posCalcEndSite(self, pos_d, strand, x):
        if strand == "+":
            return(pos_d + x)
        elif strand == "-":
            return(pos_d + x)
        else:
            raise("Strand not known %s" % strand)

   '''
    Returns GenomicPosition for desired site based on start site
    '''
    def getStartSiteAsBed_firstread(self, almnt):
        dp = self.choice.split("s")
        
        ignore = False
        if dp[1] == '':
            dp = 0
        else:
            if "i" in dp[1]:
                ignore = True     
                dp = dp[1].split("i")
                dp = int(dp[0])
            else:
                dp = int(dp[1])

        if almnt.iv.strand == "+":
            x = almnt.iv.start_d + dp
        elif almnt.iv.strand == "-":
            x = almnt.iv.start_d - dp
        else:
            raise("Strand not known %s" % almnt.iv.strand)

        if x < 0:
            if ignore == False:
                error = "Value Error: Start position cannot be less than zero. Alignment: " + str(almnt.iv) + ", Read: " + almnt.read.name +  ". Check your data!"
                raise ValueError(error)
            
            return None
        else:
        
            y = x+1
            
            try:
                yb = almnt.optional_field('YB')
            except Exception:
                yb = 1
                pass
            
            seq = (almnt.iv.chrom, str(min(x,y)), str(max(x,y)), almnt.read.name+"|"+str(len(almnt.read.seq)), str(yb), almnt.iv.strand)
    
            return(str("\t").join(seq))
	

    def getStartSiteAsBed_secondread(self, almnt):
        dp = self.choice.split("s")
        
        ignore = False
        if dp[1] == '':
            dp = 0
        else:
            if "i" in dp[1]:
                ignore = True
                dp = dp[1].split("i")
                dp = int(dp[0])
            else:
                dp = int(dp[1])

        if almnt.iv.strand == "+":
            x = almnt.iv.end_d - 1
        elif almnt.iv.strand == "-":
            x = almnt.iv.end_d + 1
        else:
            raise("Strand not known %s" % almnt.iv.strand)

        if x < 0:
            if ignore == False:
                error = "Value Error: Start position cannot be less than zero. Alignment: " + str(almnt.iv) + ", Read: " + almnt.read.name +  ". Check your data!"
                raise ValueError(error)
            return None
        else:
            y = x+1
            try:
                yb = almnt.optional_field('YB')
            except Exception:
                yb = 1
                pass

            seq = (almnt.iv.chrom, str(min(x,y)), str(max(x,y)), almnt.read.name+"|"+str(len(almnt.read.seq)), str(yb), almnt.iv.strand)

            return(str("\t").join(seq))
 
    '''
    Returns GenomicPosition for middle site and supports Gapped alignments!
    '''
    def determineMiddleSite(self, almnt):

        pos = round(decimal.Decimal(len(almnt.read.seq)) / 2, 0)

        cigarList = almnt.cigar

        if almnt.iv.strand == "-":
            cigarList = reversed(cigarList)

        for cig in cigarList:

            if cig.type == "M" or cig.type == "I":
                if cig.size - pos > 0:
                    pos = self.posCalcMiddleSite(cig.ref_iv.start_d, cig.ref_iv.strand, pos)
                    break
                else:
                    pos = pos - cig.size
        return(HTSeq.GenomicPosition(almnt.iv.chrom, pos, almnt.iv.strand))

    def getMiddleSiteAsBed(self, almnt):
        pos = self.determineMiddleSite(almnt)
        x = pos.start_d
        y = x + 1

        try:
            yb = almnt.optional_field('YB')
        except Exception:
            yb = 1
            pass
        
        seq = (almnt.iv.chrom, str(min(x,y)), str(max(x,y)), almnt.read.name+"|"+str(len(almnt.read.seq)), str(yb), almnt.iv.strand)
        
        return(str("\t").join(seq))
 
    '''
    Returns GenomicPosition for end site
    '''
    def determineEndSite(self, iv):
        return(HTSeq.GenomicPosition(iv.chrom, iv.end_d, iv.strand))

    def getEndSiteAsBed(self, almnt):
        x = almnt.iv.end_d
        y = x + 1
        
        try:
            yb = almnt.optional_field('YB')
        except Exception:
            yb = 1
            pass
        
        seq = (almnt.iv.chrom, str(min(x,y)), str(max(x,y)), almnt.read.name + "|" + str(len(almnt.read.seq)), str(yb), almnt.iv.strand)

        return(str("\t").join(seq))   
    

    '''
    Writes deletions of an alignment parsed by CIGAR string  as bed line
    '''
    def writeDeletionSiteAsBed(self, almnt, fOutput):
           
        seq = ()
        deletion = []
        
        #Goes through cigar string and only searches for deletions
        #if no deletion is found the read has none
        for i in almnt.cigar:
            if i.type == 'D':
                deletion.append(i)
        
        #if the read has a deletion it will be written out in the bed file
        #for deletions query_from and query_to are always the same
        if deletion:
            
            for i in range(len(deletion)):
                x = deletion[i].ref_iv.start_d
                y = deletion[i].ref_iv.end_d
                if almnt.iv.strand == '-':
                    b = x
                    x = y
                    y = b
                seq = (almnt.iv.chrom, str(x), str(y), almnt.read.name + "|" + str(len(almnt.read.seq)), str(deletion[i].query_from), almnt.iv.strand)
                seq = (str("\t").join(seq)) 
                fOutput.write(seq + "\n")    
   
    '''
    Writes insertions of an alignment parsed by CIGAR string as bed line
    '''
    def writeInsertionSiteAsBed(self, almnt, fOutput):
           
        seq = ()
        insertion = []
        
        #Goes through cigar string and only searches for insertions
        #if no insertion is found the read has none
        for i in almnt.cigar:
            if i.type == 'I':
                insertion.append(i)
        
        #if the read has an insertion it will be written out in the bed file
        if insertion:
            
            for i in range(len(insertion)):
                x = insertion[i].ref_iv.start_d
                y = insertion[i].ref_iv.end_d
                if almnt.iv.strand == '-':
                    b = x
                    x = y
                    y = b
                seq = (almnt.iv.chrom, str(x), str(y), almnt.read.name+"|"+str(len(almnt.read.seq)), str(insertion[i].query_from), almnt.iv.strand)
                seq = (str("\t").join(seq)) 
                fOutput.write(seq + "\n")     


    '''
    Returns start site as bed line 
    '''
    def getStartSiteAsBed(self, almnt):
        if self.mate == 1:
            return(self.getStartSiteAsBed_firstread(almnt))
        elif self.mate == 2:
            return(self.getStartSiteAsBed_secondread(almnt))
        else:
            raise ValueError("mate argument must be 1 or 2")

    '''
    Extract start sites
    '''
    def extract_StartSites(self):
        almnt_file = HTSeq.BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)

        for almnt in almnt_file:
                if self.readFullfillsQualityCriteria(almnt):
                    out = self.getStartSiteAsBed(almnt)
                    if not out == None:
                        fOutput.write(out + "\n")

        fOutput.close()

    '''
    Extract middle sites
    '''   
    def extract_MiddleSites(self):
        almnt_file = HTSeq.BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
    
        for almnt in almnt_file:
            if self.readFullfillsQualityCriteria(almnt):
                fOutput.write(self.getMiddleSiteAsBed(almnt) + "\n")
       
        fOutput.close() 
           
    '''
    Extract end sites. 
    '''   
    def extract_EndSites(self):
        almnt_file = HTSeq.BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
    
        for almnt in almnt_file:
            if self.readFullfillsQualityCriteria(almnt):
                fOutput.write(self.getEndSiteAsBed(almnt) + "\n")

        fOutput.close()
        
    '''
    Extract deletion sites.
    Deletion sites are determined by parsing the CIGAR string.
    ''' 
    def extract_DeletionSites(self):
        almnt_file = HTSeq.BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
    
        for almnt in almnt_file:
            if self.readFullfillsQualityCriteria(almnt):
                self.writeDeletionSiteAsBed(almnt, fOutput)
       
        fOutput.close()
        
    '''
    Extract insertion sites.
    Insertion sites are determined by parsing the CIGAR string.
    ''' 
    def extract_InsertionSites(self):
        almnt_file = HTSeq.BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
    
        for almnt in almnt_file:
            if self.readFullfillsQualityCriteria(almnt):
                self.writeInsertionSiteAsBed(almnt, fOutput)
       
        fOutput.close()      
