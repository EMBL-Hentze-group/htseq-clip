# --------------------------------------------------
# fastaCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: February 2016
# --------------------------------------------------

import os, gzip

try:
    from Bio import SeqIO
except Exception:
    print "Please install the biopython framework e.g. like this"
    print "pip install biopython"
    print "pip install biopython --user"
    os._exit(1)

try:
    import HTSeq
except Exception:
    print "Please install the HTSeq framework e.g. like this"
    print "pip install HTSeq"
    print "pip install HTSeq --user"
    os._exit(1)
       
class fastaCLIP:
    
    input = ""
    output = ""
    readLength = 0

    def __init__(self, options):
        
        if hasattr(options, 'input'):
            self.input = options.input
        
        if hasattr(options, 'output'):
            self.output = options.output
            
        if hasattr(options, 'maxReadLength'):
            self.readLength = options.maxReadLength
    
    '''
    Function that formats a .fa or .fasta whole genome file into a .fastq
    by a given sequence length
    '''        
    def genomeToReads(self):
        
        qcString = ""
        
        #Quality string for fastq format
        for i in range(0,self.readLength):
            qcString += "~"
            
        genome = {}
        
        #get Sequences by Chromosomes    
        handle = open(self.input, "rU")
        for record in SeqIO.parse(handle, "fasta") :
            if not genome.has_key(record.id):
                genome[record.id] = record.seq
        handle.close()
        
        print "Data preprocessing done!"
        
        #Write out in fastq format   
        if self.output.endswith(".gz"):
            fqOutput = gzip.open(self.output, 'w') 
        else:        
            fqOutput = open(self.output, 'w') 
        
        #foreach chromosome in the genome
        for chrom in genome:
            
            count = 1
            seq = ""
            
            #go through the sequence and write out sequences
            #with the given length
            for i in range(0, len(genome[chrom])):
                
                #if len(seq) smaller, then a sequence by the given length is generated
                #if it the length is equal to the given length write it out
                #and delete the first base to go one position further
                #then u always add one base write out and delete the first one   
                if len(seq) < self.readLength:
                    seq += genome[chrom][i]
                        
                    if len(seq) == self.readLength:
      
                        fqOutput.write("@"+chrom+":"+str(count)+":AAAAAAAAA"+"\n")
                        fqOutput.write(seq+"\n")
                        fqOutput.write("+"+"\n")
                        fqOutput.write(qcString+"\n")
                                
                        seq = seq[1:]
                        rlCount = 0
                        count += 1
                        
                    if count % 2000000 == 0:
                        print chrom + " proceeding ..."
                        
            print chrom + " finished!"
            
        print "Genome to reads completed!"
           
        fqOutput.close()
            
