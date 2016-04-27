# --------------------------------------------------
# fastaCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: February 2016
# --------------------------------------------------

import gzip, HTSeq, itertools
from Bio import SeqIO
       
class fastaCLIP:
    
    fInput = ""
    fOutput = ""
    mm = 0
    bcp = ""
    readLength = 0

    def __init__(self, options):
        
        if hasattr(options, 'input'):
            self.fInput = options.input
        
        if hasattr(options, 'output'):
            self.fOutput = options.output
            
        if hasattr(options, 'maxReadLength'):
            self.readLength = options.maxReadLength
            
        if hasattr(options, 'barcode'):
            self.bcp = options.barcode
            
        if hasattr(options, 'mismatches'):
            self.mm = options.mismatches
    
    #===================================================================================  
    #===================================================================================            
    '''
    Function that removes the multiplex barcode from 
    the full barcode for the removal of duplicates
    '''   
    def get_BC_Dictionary(self, fastq, bcp):   
        d = {}
        #Convert the barcode position string into 4 integers
        bcp = [int(i) for i in bcp]
        
        #Get all barcodes without the multiplex
        for read in itertools.islice(fastq, None):
            full_barcode = read.name.split(":")[-1]
            bc = full_barcode[bcp[0]:bcp[1]] + full_barcode[bcp[2]:bcp[3]] 
            l = len(read)
            key = (l, bc)
            if not d.has_key(key):
                d[key] = []
            d[key].append(read)
        
        return d
    
    #=================================================================================== 
    '''
    Function that removes the duplicates from a read
    if necessary
    '''
    def remove_Duplicates(self, d, key, mm, out_fastq_file):    
        l = d[key]
        reads = []
        pos = []
                
        #Append the read sequences to a list for better comparison
        for read in l:
            reads.append(read.seq)
        
        #Checking if there are duplicates which should be removed 
        #and writing the others into the output file
        for i in range(len(reads)):
            s1 = reads[i]
            count = 0
            for j in range(i+1, len(reads)):
                s2 = reads[j]     
                mismatches = [x for x in xrange(len(s1)) if s1[x] != s2[x]]
                if len(mismatches) <= mm:
                    count += 1
                    if j not in pos:
                        pos.append(j)
                       
            if count > 0 and i not in pos:
                out_fastq_file.write("@"+l[i].name+":"+str(count)+"\n")
                out_fastq_file.write(l[i].seq+"\n")
                out_fastq_file.write("+\n")
                out_fastq_file.write(l[i].qualstr+"\n")                 
            elif len(pos) == 0:
                out_fastq_file.write("@"+l[i].name+":0\n")
                out_fastq_file.write(l[i].seq+"\n")
                out_fastq_file.write("+\n")
                out_fastq_file.write(l[i].qualstr+"\n") 
     
    
    #=================================================================================== 
    '''
    Function that removes the random barcode duplicate
    '''               
    def remove(self):  
        
        #Getting arguments from user
        in_fastq_file = HTSeq.FastqReader(self.fInput)
        outFileName = self.fOutput
        mm = self.mm
        bcp = self.bcp
        
        #Generating output file
        out_fastq_file = open(outFileName, "w")
        
        #Generate a barcode dictionary for the removal of duplicates
        d = self.get_BC_Dictionary(in_fastq_file, bcp)
        
        # For each barcode, chose one read to represent
        # and store the number of reads there was originally
        # in the representant's fastq file in a new fastq file
        for key in d:
            if len(d[key]) < 2:
                read = d[key][0]
                out_fastq_file.write("@"+read.name+":0\n")
                out_fastq_file.write(read.seq+"\n")
                out_fastq_file.write("+\n")
                out_fastq_file.write(read.qualstr+"\n")
            else:
                self.remove_Duplicates(d, key, mm, out_fastq_file)      
            
        out_fastq_file.close()       
    #===================================================================================
    #===================================================================================          
    
    #===================================================================================         
    #===================================================================================
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
        handle = open(self.fInput, "rU")
        for record in SeqIO.parse(handle, "fasta") :
            if not genome.has_key(record.id):
                genome[record.id] = record.seq
        handle.close()
        
        print "Data preprocessing done!"
        
        #Write out in fastq format   
        if self.fOutput.endswith(".gz"):
            fqOutput = gzip.open(self.fOutput, 'w') 
        else:        
            fqOutput = open(self.fOutput, 'w') 
        
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
                        count += 1
                        
                    if count % 2000000 == 0:
                        print chrom + " proceeding ..."
                        
            print chrom + " finished!"
            
        print "Genome to reads completed!"
           
        fqOutput.close()
    #=================================================================================== 
    #=================================================================================== 
            
