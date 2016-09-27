__author__ = 'Nadia'

# --------------------------------------------------
# featureCLIP class
# Authors: Nadia Ashraf, nadia.ashraf@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: May 2016
# --------------------------------------------------

import gzip, HTSeq
from collections import OrderedDict

class feature:

    fInput = ""   #crosslink sites file .bed format
    fOutput = ""  #output file .txt format
    fCompare = ""  #input file (repeats file) .bed format
    choice = ""
    dist = 200     #upstream/downstream search window size
    strand = 'y'   #if strand info is provided or not
    score = 'y'    #if score is provided or not

    def __init__(self, options):

        if hasattr(options, 'input'):
            self.fInput = options.input[0]

        if hasattr(options, 'output'):
            self.fOutput = options.output[0]

        if hasattr(options, 'compare'):
            self.fCompare = options.compare

        if hasattr(options, 'choice'):
            self.choice = options.choice

        if hasattr(options, 'dist'):
            self.dist = options.dist

        if hasattr(options, 'strand'):
            self.strand = options.strand

        if hasattr(options, 'score'):
            self.score = options.score



    #=====================================================
    """
    Store data in dictionary for comparison. Assumes that each feature name is unique.
    Format: { chromosome : { strand : [(Start postion, end postion, name), (...), ...] }}
    """

    def build_dict(self,data_file):

        d = {}

        for almnt in data_file:

            if self.score == 'y' and self.strand == 'n':

                if not d.has_key(almnt.iv.chrom):
                        d[almnt.iv.chrom] = {'.' : [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)]]}
                else:
                    if not d[almnt.iv.chrom].has_key('.'):
                        d[almnt.iv.chrom]['.'] = [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)]]
                    else:
                        d[almnt.iv.chrom]['.'].append([almnt.iv.start_d, almnt.iv.end_d, almnt.name, int(almnt.score)])

            elif self.score == 'n' and self.strand == 'y':

                if almnt.iv.strand == '+':
                    if not d.has_key(almnt.iv.chrom):
                        d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, 0]]}

                    else:
                        if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                            d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, 0]]
                        else:
                            d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.start_d, almnt.iv.end_d, almnt.name, 0])

                else:
                    if not d.has_key(almnt.iv.chrom):
                        d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.end_d, almnt.iv.start_d, almnt.name, 0]]}

                    else:
                        if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                            d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.end_d, almnt.iv.start_d, almnt.name, 0]]
                        else:
                            d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.end_d, almnt.iv.start_d, almnt.name, 0])


            elif self.score == 'y' and self.strand == 'y':
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

            else:
                if not d.has_key(almnt.iv.chrom):
                        d[almnt.iv.chrom] = {'.' : [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, 0]]}
                else:
                    if not d[almnt.iv.chrom].has_key('.'):
                        d[almnt.iv.chrom]['.'] = [[almnt.iv.start_d, almnt.iv.end_d, almnt.name, 0]]
                    else:
                        d[almnt.iv.chrom]['.'].append([almnt.iv.start_d, almnt.iv.end_d, almnt.name, 0])

        return d

    def count_all(self):

        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w')
        else:
            output = open(self.fOutput, 'w')
            seq = ('Chromosome','Region start pos','Region end pos','Region name','Strand','Length(nt)', 'Toatal cross-link sites in region',
                   'Number of pos where cross-links are located', 'Max count in one pos','Count density','Total remove duplicates','Max no.of removed duplicates')
            output.write(str("\t").join(seq) + "\n")

        #Get the information for normalisation of the plots

        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)

        d1 = self.build_dict(almnt_file1)
        d2 = self.build_dict(almnt_file2)

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

                        name = b[2]


                        seq = (chrom, str(b[0]+1), str(b[1]+1), name, strand, str(length), str(0), str(0), str(0), str(0), str(0), str(0))
                        output.write(str("\t").join(seq) + "\n")

                    continue

                #if the input file contains reads on the current chromosome but not on the same
                #strand, then write out all positions with zero
                elif not d1[chrom].has_key(strand):
                    for b in B:

                        length = b[1] - b[0]

                        name = b[2]

                        seq = (chrom, str(b[0]+1), str(b[1]+1), name, strand,str(length), str(0), str(0), str(0), str(0), str(0), str(0))
                        output.write(str("\t").join(seq) + "\n")
                    continue

                A = d1[chrom][strand]

                self.calculateCount(A, B, chrom, strand, output)

        output.close()

    def count_only(self):

        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w')
        else:
            output = open(self.fOutput, 'w')
            seq = ('Chromosome','Region start pos','Region end pos','Region name','Region score','Strand','Length(nt)',
                  'Total cross-link sites in region','Number of pos where cross-links are located', 'Max count in one pos','Count density',
                  'Total remove duplicates','Max no.of removed duplicates')

            output.write(str("\t").join(seq) + "\n")

        #Get the information for normalisation of the plots


        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)

        #self.process(almnt_file2,output)

        d1 = self.build_dict(almnt_file1)
        d2 = self.build_dict(almnt_file2)

        if (not d1.has_key('chr1')) or (not d2.has_key('chr1')):
            error = "Chromosome number format is not same, please make sure it's same for both files!"
            raise ValueError(error)

        else:
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

                        check = False
                        continue
                    elif a[0] > b_Last[1]:

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
                                check = False


    def writeOut(self, chrom, strand, b, d_count, d_dup, length, output):

        name = b[2]


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

            seq = (chrom, str(b[0]+1), str(b[1]+1), name, str(b[3]), strand, str(length), str(sum(d_count.values())), str(counts), str(d_count[m_count]), str(density), str(dup_counts), str(d_dup[m_dup]))
            #seq = (chrom, str(b[0]+1), str(b[1]+1), name, strand, str(length), str(sum(d_count.values())), str(counts), str(d_count[m_count]), str(density), str(dup_counts), str(d_dup[m_dup]))
            output.write(str("\t").join(seq) + "\n")
        else:
            seq = (chrom, str(b[0]+1), str(b[1]+1), name, str(b[3]), strand, str(length), str(0), str(0), str(0), str(0), str(0), str(0))
            #seq = (chrom, str(b[0]+1), str(b[1]+1), name, strand, str(length), str(0), str(0), str(0), str(0), str(0), str(0))
            output.write(str("\t").join(seq) + "\n")

    """
    Function for processing repeats file. Get track information
    """

    def process(self, rmsk_file,output):

        typeFooter = {}
        chromFooter = {}

        for repeat in rmsk_file:
            repeat.iv.strand = '.'
            if not chromFooter.has_key(repeat.iv.chrom):
                chromFooter[repeat.iv.chrom] = (repeat.iv.end - repeat.iv.start)
            else:
                chromFooter[repeat.iv.chrom] += (repeat.iv.end - repeat.iv.start)

            if not typeFooter.has_key(repeat.name):
                typeFooter[repeat.name] = (repeat.iv.end - repeat.iv.start)
            else:
                typeFooter[repeat.name] += (repeat.iv.end - repeat.iv.start)

        chromFooter = OrderedDict(sorted(chromFooter.items(), key=lambda x: (-x[1], x[0])))
        for k in chromFooter:
            output.write("#track chr "+str(k)+" "+str(chromFooter[k])+"\n")

        typeFooter = OrderedDict(sorted(typeFooter.items(), key=lambda x: (-x[1], x[0])))

        for k in typeFooter:
            output.write("#track type "+str(k)+" "+str(typeFooter[k])+"\n")


    """
    Function for identifyig nearest cross-link site to a feature
    """

    def dist_cl(self):
        if self.fOutput.endswith(".gz"):
            output = gzip.open(self.fOutput, 'w')
        else:
            output = open(self.fOutput, 'w')
            seq = ('Chromosome','Region start pos','Region end pos','Region name','Region score','Strand','Length in nt',
                   'Total cross-link sites in region','Number of pos where cross-links are located', 'Max count in one pos','Position of nearest cross link site','Distance to nearest cross-link site','Count density',
                   'Total remove duplicates','Max no.of removed duplicates')
            output.write(str("\t").join(seq) + "\n")

        #Get the information for normalisation of the plots


        almnt_file1 = HTSeq.BED_Reader(self.fInput)
        almnt_file2 = HTSeq.BED_Reader(self.fCompare)

        d1 = self.build_dict(almnt_file1)
        d2= self.built_dict_SW(almnt_file2,self.dist)

        for chrom in d1:
            if not d2.has_key(chrom):
                continue
            for strand in d1[chrom]:
                if not d2[chrom].has_key(strand):
                    continue

                A = d1[chrom][strand]
                B = d2[chrom][strand]

                self.get_cl_dist(A, B, chrom, strand, output)

        output.close()


    """
    Function for getting dictionary for upstream and downstream regions of features
    """
    def built_dict_SW(self,data_file,window):
        d = {}

        for almnt in data_file:

            if almnt.iv.strand == '+':
                if not d.has_key(almnt.iv.chrom):
                    d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.start_d-window, almnt.iv.start_d, 'upstream_'+almnt.name, int(almnt.score)]]}
                    d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.end_d, almnt.iv.end_d+window, 'downstream_'+almnt.name, int(almnt.score)]]}
                else:
                    if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                        d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.start_d-window, almnt.iv.start_d, 'upstream_'+almnt.name, int(almnt.score)]]
                        d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.end_d, almnt.iv.end_d+window, 'downstream_'+almnt.name, int(almnt.score)]]
                    else:
                        d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.start_d-window, almnt.iv.start_d, 'upstream_'+almnt.name, int(almnt.score)])
                        d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.end_d, almnt.iv.end_d+window, 'downstream_'+almnt.name, int(almnt.score)])

            else:
                if not d.has_key(almnt.iv.chrom):
                    d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.end_d-window, almnt.iv.end_d, 'downstream_'+almnt.name, int(almnt.score)]]}
                    d[almnt.iv.chrom] = {almnt.iv.strand : [[almnt.iv.start_d, almnt.iv.start_d+window, 'upstream_'+almnt.name, int(almnt.score)]]}
                else:
                    if not d[almnt.iv.chrom].has_key(almnt.iv.strand):
                        d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.end_d-window, almnt.iv.end_d, 'downstream_'+almnt.name, int(almnt.score)]]
                        d[almnt.iv.chrom][almnt.iv.strand] = [[almnt.iv.start_d, almnt.iv.start_d+window, 'upstream_'+almnt.name, int(almnt.score)]]
                    else:
                        d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.end_d-window, almnt.iv.end_d, 'downstream_'+almnt.name, int(almnt.score)])
                        d[almnt.iv.chrom][almnt.iv.strand].append([almnt.iv.start_d, almnt.iv.start_d+window, 'upstream_'+almnt.name, int(almnt.score)])

        return d

    """
    Function for counting CLs outside feature and determining nearest CL site to the feature
    """
    def get_cl_dist(self,A,B,chrom,strand,output):

        #A is dictionary of all cross link (cl) sites
        #B is dictionary of annotation/features/regions. Contains genomic start end positions of all regions, strand information (if available), chromosome number, alignment score(if available),
        #and name of the region .

        #get first region and count how many cl sites fall in that region
        #loop over cls and compare their position with the region coordinates.
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



            #length of feature
            length = 0
            #data structure for analysis
            d_count = {}
            d_dup = {}
            cl_dist = {}

            for a in A:

                check = True

                while check:
                    #if smaller then intergenic region before first region
                    #if bigger intergenic region after last region
                    #else its in the regions

                    if a[1] < b_First[0]:

                        check = False
                        continue
                    elif a[0] > b_Last[1]:

                        check = False
                    else:
                        b_Curr = B[bi]

                        #if bigger search go for the next region and write out
                        #the last region if zero nothing is found if not zero there are positions found
                        #else count in the current region where the cross-link site is in the positions
                        if a[0] > b_Curr[1] or finished == True:

                            length = b_Curr[1] - b_Curr[0]

                            if len(d_count) != 0:
                                self.write_out(chrom, strand, b_Curr, d_count, d_dup, length,cl_dist, output)

                            if not b_Curr == b_Last:
                                bi = bi + 1
                            else:
                                check = False
                            d_count = {}
                            d_dup = {}
                            cl_dist = {}

                        else:

                            if a[0] >= b_Curr[0] and a[1] <= b_Curr[1]:

                                #if nothing is in the feature generate the data structure
                                if len(d_count) == 0:
                                    length = b_Curr[1] - b_Curr[0]

                                    #create dictionary instances for storing cout,duplication, and distance information
                                    for i in range(length+1):
                                        d_count[i+b_Curr[0]] = 0
                                        d_dup[i+b_Curr[0]] = 0
                                        cl_dist[i+b_Curr[0]] = 0

                                #for count add 1 every time same cross link position is hit
                                #for duplicate add the duplicate information from cross links bed file
                                d_count[a[0]] += 1
                                d_dup[a[0]] += a[3]

                                #for nearest distance calculte start - cl position and end - cl position and keep the minimum of the two.
                                #separate calculations for forward and reverse strands.

                                #upstream region; get distance of cl from start pos of region and end pos of the previous region and keep the minimum of the two
                                #downstream region; get distance of cl form end pos of region and start pos of next region and keep minimum of the two.


                                if 'upstream' in b_Curr[2]:
                                    cl_dist[a[0]] += min(abs(b_Curr[0]-a[0]), abs(a[0]-B[bi-1][1]))
                                else:
                                    cl_dist[a[0]] += min(abs(a[0]-b_Curr[1]), abs(B[bi+1][0]-a[0]))


                                #if last cross-link site found write out
                                if a == A[-1]:
                                    finished = True

                                    length = b_Curr[1] - b_Curr[0]

                                    self.write_out(chrom, strand, b_Curr, d_count, d_dup, length,cl_dist, output)

                                    if not b_Curr == b_Last:
                                        bi = bi + 1
                                    else:
                                        check = False

                                    d_count = {}
                                    d_dup = {}
                                    cl_dist = {}
                                else:
                                    check = False
                            else:
                                check = False

    def write_out(self, chrom, strand, b, d_count, d_dup, length, dist,output):

        name = b[2]


        if len(d_count) > 0:
            min_dist = 0
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

            min_pos = min(dist.keys(), key=(lambda k: d_dup[k]))

            density = float(counts) / float(length)
            min_dist= min(abs((b[0]+1) - min_pos),abs((b[1]+1)-min_pos))

            seq = (chrom, str(b[0]+1), str(b[1]+1), name, str(b[3]), strand, str(length), str(sum(d_count.values())), str(counts), str(d_count[m_count]),str(min_pos),str(min_dist), str(density), str(dup_counts), str(d_dup[m_dup]))
            output.write(str("\t").join(seq) + "\n")
        else:
            seq = (chrom, str(b[0]+1), str(b[1]+1), name, str(b[3]), strand, str(length), str(0), str(0), str(0),str(0), str(0), str(0), str(0), str(0))
            output.write(str("\t").join(seq) + "\n")


