# --------------------------------------------------
# gtfCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
# Institution: EMBL Heidelberg
# Date: May 2016
# --------------------------------------------------

import gzip, HTSeq
from collections import OrderedDict


class gtfClip:

    features = {}
    fInput = ''
    gtfFile = ''
    fOutput = ''
    geneType = ''
    geneName= ''
    windowSize = 50
    windowStep = 20

    def __init__(self, options):

        if hasattr(options, 'gtf'):
            self.gtfFile = options.gtf

        if hasattr(options, 'input'):
            self.fInput = options.input

        if hasattr(options, 'output'):
            self.fOutput = options.output

        if hasattr(options, 'windowSize'):
            self.windowSize = options.windowSize

        if hasattr(options, 'windowStep'):
            self.windowStep = options.windowStep

        if hasattr(options, 'type'):
            self.geneType = options.type

        if hasattr (options,'name'):
            self.geneName = options.name

    '''
    This method processes through the gtf file and determines the positions of exons and introns
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
        tss = 0
        sc_start = 0
        sc_stop = 0
        ec = 0
        ic = 0
        ps = [] #list of tuples with each tuple storing start/end positions of exons
        ps_cds = []
        eFlag = None
        fs = []
        ts = []
        sc = []
        ss = []
        eCount = 0
        iCount = 0

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
                    gene = list(gas[HTSeq.GenomicInterval(chrom, start, end, strand)].steps())
                    #Sort the list where all exons are in
                    ps.sort(key=lambda tup: tup[0])

                    # if feature is protein codin process forward and reverse strand individually
                    if t == 'protein_coding':
                        out, ecount = self.process(gene, ps_cds, typeFooter, tss, sc_start, sc_stop, name, ec, ic)
                    else:
                        eCount = (len(gene)/2)+1
                        iCount = len(gene)/2

                        for g in gene:
                            if not g[1]:
                                if not typeFooter.has_key(t):
                                    typeFooter[t] = (g[0].end - g[0].start)
                                else:
                                    typeFooter[t] += (g[0].end - g[0].start)

                                seq = [g[0].chrom, str(g[0].start), str(g[0].end), name+'@intron@'+str(ic)+'/'+str(iCount), 1, g[0].strand, "i"]
                                fs.append(seq)
                                ic += 1

                                #Depending on if current feature is exon or intron calculate the correct position
                                #and flag the exons
                            else:
                                while len(ps) > 0 and ps[0][0] < g[0].end:
                                    if ps[0][0] == g[0].start and ps[0][1] == g[0].end:
                                        eFlag = 3
                                    elif ps[0][0] == g[0].start and ps[0][1] < g[0].end and eFlag != 0:
                                        eFlag = 2
                                    elif ps[0][0] > g[0].start and ps[0][1] == g[0].end and eFlag != 0:
                                        eFlag = 1
                                    elif  ps[0][0] > g[0].start and ps[0][1] < g[0].end and eFlag != 0:
                                        eFlag = 0
                                    del ps[0]

                                if not typeFooter.has_key(t):
                                    typeFooter[t] = (g[0].end - g[0].start)
                                    sec_tag = 'exon'
                                else:
                                    typeFooter[t] += (g[0].end - g[0].start)
                                    sec_tag = 'exon'
                                seq = [g[0].chrom, str(g[0].start), str(g[0].end), name+'@'+sec_tag+'@'+str(ec)+'/'+str(eCount), eFlag, g[0].strand, "e"]
                                fs.append(seq)
                                ec += 1

                    #Flagging of introns and writing out in output file
                    #else there is only one exon
                    if t != 'protein_coding':
                        if not len(fs) == 1:
                            # Initialisation of Flags

                            if len(fs) == 2:
                                for line in fs:
                                    if line[6] == "e":
                                        seq = (line[0], line[1], line[2], line[3], str(line[4]), line[5], "\n")
                                        output.write(str('\t').join(seq))
                                    elif line[6] == "i":
                                        seq = (line[0], line[1], line[2], line[3], str(line[4]), line[5], "\n")
                                        output.write(str('\t').join(seq))
                            else:
                                for val in range(len(fs)):
                                    if fs[val][6] == 'i':

                                        p = val-1
                                        eFlag1 = fs[p][4]
                                        p = p + 1
                                        eFlag2 = fs[p][4]
                                    #Simple binary operations to calculate the flag for the introns
                                    #The intron flag is calculated by using the left exon flag and the
                                    #right exon flag by using the bits of the binary coded flag
                                    # for further information how they are calculated look up the
                                    # documentation
                                        fs[val][4] = ((eFlag1 & 1) << 1) | (eFlag2 >> 1)

                                        seq = (fs[val][0], fs[val][1], fs[val][2], fs[val][3], str(fs[val][4]), fs[val][5], "\n")
                                        output.write(str('\t').join(seq))


                                    elif fs[val][6] == "e":
                                        seq = (fs[val][0], fs[val][1], fs[val][2], fs[val][3], str(fs[val][4]), fs[val][5], "\n")
                                        output.write(str('\t').join(seq))

                        else:
                                    #if gene type is is something else, you have to flag it with 3 because there only exists one isoform
                            if fs[0][4] == None:
                                fs[0][4] = 3
                            seq = (fs[0][0], fs[0][1], fs[0][2], fs[0][3], str(fs[0][4]), fs[0][5], "\n")
                            output.write(str('\t').join(seq))
                    else:
                        if not len(out) == 1:

                            for val in range(len(out)):
                                if out[val][6] == 'i':

                                    p = val-1
                                    eFlag1 = out[p][4]
                                    p = p + 1
                                    eFlag2 = out[p][4]
                                    #Simple binary operations to calculate the flag for the introns
                                    #The intron flag is calculated by using the left exon flag and the
                                    #right exon flag by using the bits of the binary coded flag
                                    # for further information how they are calculated look up the
                                    # documentation
                                    out[val][4] = ((eFlag1 & 1) << 1) | (eFlag2 >> 1)

                                    seq = (out[val][0], out[val][1], out[val][2], out[val][3]+'/'+str(ecount-1), str(out[val][4]), out[val][5], "\n")
                                    output.write(str('\t').join(seq))


                                elif out[val][6] == "e":
                                    seq = (out[val][0], out[val][1], out[val][2], out[val][3]+'/'+str(ecount), str(out[val][4]), out[val][5], "\n")
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

                if self.geneName:
                    if feature.attr.has_key(self.geneName):
                        gene_name = str(feature.attr[str(self.geneName)])
                    else:
                        error = "Wrong gene Name: "+self.geneName+". Check your annotation file!!"
                        raise KeyError(error)

                    name = name + '@' + gene_name + '@' + attribute
                else:
                    name = name+'@'+attribute
                    
                t = attribute
                start = feature.iv.start
                end = feature.iv.end
                sc_start = 0
                sc_stop = 0
                tss = 0
                sc = []
                ss = []
                ts =[]

                ec = 1
                ic = 1
                ps = []
                ps_cds = []
                fs = []

            elif feature.type == "start_codon":
                tss = feature.iv.start_d

            elif feature.type == "stop_codon":
                sc_start = feature.iv.start_d
                sc_stop = feature.iv.end_d

            elif feature.type == "exon" or feature.type == "exonic":
                giv = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.strand)
                ps.append((feature.iv.start, feature.iv.end))
                tag = 'exon'
                gas[giv] += tag

            elif feature.type == "CDS":
                giv = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.strand)
                ps_cds.append((feature.iv.start, feature.iv.end))
                tag = 'cds'
                gas[giv] += tag

            elif feature.type == "UTR":
                giv = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.strand)
                ps_cds.append((feature.iv.start, feature.iv.end))
                tag = 'utr'
                gas[giv] += tag
                
            if not ts:
                tss = 0
            else:
                tss = int(min(ts))
            if not sc:
                sc_start = 0
            else:
                sc_start = int(max(sc))
            if not ss:
                sc_stop = 0
            else:
                sc_stop = int(min(ss))

        #Write out information for normalization of plots
        chromFooter = OrderedDict(sorted(chromFooter.items(), key=lambda x: (-x[1], x[0])))

        for k in chromFooter:
            output.write("track chr "+str(k)+" "+str(chromFooter[k])+"\n")

        typeFooter = OrderedDict(sorted(typeFooter.items(), key=lambda x: (-x[1], x[0])))

        for k in typeFooter:
            output.write("track type "+str(k)+" "+str(typeFooter[k])+"\n")
            
        if self.geneName:
            output.write("track"+" "+"gene_name"+"\n")
        else:
            output.write("track"+"\n")
        print 'Finished!'
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

    """
    Function for flagging the strands
    """
    def flag_region(self,gene,pos_list,flag):

         #flag regions(utrs and cds).

            #loopove through positions of cds and utrs for each transcript and flag the positions accordingly
            # Flag = 3 means starting and ending position for the regions are same for all transcripts
            #Flag = 1 means in some transcripts ending positions differ
            #Flag = 2 means in some transcripts starting positions differ
            #Flag = 0 means neither starting nor ending position match for any transcripts (unreliable positions)

        for ps in range(0,len(pos_list)):

            if ps < len(pos_list) and pos_list[ps][1] <= gene[0].end:
                #print 'entering flagging loop'

                if pos_list[ps][0] == gene[0].start and pos_list[ps][1] == gene[0].end:

                    flag = 3
                    del pos_list[ps]

                elif gene[0].strand == '-' and pos_list[ps][0] == gene[0].start and pos_list[ps][1] != gene[0].end and flag != 0:

                    flag = 1
                    del pos_list[ps]

                elif gene[0].strand == '+' and pos_list[ps][0] == gene[0].start and pos_list[ps][1] != gene[0].end and flag != 0:

                    flag = 2
                    del pos_list[ps]

                elif gene[0].strand == '-' and pos_list[ps][0] != gene[0].start and pos_list[ps][1] == gene[0].end and flag != 0:

                    flag = 2
                    del pos_list[ps]

                elif gene[0].strand == '+' and pos_list[ps][0] != gene[0].start and pos_list[ps][1] == gene[0].end and flag != 0:

                    flag = 1
                    del pos_list[ps]

                elif pos_list[ps][0] != gene[0].start and pos_list[ps][1] != gene[0].end and flag != 0:

                    flag = 0
                    del pos_list[ps]

            else:
                break

        return pos_list, flag

    """
    Function for separating exons into UTRs and CDSs
    """
    def regions(self, gene, tss, sc_start, sc_stop, typeFooter, name, ec, eflag, fs):

        #separating utr and cds on forward strand
        #for identifying cds and utrs compare starting and ending positions of the region with
        #the start codon start position and end codon start position
        #5utr ending should be less than start position of start codon
        #3utr start position should be greater than or equal to start position of stop codon
        #CDS region should have positions greater than equal to start position of start codon
        #and less then the stop codon start position.

        if gene[0].strand == '+':
            print '+ strand'
            print gene
            print tss,sc_start
            if "utr" in (list(gene[1])) and len(list(gene[1])) == 2:

                if tss >= gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)
                    print '5utr done'

                elif sc_start <= gene[0].start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)
                    print '3utr done'

                elif tss == 0:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)
                    print'2 5utr done'

                    #determining length of CDS
            elif "cds" in (list(gene[1])) and len(list(gene[1])) == 2:
                if gene[0].start >= tss and gene[0].end <= sc_start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)
                    print 'cds condition 1'

                elif tss == 0 and sc_start == 0:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)
                    print 'cds condition 2'

                elif tss > gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)
                    print 'cds condition 3'

                elif sc_start < gene[0].start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)
                    print 'cds condition 4'

                else:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)
                    print 'cds condition 5'
                    #if set contains utr,exon,cds
            elif len(list(gene[1])) > 2:
                if gene[0].start < tss and gene[0].end <= tss:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,tss)
                    print 'overlapping condition 1'

                elif gene[0].start > sc_start :
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)
                    print 'overlapping condition 2'

                elif gene[0].start >= tss and gene[0].end <= sc_start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)
                    print 'overlapping condition 3'

                elif gene[0].start >= tss and gene[0].end >= sc_start:
                    if not typeFooter.has_key("CDS"):
                        typeFooter["CDS"] = (sc_start - gene[0].start)
                        if not typeFooter.has_key("3UTR"):
                            typeFooter["3UTR"] = (gene[0].end - sc_start)
                    else:
                        typeFooter["CDS"] += (gene[0].end - gene[0].start)
                        typeFooter["3UTR"] += (gene[0].end - sc_start)
                    sec_tag = 'CDS'
                    print 'overlapping condition 3'
        else:
            print '- strand'
            print gene
            if "utr" in (list(gene[1])) and len(list(gene[1])) == 2:
                if tss < gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)

                elif sc_start <= gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)

                elif sc_start+1 >= gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,sc_start)

                elif tss == 0:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)

            #determining length of CDS
            elif "cds" in (list(gene[1])) and len(list(gene[1])) == 2:
                if gene[0].end <= tss and gene[0].start >= sc_start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)

                elif tss == 0 and sc_start == 0:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)

                elif tss < gene[0].start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)

                elif sc_start < gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)

                else :
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)

            elif len(list(gene[1])) > 2:
                if tss < gene[0].start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'5UTR',gene[0].start,gene[0].end)

                elif sc_stop > gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)

                elif sc_start+1 <= gene[0].end:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'3UTR',gene[0].start,gene[0].end)

                elif gene[0].end <= tss and gene[0].start >= sc_start:
                    sec_tag, typeFooter = self.get_tag(typeFooter,'CDS',gene[0].start,gene[0].end)
        #Calculate the values which are need for the normalization of the plot
        # append regions to the list
        # if any region has flag None turn it to -1 to differentiate it as lone region.
        if eflag == None:
            eflag = -1
            seq = [gene[0].chrom, str(gene[0].start), str(gene[0].end), name+'@'+sec_tag+'@'+str(ec), eflag, gene[0].strand, "e"]
        else:
            seq = [gene[0].chrom, str(gene[0].start), str(gene[0].end), name+'@'+sec_tag+'@'+str(ec), eflag, gene[0].strand, "e"]
        fs.append(seq)
        return fs


    """
    Function for processing each gene
    """
    def process(self,gene,ps_cds,typeFooter,tss,sc_start,sc_stop,name,ec,ic):
        fs = []
    #confidence info about region.
        eflag = None

    # contains sorted positions for all transcripts
        ps_cds.sort(key=lambda tup: tup[0])

    # for each gene (list of lists)
    # g gene (list) strand info, start and end position
        for g in gene:

        # if set is empty, then it is an intron,
        # increment exon count and append intron data to list
            if not g[1]:
            #ec keeps count of exons
                ec += 1

        #calculate the intron length
                if not typeFooter.has_key("protein_coding_intron"):
                    typeFooter["protein_coding_intron"] = (g[0].end - g[0].start)
                else:
                    typeFooter["protein_coding_intron"] += (g[0].end - g[0].start)

        #append intron info to fs
                seq = [g[0].chrom, str(g[0].start), str(g[0].end), name+'@intron@'+str(ic), 1, g[0].strand, "i"]
                fs.append(seq)
                ic += 1

        #if set is of length 1, i.e region does not contain cds or utr then continue
            elif len(list(g[1])) == 1:
                if not typeFooter.has_key("protein_coding_exon"):
                    typeFooter["protein_coding_exon"] = (g[0].end - g[0].start)
                else:
                    typeFooter["protein_coding_exon"] += (g[0].end - g[0].start)
                seq = [g[0].chrom, str(g[0].start), str(g[0].end), name+'@exon@'+str(ic), -1, g[0].strand, "e"]
                fs.append(seq)
            else:
                #flag the region first
                ps_cds,eflag = self.flag_region(g,ps_cds,eflag)
                fs = self.regions(g, tss, sc_start,sc_stop,typeFooter,name,ec,eflag,fs)
        return fs, ec

    """
    Function for getting region tag
    """
    def get_tag(self,typeFooter,region,start,end):
        if not typeFooter.has_key(region):
            typeFooter[region] = (end - start)
        else:
            typeFooter[region] += (end - start)
        sec_tag = region
        return sec_tag, typeFooter
