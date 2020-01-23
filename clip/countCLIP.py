from builtins import str
from builtins import object

import gzip
import logging
import os
import sys
import tempfile
from decimal import Decimal
from shutil import copyfile

from HTSeq import BED_Reader, GenomicArray, GenomicFeature, GenomicInterval

from .output import Output

class TempBed(object):
    '''
    Temp file object for annotaiton files
    This is a workaround to weird issues with reading gziped bed file
    with workflows using slurm/snakemake
    '''
    def __init__(self,bedfile):
        self.bedfile = bedfile
    
    def __enter__(self):
        self.tmpBed = tempfile.mkstemp()[1]
        copyfile(self.bedfile,self.tmpBed)
        return self.tmpBed
    
    def __exit__(self,exec_type,exec_val,exec_tback):
        os.remove(self.tmpBed)

class countCLIP(object):
    # indices
    __indexGeneID__ = 0 
    __indexGeneSymbol__ = 1
    __indexGeneType__ = 2
    __indexGeneRegion__ = 3
    __indexRegionPos__ = 4
    __indexUID__ = 5
    __indexWindowNumber__ = 6
    # separators
    __dividerRegion__ = "/"
    __dividerName__  = "@"
    # zero padding
    __zeroFillFeature__ = 4
    __zeroFillWindow__ = 6
    # expected number of columns
    __numberColumnsAnnotationWithWindows__ = 7
    # max number of entries to read in for sanity checking
    __iMax__ = 100

    def __init__(self, options):
        self.annotation = options.annotation
        self._nameColCount = 0
        self._isWindowed  = None
        self._decoder = None
        if hasattr(options,'input'):
            self.sites = options.input
        self.output = Output(options.output)
        
        self.fo = self._annotationParser()

        self._annotationSanityCheck()
    
    def _toStr(self,line):
        '''
        helper function
        given a string return it as it is
        '''
        return line
    
    def _byteToStr(self,line):
        '''
        helper function
        given bytes decode to string
        '''
        return line.decode('utf-8')

    def _annotationParser(self):
        if self.annotation.lower().endswith((".gz",".gzip")):
            self._decoder = self._byteToStr
            return gzip.open
        else:
            self._decoder = self._toStr
            return open
        
    def _annotationSanityCheck(self):
        '''
        Helper function
        Check sanity of input annotations
        '''
        nameCountSet = set()
        i = 0
        with TempBed(self.annotation) as tann:
            with self.fo(tann) as _ah:
                for line in _ah:
                    line = self._decoder(line)
                    if line.startswith( "track" ):
                        continue
                    afields = line.strip("\n").split("\t")
                    if len(afields) < 6:
                        raise ValueError("BED file line contains less than 6 fields")
                    if len(afields) > 9:
                        raise ValueError("BED file line contains more than 9 fields")
                    nameCountSet.add(len(afields[3].split(self.__dividerName__)))
                    i+=1
                    if i>=self.__iMax__:
                        break
        nameCountSet = list(nameCountSet)
        if len(nameCountSet)!=1:
            raise ValueError("The 'name' column in {} is incorrectly formatted. The number of 'name' entries must be equal in all rows of the file".format(self.annotation))
        if nameCountSet[0]== self.__numberColumnsAnnotationWithWindows__:
            logging.debug("{} format: sliding window annotation file")
            self._isWindowed = True
            self._nameColCount = self.__numberColumnsAnnotationWithWindows__
        elif nameCountSet[0]== (self.__numberColumnsAnnotationWithWindows__-1):
            logging.debug("{} format: gene feature annotation file")
            self._isWindowed = False
            self._nameColCount = self.__numberColumnsAnnotationWithWindows__-1
        else:
            sw = self.__numberColumnsAnnotationWithWindows__
            nosw = sw-1
            raise ValueError("{} unknown format. There must be either {} or {} fields separated by {} in 'name' column".format(self.annotation,nosw,sw,self.__dividerName__))

    def _splitName(self, name):
        """
        Helper function
        This function splits the name of an htseq-clip processed annotation file
        which can be windowed or not windows and returns the processed values
        """
        # splitting up the name
        try:
            arr = name.split(self.__dividerName__)
        except ValueError as err:
            logging.exception('Error processing annotation name: {}'.format(str(err)))
            logging.exception('Could not split name with divider {} from {}'.format(self.__dividerName__,name))

        if len(arr)!= self._nameColCount:
            raise ValueError("Malformed name entry '{}'. There must be {} fields separated by {}".format(name,self._nameColCount,self.__dividerName__))
            
        geneID = arr[self.__indexGeneID__]
        geneSymbol = arr[self.__indexGeneSymbol__]
        UID = arr[self.__indexUID__]
        geneType = arr[self.__indexGeneType__]
        geneRegion = arr[self.__indexGeneRegion__]
        
        try:
            (genePositionNumber, genePositionTotal) = arr[self.__indexRegionPos__].split(self.__dividerRegion__)
        except ValueError as err:
            logging.exception('Error processing annotation name: {}'.format(str(err)))
            logging.exception('Could not split value with index {} with divider {} from {}'.format(self.__indexRegionPos__, 
                                                                                       self.__dividerRegion__,
                                                                                       arr))
            
        # # empty string if not windowed. 
        # annotationWindowNumber = windowNumber.zfill(self.__zeroFillWindow__) if isWindowed else ""
        # annotationWindowSeparator = "W" if isWindowed else ""
        
        if self._isWindowed:
            windowNumber = arr[self.__indexWindowNumber__]
            return (UID, geneID, geneSymbol, geneType, geneRegion, genePositionNumber, genePositionTotal, windowNumber)
        else:
            return (UID, geneID, geneSymbol, geneType, geneRegion, genePositionNumber, genePositionTotal)
    
    def annotationToIDs(self):
        '''
        given an annotation file convert it to a mapping file
        '''
        colHeader = ['unique_id','chromosome','begin','end','strand','gene_id','gene_name','gene_type','gene_region','Nr_of_region','Total_nr_of_region']
        if self._isWindowed:
            colHeader.append('window_number')
        self.output.write("\t".join(colHeader)+"\n")
        for anno in BED_Reader(self.annotation):
            # annData = (anno.iv.chrom,str(anno.iv.start),str(anno.iv.end),anno.iv.strand)
            nameAnn = self._splitName(anno.name)
            self.output.write("\t".join((nameAnn[0],anno.iv.chrom,str(anno.iv.start),str(anno.iv.end),anno.iv.strand)+nameAnn[1:])+"\n")
        self.output.close()
    
    def _countSites(self,strandedCounting):
        '''
        Helper function, count sites and return as an object of GenomicArray
        '''
        # data structures to hold sites of interest (crosslink sites, deletion sites, insertion sites)
        sitesga = GenomicArray(chroms="auto", stranded=strandedCounting, typecode='i', storage='step')
        
        if self.sites.lower().endswith((".gz",".gzip")):
            sfo = gzip.open
            decoder = self._byteToStr
        else:
            sfo = open
            decoder = self._toStr
        with TempBed(self.sites) as tb:
            with sfo(tb) as _sh:
                for line in _sh:
                    line = decoder(line)
                    if line.startswith("track"):
                        continue
                    sFields = line.strip('\n').split("\t")
                    if len(sFields) < 3:
                        raise ValueError("BED file line contains less than 3 fields")
                    if len(sFields) > 9:
                        raise ValueError("BED file line contains more than 9 fields")
                    iv = GenomicInterval( sFields[0], int(sFields[1]), int(sFields[2]), sFields[5] if len(sFields) > 5 else "." )
                    # add the sites to the data structure
                    sitesga[iv]+=1
        return sitesga

    def count(self, strandedCounting = True):
        sites_ga =  self._countSites(strandedCounting)
        # for site in BED_Reader(self.sites):
        #     sites_ga[ site.iv ] += 1
        # write column headers for the output
        # sys.stderr.write("{}".format(self._isWindowed))
        if self._isWindowed:
            colHeader = ['unique_id','window_number','window_length','crosslink_count_total','crosslink_count_position_nr','crosslink_count_position_max','crosslink_density']
        else:
            colHeader = ['unique_id','region_length','crosslink_count_total','crosslink_count_position_nr','crosslink_count_position_max','crosslink_density']
        self.output.write("\t".join(colHeader)+'\n')
        # go through annotation and count the sites in the specific annotation
        # the sites_ga will return steps, which are already processed ranges with the amount of sites
        # e.g. an annotation element of length 10 is split up in elements of length 5 with 0 counts, length 2 with 1 count, 
        # length 2 with 3 counts, and length 1 with 0 counts. in this example the total sum of counts would be 
        # 5*0 + 2*1 + 2*3 + 1*0 which would be 7. 
        # length of positions occupied would be 2 + 2, therefore 4. 
        # the maximum counts, the height would be 3.
        # this is a workaround to annotation file crashing on snakemake runs
        with TempBed(self.annotation) as ta:
            with self.fo(ta) as _ah: # annotation handler
                # these lines are a work around to odd splitting behavior
                for line in _ah:
                    line = self._decoder(line)
                    if line.startswith( "track" ):
                        continue
                    bfields = line.strip('\n').split("\t")
                    if len(bfields) < 3:
                        raise ValueError("BED file line contains less than 3 fields")
                    if len(bfields) > 9:
                        raise ValueError("BED file line contains more than 9 fields")
                    iv = GenomicInterval( bfields[0], int(bfields[1]), int(bfields[2]), bfields[5] if len(bfields) > 5 else "." )
                    anno = GenomicFeature( bfields[3] if len(bfields) > 3 else "unnamed", "BED line", iv )
                    anno.score = float( bfields[4] ) if len(bfields) > 4 else None
                    anno.thick = GenomicInterval( iv.chrom, int( bfields[6] ), int( bfields[7] ), iv.strand ) if len(bfields) > 7 else None
                    anno.itemRgb = [ int(a) for a in bfields[8].split(",") ]  if len(bfields) > 8 else None
                    # resume working with annotation
                    total_counts = 0   # variable for storing the total sum of site counts
                    max_counts = 0     # maximum site height, maximum site counts
                    total_length = anno.iv.length  # length of the annotation element
                    occupied_positions_length = 0  # length of the positions with site counts
                    nameAnns = self._splitName(anno.name)

                    for giv, count in sites_ga[ anno.iv ].steps():
                        if count == 0:
                            continue
                        # add the number of sites to the total counts 
                        total_counts += count * giv.length
                        
                        # store maximal count 
                        if count > max_counts:
                            max_counts = count
                        
                        # add length if sites are present in this interval
                        occupied_positions_length += giv.length
                
                    if total_counts == 0:
                        continue
                    elif total_counts < occupied_positions_length:
                        raise Exception("Total site counts {} can not be smaller than occupied_positions_length {} at {}\n{}".format(total_counts, occupied_positions_length, anno.name, anno.iv)) 
                    # site density is the number of occupied positions divided the total length of the annotation element
                    site_density = round(float(occupied_positions_length) / float(total_length),5)
                    # print count data
                    countString = [str(total_length),str(total_counts),str(occupied_positions_length),str(max_counts),str(site_density)]
                    if self._isWindowed:
                        self.output.write("\t".join([nameAnns[0],nameAnns[-1]]+countString)+'\n')
                    else:
                        self.output.write("\t".join([nameAnns[0]]+countString)+'\n')
        self.output.close()
