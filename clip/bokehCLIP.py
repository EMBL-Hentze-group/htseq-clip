# --------------------------------------------------
# bedCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: December 2015
# --------------------------------------------------

import os, gzip
from decimal import *
import numpy as np
from math import log10
from collections import OrderedDict
from random import randint
from scipy import log2

try:
    from bokeh.plotting import figure, output_file, show, save, vplot, hplot
    from bokeh.models import Range1d
    from bokeh.charts import Histogram, Line, Scatter, BoxPlot, output_file, show, save, vplot, hplot, Bar, defaults
    from bokeh.models.widgets import DataTable, TableColumn, Panel, Tabs
    from bokeh.models import ColumnDataSource
    from bokeh.charts.attributes import CatAttr
    
except Exception:
    print "Please install the bokeh framework e.g. like this"
    print "pip install bokeh"
    print "pip install bokeh --user"
    os._exit(1)

class bokehCLIP:
    
    input = ''
    fCompare = ''
    output = ''
    
    def __init__(self, options):
        
        if hasattr(options, 'input'):
            self.input = options.input
            
        if hasattr(options, 'compare'):
            self.fCompare = options.compare
            
        if hasattr(options, 'output'):
            self.output = options.output
        
    #===================================================================================
    #===================================================================================
    def plot_read(self):
        
        if self.input.endswith(".gz"):
            almnt_file = gzip.open(self.input, 'r') 
        else:        
            almnt_file = open(self.input, 'r')
         
        heading = self.output
        
        heading = heading.replace('\\','/')
        heading = heading.split(".html")
        path    = heading[0]
        heading = heading[0].split("/")
        heading = heading[-1]
        
        readLength  = {}
        duplicates  = {}
        chromosomes = {}
        strands     = {}
        
        rlScatter = []
        dpScatter = []
        
        #Extracting data out of file
        for line in almnt_file:
            
            line = line.split("\n")
            line = line[0].split("\t")
            rl = line[3].split("|")
            rl = rl[1]
                   
            if not readLength.has_key(int(rl)):
                readLength[int(rl)] = 1
            else:
                readLength[int(rl)] += 1

            if not duplicates.has_key(int(line[4])):
                duplicates[int(line[4])] = 1
            else:
                duplicates[int(line[4])] += 1
            
            if not chromosomes.has_key(line[0]):
                chromosomes[line[0]] = 1
            else:
                chromosomes[line[0]] += 1
                
            if not strands.has_key(line[5]):
                strands[line[5]] = 1
            else:
                strands[line[5]] += 1  
                
            rlScatter.append(int(rl)) 
            dpScatter.append(int(line[4]))
        
           
#         rlScatterData = []
#         dpScatterData = []
#         
#         #Generate Data for scatter plot
#         if len(rlScatter) > 50000:    
#             sd = [randint(0,len(rlScatter)-1) for p in range(0,50000)]
#         
#             for i in sd:
#                 rlScatterData.append(rlScatter[i])
#                 dpScatterData.append(dpScatter[i])
#         else:
#             rlScatterData = rlScatter
#             dpScatterData = dpScatter
               
        #Data dicts              
        rlData     = {}
        dpData     = {}
        chromData  = {}
        strandData = {}
        
        #Data Arrays
        rlKeys   = []
        rlCounts = []
        
        dpKeys   = []
        dpCounts = []
           
        chromKeys   = []
        chromCounts = []
        
        strandKeys   = []
        strandCounts = []
        
        #Preparing data for plotting
        for r in readLength:
            rlKeys.append(r)
            rlCounts.append(readLength[r])
            
        for d in duplicates:
            dpKeys.append(d)
            dpCounts.append(duplicates[d])
            
        chromosomes = OrderedDict(sorted(chromosomes.items(), key=lambda x: (-x[1], x[0])))
              
        for c in chromosomes:
            chromKeys.append(c)
            chromCounts.append(chromosomes[c])
            
        for s in strands:
            strandKeys.append(s)
            strandCounts.append(strands[s])
            
        #Generate dataframes for plots    
        rlData = {
            'rlKeys' : rlKeys,
            'rlCounts' : rlCounts                         
        }
        
        dpData = {
            'dpKeys' : dpKeys,
            'dpCounts' : dpCounts        
        }
        
        chromData = {
            'chromKeys' : chromKeys,
            'chromCounts' : chromCounts
        }
        
        strandData = {
            'strandKeys' : strandKeys,
            'strandCounts' : strandCounts
        }
        
#         scatterData = {
#             'rl' : rlScatter,
#             'dp' : dpScatter
#         }
              
        #Toolbar for each plot       
        TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save"  
               
        #Plotting
        output_file(path+".html")
        
        #Read length         
        rlP = Bar(rlData, title=heading+"_readLength", values='rlCounts', label='rlKeys', tools=TOOLS, xlabel="Read length", ylabel="Number of reads")
        
        rlData_table = DataTable(source=ColumnDataSource(rlData), columns=[TableColumn(field="rlKeys", title="Read Length"), TableColumn(field="rlCounts", title="Number of reads")], width=600, height=600)
        
        rlPlot = vplot(
                hplot(rlP, rlData_table)
        )
        tab1 = Panel(child=rlPlot, title="Read length")
        
        rlDataOutput = open(path+"_rlData.txt", "w")
        rlDataOutput.write("Read Length\tNumber of Reads\n")
        for i in range(0, len(rlKeys)):
            rlDataOutput.write(str(rlKeys[i]) + "\t" + str(rlCounts[i]) + "\n")
        
        rlDataOutput.close()
         
        #Duplicates         
        dpData_table = DataTable(source=ColumnDataSource(dpData), columns=[TableColumn(field="dpKeys", title="Duplicates"),TableColumn(field="dpCounts", title="Number of duplicates")], width=600, height=600)
        
        tab2 = Panel(child=dpData_table, title="Duplicates")
        
        dpDataOutput = open(path+"_dpData.txt", "w")
        dpDataOutput.write("Duplicates\tNumber of duplicates\n")
        for i in range(0, len(dpKeys)):
            dpDataOutput.write(str(dpKeys[i]) + "\t" + str(dpCounts[i]) + "\n")
        
        dpDataOutput.close()
               
        #Chromosomes
        chromP = Bar(chromData, title=heading+"_chromosomes", values='chromCounts', label=CatAttr(columns=['chromKeys'], sort=False),tools=TOOLS, xlabel="Chromosomes", ylabel="Count") 
        chromData_table = DataTable(source=ColumnDataSource(chromData), columns=[TableColumn(field="chromKeys", title="Chromosomes"),TableColumn(field="chromCounts", title="Count")], width=600, height=600)

        chromPlot = vplot(
                hplot(chromP, chromData_table)
        )
        
        tab3 = Panel(child=chromPlot, title="Chromosomes")
        
        chromDataOutput = open(path+"_chromData.txt", "w")
        chromDataOutput.write("Chromosomes\tCounts\n")
        for i in range(0, len(chromKeys)):
            chromDataOutput.write(str(chromKeys[i]) + "\t" + str(chromCounts[i]) + "\n")
        
        chromDataOutput.close()
           
        #Strands
        strandP = Bar(strandData, title=heading+"_strands", values='strandCounts', label='strandKeys',tools=TOOLS, xlabel="Strands", ylabel="Count")
        strandData_table = DataTable(source=ColumnDataSource(strandData), columns=[TableColumn(field="strandKeys", title="Strands"),TableColumn(field="strandCounts", title="Count")], width=600, height=600)
        
        strandPlot = vplot(
                hplot(strandP, strandData_table)
        )
        
        tab4 = Panel(child=strandPlot, title="Strands")
        
        strandDataOutput = open(path+"_strandData.txt", "w")
        strandDataOutput.write("Strands\tCounts\n")
        for i in range(0, len(strandKeys)):
            strandDataOutput.write(str(strandKeys[i]) + "\t" + str(strandCounts[i]) + "\n")
        
        strandDataOutput.close()
        
        #bxPlot = BoxPlot(scatterData, label='rl', values='dp', title="Duplicates vs Read length (N = 50.000)", tools="resize, save", xlabel="Read length", ylabel="Duplicates")
        #tab5 = Panel(child=bxPlot, title="Read length vs Duplicates")
        
        tabs = Tabs(tabs=[ tab1, tab2, tab3, tab4])
        
        save(tabs)
        
    #===================================================================================
    #===================================================================================
    def plot_count(self):
        
        if self.input.endswith(".gz"):
            almnt_file = gzip.open(self.input, 'r') 
        else:        
            almnt_file = open(self.input, 'r')
        
        heading = self.output
        
        heading = heading.replace('\\','/')
        heading = heading.split(".html")
        path    = heading[0]
        heading = heading[0].split("/")
        heading = heading[-1]
        
        chrNorm     = {}
        featureNorm = {}
        
        counts      = {}
        countKeys   = []
        countCounts = []

        
        density = {}
        
        duplicates = {}
        dpKeys     = []
        dpCounts   = []
        
        
        chromosomes  = {}   
        chromKeys    = []
        chromCounts  = []
        chromNKeys   = []
        chromNCounts = []
        
        strands      = {}
        strandKeys   = []
        strandCounts = []
        
        types     = {}
        type      = []
        typeCount = []
        
        regions     = {}
        region      = []
        regionCount = [] 
        
        #duplicate/count quality control
        dcQC      = {}
        dcQCKeys  = []
        dcQCCount = []
        
        dupPerType = {}
        dptKeys    = []
        dptCount   = []
        dptNKeys   = []
        dptNCount  = []

        
        countsPerType = {}
        cptKeys       = []
        cptCount      = []
        cptNKeys      = []
        cptNCount     = []
    
        for line in almnt_file:
            
            if line.startswith("track"):
                line = line.split("\n")
                line = line[0].split(" ")
                
                if line[1] == "chr":
                    chrNorm[line[2]] = int(line[3])
                elif line[1] == "type":
                    featureNorm[line[2]] = int(line[3])
                else:
                    error = "Unknown track annotation: "+line[1]+". Check your data!!"
                    raise ValueError(error)
            elif "intergenic" in line:
                line = line.split("\n")
                line = line[0].split("\t")
                
                if not countsPerType.has_key("intergenic"):
                    countsPerType["intergenic"] = int(line[11])
                else:
                    countsPerType["intergenic"] += int(line[11])        
            else:
                       
                line = line.split("\n")
                line = line[0].split("\t")
                
                if not chromosomes.has_key(line[0]):
                    chromosomes[line[0]] = 1
                else:
                    chromosomes[line[0]] += 1
                    
                if not counts.has_key(int(line[11])):
                    counts[int(line[11])] = 1
                else:
                    counts[int(line[11])] += 1
                    
                if not density.has_key(float(line[14])):
                    density[float(line[14])] = 1
                else:
                    density[float(line[14])] += 1
                    
                if not duplicates.has_key(int(line[15])):
                    duplicates[int(line[15])] = 1
                else:
                    duplicates[int(line[15])] += 1
                    
                if line[6] == "exon" and line[9] == "protein_coding":
                    line[9] = "protein_coding_exon"
                    if not types.has_key(line[9]):
                        types[line[9]] = 1
                    else:
                        types[line[9]] += 1  
                        
                    if not dcQC.has_key(line[9]):
                        dcQC[line[9]] = [int(line[15]), int(line[11])]
                    else:
                        dcQC[line[9]][0] += int(line[15])
                        dcQC[line[9]][1] += int(line[11])
                        
                        
                elif line[6] == "intron" and line[9] == "protein_coding":
                    line[9] = "protein_coding_intron"
                    if not types.has_key(line[9]):
                        types[line[9]] = 1
                    else:
                        types[line[9]] += 1   
                        
                    if not dcQC.has_key(line[9]):
                        dcQC[line[9]] = [int(line[15]), int(line[11])]
                    else:
                        dcQC[line[9]][0] += int(line[15])
                        dcQC[line[9]][1] += int(line[11])    
                            
                else:
                    if not types.has_key(line[9]):
                        types[line[9]] = 1
                    else:
                        types[line[9]] += 1
                        
                    if not dcQC.has_key(line[9]):
                        dcQC[line[9]] = [int(line[15]), int(line[11])]
                    else:
                        dcQC[line[9]][0] += int(line[15])
                        dcQC[line[9]][1] += int(line[11])
                        
                if not regions.has_key(line[6]):
                    regions[line[6]] = 1
                else:
                    regions[line[6]] += 1
                    
                if not strands.has_key(line[5]):
                    strands[line[5]] = 1
                else:
                    strands[line[5]] += 1
                
        for k in dcQC:
            dupPerType[k] = dcQC[k][0]
            countsPerType[k] = dcQC[k][1]
            dcQC[k] = dcQC[k][0]/float(dcQC[k][1])
            
        chrN = {}
        cn = {}
        dn = {}
            
        for k in chromosomes:
            chrN[k] = (chromosomes[k]/ float(chrNorm[k]))
            
        for k in countsPerType:
            if not k == "intergenic":
                cn[k] = (countsPerType[k]/ float(featureNorm[k]))
            
        for k in dupPerType:
            dn[k] = (dupPerType[k] / float(featureNorm[k]))
            
        #Order the dicts for plotting by highest value of key
        counts      = OrderedDict(sorted(counts.items(), key=lambda x: (-x[1], x[0])))
        chromosomes = OrderedDict(sorted(chromosomes.items(), key=lambda x: (-x[1], x[0])))
        duplicates  = OrderedDict(sorted(duplicates.items(), key=lambda x: (-x[1], x[0])))
        regions     = OrderedDict(sorted(regions.items(), key=lambda x: (-x[1], x[0])))
        types       = OrderedDict(sorted(types.items(), key=lambda x: (-x[1], x[0])))
        dcQC    = OrderedDict(sorted(dcQC.items(), key=lambda x: (-x[1], x[0])))
        dupPerType  = OrderedDict(sorted(dupPerType.items(), key=lambda x: (-x[1], x[0])))
        countsPerType  = OrderedDict(sorted(countsPerType.items(), key=lambda x: (-x[1], x[0])))
        chrN  = OrderedDict(sorted(chrN.items(), key=lambda x: (-x[1], x[0])))
        cn  = OrderedDict(sorted(cn.items(), key=lambda x: (-x[1], x[0])))
        dn  = OrderedDict(sorted(dn.items(), key=lambda x: (-x[1], x[0])))
        
        ##########################
        #Prepare data for plotting
        ##########################
        
        #--------------------------------------------------
        for k in counts:
            countKeys.append(k)
            countCounts.append(counts[k])
            
        countData = {
            'cKeys' : countKeys,
            'cCounts' : countCounts
        }
        #--------------------------------------------------
            
        #--------------------------------------------------          
        for k in duplicates:
            dpKeys.append(k)
            dpCounts.append(duplicates[k])
            
        dpData = {
            'dpKeys' : dpKeys,
            'dpCounts' : dpCounts        
        }  
        #--------------------------------------------------
        
        #--------------------------------------------------      
        for k in chromosomes:
            chromKeys.append(k)
            chromCounts.append(chromosomes[k])
            
        chromData = {
            'chromKeys' : chromKeys,
            'chromCounts' : chromCounts
        }
        
        for k in chrN:
            chromNKeys.append(k)
            chromNCounts.append(log2(chrN[k]+1))
        
        chromNData = {
            'chromKeys' : chromNKeys,
            'chromCounts' : chromNCounts                   
        }
        #--------------------------------------------------
        
        #--------------------------------------------------    
        for k in strands:
            strandKeys.append(k)
            strandCounts.append(strands[k])
            
        strandData = {
            'strandKeys' : strandKeys,
            'strandCounts' : strandCounts
        }
        #--------------------------------------------------
        
        #--------------------------------------------------  
        for k in regions:
            region.append(k)
            regionCount.append(regions[k])
            
        regionData = {
            'region' : region,
            'regionCounts' : regionCount            
        }          
        #--------------------------------------------------
        
        #--------------------------------------------------   
        for k in types:
            type.append(k)
            typeCount.append(types[k])
        
        typeData = {
            'type' : type,
            'typeCounts' : typeCount       
        }
            
        for k in dcQC:
            dcQCKeys.append(k)
            dcQCCount.append(dcQC[k])
            
        dcQCData = {
            'dcQCKeys' : dcQCKeys,
            'dcQCCount' : dcQCCount
        }
            
        for k in dupPerType:
            dptKeys.append(k)
            dptCount.append(dupPerType[k])
            
        dptData = {
            'dptKeys' : dptKeys,
            'dptCount' : dptCount                 
        }
        
        nDSum = 0.0
        
        for k in dn:
            dptNKeys.append(k)
            dptNCount.append(dn[k])
            nDSum += dn[k]
            
        #Normalize the data by Xi / Sum(Xi)
        for i in range(0, len(dptNCount)-1):
            dptNCount[i] = dptNCount[i] / nDSum
            
        dptNData = {
            'dptKeys' : dptNKeys,
            'dptCount' : dptNCount                 
        }
            
        for k in countsPerType:
            cptKeys.append(k)
            cptCount.append(countsPerType[k])
            
        cptData = {
            'cptKeys' : cptKeys,
            'cptCount' : cptCount                 
        }
        
        nCSum = 0.0
        
        for k in cn:
            cptNKeys.append(k)
            cptNCount.append(cn[k])
            nCSum += cn[k]
            
        #Normalize the data by Xi / Sum(Xi)
        for i in range(0, len(cptNCount)-1):
            cptNCount[i] = cptNCount[i] / nCSum
         
        cptNData = {
            'cptKeys' : cptNKeys,
            'cptCount' : cptNCount                 
        }
        
        #--------------------------------------------------
                 
        #Toolbar for each plot       
        TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save"  
               
        #Plotting
        output_file(path+".html")
        
        #Types
        #--------------------------------------------------
        cptP = Bar(cptData, title=heading+"_Counts_per_Type", values='cptCount', label=CatAttr(columns=['cptKeys'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Count")
        cptData_table = DataTable(source=ColumnDataSource(cptData), columns=[TableColumn(field="cptKeys", title="Types"),TableColumn(field="cptCount", title="Count")], width=600, height=600)
        
        cptPlot = vplot(
                hplot(cptP, cptData_table)
        )
        
        tab1 = Panel(child=cptPlot, title="Counts per type")
        
        cptDataOutput = open(path+"_cptData.txt", "w")
        cptDataOutput.write("Types\tCount\n")
        for i in range(0, len(cptKeys)):
            cptDataOutput.write(str(cptKeys[i]) + "\t" + str(cptCount[i]) + "\n")
        
        cptDataOutput.close()
        
        cptNP = Bar(cptNData, title=heading+"_Counts_per_Type", values='cptCount', label=CatAttr(columns=['cptKeys'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Normalized count")
        cptNData_table = DataTable(source=ColumnDataSource(cptNData), columns=[TableColumn(field="cptKeys", title="Types"),TableColumn(field="cptCount", title="Normalized count")], width=600, height=600)
        
        cptNPlot = vplot(
                hplot(cptNP, cptNData_table)
        )
        
        tab2 = Panel(child=cptNPlot, title="Counts per type normalized")
        
        cptNDataOutput = open(path+"_cptNData.txt", "w")
        cptNDataOutput.write("Types\tNormalized count\n")
        for i in range(0, len(cptNKeys)):
            cptNDataOutput.write(str(cptNKeys[i]) + "\t" + str(cptNCount[i]) + "\n")
        
        cptNDataOutput.close()
        
        dptP = Bar(dptData, title=heading+"_Duplicates_per_Type", values='dptCount', label=CatAttr(columns=['dptKeys'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Count")
        dptData_table = DataTable(source=ColumnDataSource(dptData), columns=[TableColumn(field="dptKeys", title="Types"),TableColumn(field="dptCount", title="Number of duplicates")], width=600, height=600)
        
        dptPlot = vplot(
                hplot(dptP, dptData_table)
        )
        
        tab3 = Panel(child=dptPlot, title="Duplicates per type")
        
        dptDataOutput = open(path+"_dptData.txt", "w")
        dptDataOutput.write("Types\tNumber of duplicates\n")
        for i in range(0, len(dptKeys)):
            dptDataOutput.write(str(dptKeys[i]) + "\t" + str(dptCount[i]) + "\n")
        
        dptDataOutput.close()
        
        dptNP = Bar(dptNData, title=heading+"_Duplicates_per_Type", values='dptCount', label=CatAttr(columns=['dptKeys'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Normalized count")
        dptNData_table = DataTable(source=ColumnDataSource(dptNData), columns=[TableColumn(field="dptKeys", title="Types"),TableColumn(field="dptCount", title="Duplicates normalized")], width=600, height=600)
        
        dptNPlot = vplot(
                hplot(dptNP, dptNData_table)
        )
        
        tab4 = Panel(child=dptNPlot, title="Duplicates per type normalized")
        
        dptNDataOutput = open(path+"_dptNData.txt", "w")
        dptNDataOutput.write("Types\tDuplicates normalized\n")
        for i in range(0, len(dptNKeys)):
            dptNDataOutput.write(str(dptNKeys[i]) + "\t" + str(dptNCount[i]) + "\n")
        
        dptNDataOutput.close()
           
        dcQCP = Bar(dcQCData, title=heading+"_Sum(Duplicates)/Sum(Counts)", values='dcQCCount', label=CatAttr(columns=['dcQCKeys'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Sum(Duplicates)/Sum(Counts)")
        dcQCData_table = DataTable(source=ColumnDataSource(dcQCData), columns=[TableColumn(field="dcQCKeys", title="Types"),TableColumn(field="dcQCCount", title="Sum(Duplicates)/Sum(Counts)")], width=600, height=600)
        
        dcQCPlot = vplot(
                hplot(dcQCP, dcQCData_table)
        )
        
        tab5 = Panel(child=dcQCPlot, title="Sum(Duplicates)/Sum(Counts)")
        
        dcQCDataOutput = open(path+"_dcQCData.txt", "w")
        dcQCDataOutput.write("Types\tSum(Duplicates)/Sum(Counts)\n")
        for i in range(0, len(dcQCKeys)):
            dcQCDataOutput.write(str(dcQCKeys[i]) + "\t" + str(dcQCCount[i]) + "\n")
        
        dcQCDataOutput.close()
        
        
        typeP = Bar(typeData, title=heading+"_types", values='typeCounts', label=CatAttr(columns=['type'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Count")
        typeData_table = DataTable(source=ColumnDataSource(typeData), columns=[TableColumn(field="type", title="Types"),TableColumn(field="typeCounts", title="Number of types")], width=600, height=600)
         
        tPlot = vplot(
                hplot(typeP, typeData_table),
        )
        
        tab6 = Panel(child=tPlot, title="Total amount of types per feature")  
        
        typeDataOutput = open(path+"_typeData.txt", "w")
        typeDataOutput.write("Types\tNumber of types\n")
        for i in range(0, len(type)):
            typeDataOutput.write(str(type[i]) + "\t" + str(typeCount[i]) + "\n")
        
        typeDataOutput.close()
        #--------------------------------------------------
        
        #Chromosomes
        #--------------------------------------------------
        chromP = Bar(chromData, title=heading+"_chromosomes", values='chromCounts', label=CatAttr(columns=['chromKeys'], sort=False),tools=TOOLS, xlabel="Chromosomes", ylabel="Count") 
        chromData_table = DataTable(source=ColumnDataSource(chromData), columns=[TableColumn(field="chromKeys", title="Chromosomes"),TableColumn(field="chromCounts", title="Count")], width=600, height=600)

        chromPlot = vplot(
                hplot(chromP, chromData_table)
        )
        
        tab7 = Panel(child=chromPlot, title="Chromosomes")
        
        chromDataOutput = open(path+"_chromData.txt", "w")
        chromDataOutput.write("Chromosomes\tCounts\n")
        for i in range(0, len(chromKeys)):
            chromDataOutput.write(str(chromKeys[i]) + "\t" + str(chromCounts[i]) + "\n")
        
        chromDataOutput.close()
        
        chromNP = Bar(chromNData, title=heading+"_chromosomes", values='chromCounts', label=CatAttr(columns=['chromKeys'], sort=False),tools=TOOLS, xlabel="Chromosomes", ylabel="Normalized count") 
        chromNData_table = DataTable(source=ColumnDataSource(chromNData), columns=[TableColumn(field="chromKeys", title="Chromosomes"),TableColumn(field="chromCounts", title="Normalized count")], width=600, height=600)

        chromNPlot = vplot(
                hplot(chromNP, chromNData_table)
        )
        
        tab8 = Panel(child=chromNPlot, title="Chromosomes normalized")
        
        chromDataOutput = open(path+"_chromNData.txt", "w")
        chromDataOutput.write("Chromosomes\tCounts\n")
        for i in range(0, len(chromKeys)):
            chromDataOutput.write(str(chromNKeys[i]) + "\t" + str(chromNCounts[i]) + "\n")
        
        chromDataOutput.close()
        
        #--------------------------------------------------
           
        #Strands
        #--------------------------------------------------
        strandP = Bar(strandData, title=heading+"_strands", values='strandCounts', label='strandKeys',tools=TOOLS, xlabel="Strands", ylabel="Count")
        strandData_table = DataTable(source=ColumnDataSource(strandData), columns=[TableColumn(field="strandKeys", title="Strands"),TableColumn(field="strandCounts", title="Count")], width=600, height=600)
        
        strandPlot = vplot(
                hplot(strandP, strandData_table)
        )
        
        tab9 = Panel(child=strandPlot, title="Strands")
        
        strandDataOutput = open(path+"_strandData.txt", "w")
        strandDataOutput.write("Srands\tCounts\n")
        for i in range(0, len(strandKeys)):
            strandDataOutput.write(str(strandKeys[i]) + "\t" + str(strandCounts[i]) + "\n")
        
        strandDataOutput.close()
        #--------------------------------------------------
        
        #Regions
        #--------------------------------------------------
        regionP = Bar(regionData, title=heading+"_regions", values='regionCounts', label=CatAttr(columns=['region'], sort=False),tools=TOOLS, xlabel="Regions", ylabel="Count")
        regionData_table = DataTable(source=ColumnDataSource(regionData), columns=[TableColumn(field="region", title="Regions"),TableColumn(field="regionCounts", title="Number of regions")], width=600, height=600)
        
        rPlot = vplot(
                hplot(regionP, regionData_table)
        )
        
        regionTab = Panel(child=rPlot, title="Regions")
        
        regionDataOutput = open(path+"_regionData.txt", "w")
        regionDataOutput.write("Regions\tNumber of Regions\n")
        for i in range(0, len(region)):
            regionDataOutput.write(str(region[i]) + "\t" + str(regionCount[i]) + "\n")
        
        regionDataOutput.close()
        #--------------------------------------------------
        
        #Counts
        cData_table = DataTable(source=ColumnDataSource(countData), columns=[TableColumn(field="cKeys", title="Counts"),TableColumn(field="cCounts", title="Number of counts")], width=600, height=600)
         
        tab10 = Panel(child=cData_table, title="Counts")
        
        countDataOutput = open(path+"_countData.txt", "w")
        countDataOutput.write("Counts\tNumber of counts\n")
        for i in range(0, len(countKeys)):
            countDataOutput.write(str(countKeys[i]) + "\t" + str(countCounts[i]) + "\n")
        
        countDataOutput.close()
           
        #Duplicates         
        dpData_table = DataTable(source=ColumnDataSource(dpData), columns=[TableColumn(field="dpKeys", title="Duplicates"),TableColumn(field="dpCounts", title="Number of duplicates")], width=600, height=600)
         
        tab11 = Panel(child=dpData_table, title="Duplicates")
        
        dupDataOutput = open(path+"_dupData.txt", "w")
        dupDataOutput.write("Duplicates\tNumber of duplicates\n")
        for i in range(0, len(dpKeys)):
            dupDataOutput.write(str(dpKeys[i]) + "\t" + str(dpCounts[i]) + "\n")
        
        dupDataOutput.close()
        #--------------------------------------------------
        
        
        tTabs = Tabs(tabs=[tab1, tab2, tab3, tab4, tab5, tab6])
        tTab = Panel(child=tTabs, title="Types")
        
        chrTabs = Tabs(tabs=[tab7, tab8])
        chrTab = Panel(child=chrTabs, title="Chromosomes")
     
        tabs = Tabs(tabs=[ tTab, chrTab, tab9, regionTab, tab10, tab11])
        
        save(tabs) 
        
    #===================================================================================
    #===================================================================================
    def plot_comparison(self):
        
        if self.input.endswith(".gz"):
            replicate1 = gzip.open(self.input, 'r') 
        else:        
            replicate1 = open(self.input, 'r')
            
        if self.fCompare.endswith(".gz"):
            replicate2 = gzip.open(self.fCompare, 'r') 
        else:        
            replicate2 = open(self.fCompare, 'r')
            
        headingR1 = self.input
        headingR1 = headingR1.replace('\\','/')
        headingR1 = headingR1.split('/')
        headingR1 = headingR1[-1]
             
        headingR2 = self.fCompare
        headingR2 = headingR2.replace('\\','/')
        headingR2 = headingR2.split('/')
        headingR2 = headingR2[-1]
         
        heading = self.output
        heading = heading.replace('\\','/')
        heading = heading.split(".html")
        path = heading[0]
        heading = heading[0].split("/")
        heading = heading[-1]
        
        r1Density = {}
        r2Density = {}
        
        r1Counts = {}
        r2Counts = {}
        
        r1ScatterCount = []
        r2ScatterCount = []
        
        r1ScatterDensity = []
        r2ScatterDensity = []
        
        #Go through replicate files and
        #build up dictionaries
        for line in replicate1:
            
            if line.startswith("track"):
                continue
            
            line = line.split("\n")
            line = line[0].split("\t")
             
            key = "c:"+str(line[0])+"s:"+str(line[1])+"e:"+str(line[2])+"id:"+str(line[3])
            r1Counts[key]  = int(line[11])
            r1Density[key] = float(line[14])
            
        for line in replicate2:
            
            if line.startswith("track"):
                continue
            
            line = line.split("\n")
            line = line[0].split("\t")
             
            key = "c:"+str(line[0])+"s:"+str(line[1])+"e:"+str(line[2])+"id:"+str(line[3])
            r2Counts[key]  = int(line[11])
            r2Density[key] = float(line[14])
        
        #Generate arrays for scatter plotting
        
        #Counts
        for key in r1Counts:
            if key in r2Counts:
                r1ScatterCount.append(log10(r1Counts[key]+1))
                r2ScatterCount.append(log10(r2Counts[key]+1))
            else:
                r1ScatterCount.append(log10(r1Counts[key]+1))
                r2ScatterCount.append(log10(0+1))
        
        for key in r2Counts:
            if not key in r1Counts:
                r1ScatterCount.append(log10(0+1))
                r2ScatterCount.append(log10(r2Counts[key]+1))
                
        #Density
        for key in r1Density:
            if key in r2Density:
                r1ScatterDensity.append(log10(r1Density[key]+1))
                r2ScatterDensity.append(log10(r2Density[key]+1))
            else:
                r1ScatterDensity.append(log10(r1Density[key]+1))
                r2ScatterDensity.append(log10(0.0+1))
        
        for key in r2Density:
            if not key in r1Density:
                r1ScatterDensity.append(log10(0.0+1.0))
                r2ScatterDensity.append(log10(r2Density[key]+1))
        
        r1SC = []
        r2SC = []
                
        r1SD = []
        r2SD = []
                
        if len(r1ScatterCount) > 50000:    
            sd = [randint(0,len(r1ScatterCount)-1) for p in range(0,50000)]
        
            for i in sd:
                r1SC.append(r1ScatterCount[i])
                r2SC.append(r2ScatterCount[i])
                
                r1SD.append(r1ScatterDensity[i])
                r2SD.append(r2ScatterDensity[i])
            
        else:
            r1SC = r1ScatterCount
            r2SC = r2ScatterCount
                
            r1SD = r1ScatterDensity
            r2SD = r2ScatterDensity      
            
        ScatterDataCount = {
            'r1' : r1SC, 
            'r2' : r2SC
        }
        
        ScatterDataDensity = {
            'r1' : r1SD,
            'r2' : r2SD
        }
        
        
        #Plotting
        output_file(path+"_replicates.html")
          
        #Toolbar for each plot       
        TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save" 
        
        #Scatterplot
        scPlot = Scatter(ScatterDataCount, x='r1', y='r2', title="Counts vs Counts (N = 50.000) ", marker='circle', tools="resize, save", xlabel=headingR1+" [log10(counts+1)]", ylabel=headingR2+" [log10(counts+1)]") 
        tab1 = Panel(child=scPlot, title="Counts vs Counts")
        
        sdPlot = Scatter(ScatterDataDensity, x='r1', y='r2', title="Density vs Density (N = 50.000)", marker='circle', tools="resize, save", xlabel=headingR1+" [log10(counts+1)]", ylabel=headingR2+" [log10(counts+1)]") 
        tab2 = Panel(child=sdPlot, title="Density vs Density")
        
        tabs = Tabs(tabs=[ tab1, tab2 ])
        save(tabs)
       
       
    #===================================================================================
    #===================================================================================  
    def plot_junction(self):
        
        if self.input.endswith(".gz"):
            almnt_file = gzip.open(self.input, 'r') 
        else:        
            almnt_file = open(self.input, 'r')
         
        heading = self.output
        
        heading = heading.replace('\\','/')
        heading = heading.split(".html")
        path = heading[0]
        heading = heading[0].split("/")
        heading = heading[-1]
            
        ieJunc = {}
        eiJunc = {}
        eeJunc = {}
        iiJunc = {}
        
        eiJuncDist  = []
        eiJuncCount = []
        
        ieJuncDist  = []
        ieJuncCount = []
        
        eeJuncDist  = []
        eeJuncCount = []
        
        iiJuncDist  = []
        iiJuncCount = []
        
        ieJuncF = {}
        eiJuncF = {}
        eeJuncF = {}
        iiJuncF = {}
        
        eiJuncDistF  = []
        eiJuncCountF = []
        
        ieJuncDistF  = []
        ieJuncCountF = []
        
        eeJuncDistF  = []
        eeJuncCountF = []
        
        iiJuncDistF  = []
        iiJuncCountF = []
        
        types     = {}
        type      = []
        typeCount = []
        
        regions     = {}
        region      = []
        regionCount = []
        
        bv = None
        
        chromNorm   = {}
        featureNorm = {}
               
        #Adding information to the arrays for the plots
        for line in almnt_file:
            
            if line.startswith("track"):
                line = line.split("\n")
                line = line[0].split(" ")
                
                if line[1] == "chr":
                    chrNorm[line[2]] = int(line[3])
                elif line[1] == "type":
                    featureNorm[line[2]] = int(line[3])
                else:
                    error = "Unknown track annotation: "+line[1]+". Check your data!!"
                    raise ValueError(error)
                
            else:
                line = line.split("\n")
                line = line[0].split("\t")
                
                if line[9] == "exon":
                    if not int(line[4]) > 4000: 
                        
                        #Unfiltered
                        if not ieJunc.has_key(int(line[4])):
                            ieJunc[int(line[4])] = 1
                        else:
                            ieJunc[int(line[4])] += 1  
                            
                        if not eeJunc.has_key(int(line[4])):
                            eeJunc[int(line[4])] = 1
                        else:
                            eeJunc[int(line[4])] += 1  
                            
                        #Filtered
                        if int(line[6]) == 3 or int(line[6]) == 2:
                            if not ieJuncF.has_key(int(line[4])):
                                ieJuncF[int(line[4])] = 1
                            else:
                                ieJuncF[int(line[4])] += 1  
                            
                            if not eeJuncF.has_key(int(line[4])):
                                eeJuncF[int(line[4])] = 1
                            else:
                                eeJuncF[int(line[4])] += 1   
                                                   
                    if not int(line[5]) < -4000:
                        
                        #Unfiltered
                        if not eiJunc.has_key(int(line[5])):
                            eiJunc[int(line[5])] = 1
                        else:
                            eiJunc[int(line[5])] += 1
                            
                        if not eeJunc.has_key(int(line[5])):
                            eeJunc[int(line[5])] = 1
                        else:
                            eeJunc[int(line[5])] += 1 
                        
                        #Filtered
                        if int(line[6]) == 3 or int(line[6]) == 1:
                            if not eiJuncF.has_key(int(line[5])):
                                eiJuncF[int(line[5])] = 1
                            else:
                                eiJuncF[int(line[5])] += 1
                            
                            if not eeJuncF.has_key(int(line[5])):
                                eeJuncF[int(line[5])] = 1
                            else:
                                eeJuncF[int(line[5])] += 1 
              
                    if line[8] == "protein_coding":
                        line[8] = line[8]+"_exon"
                        
                    if not regions.has_key("exon"):
                        regions["exon"] = 1
                    else:
                        regions["exon"] += 1
                    
                elif line[9] == "intron":
                    if not int(line[4]) > 4000: 
                        
                        #Unfiltered
                        if not eiJunc.has_key(int(line[4])):
                            eiJunc[int(line[4])] = 1
                        else:
                            eiJunc[int(line[4])] += 1    
                            
                        if not iiJunc.has_key(int(line[4])):
                            iiJunc[int(line[4])] = 1
                        else:
                            iiJunc[int(line[4])] += 1 
                            
                        #Filtered
                        if int(line[6]) == 3 or int(line[6]) == 2:
                            if not eiJuncF.has_key(int(line[4])):
                                eiJuncF[int(line[4])] = 1
                            else:
                                eiJuncF[int(line[4])] += 1  
                            
                            if not iiJuncF.has_key(int(line[4])):
                                iiJuncF[int(line[4])] = 1
                            else:
                                iiJuncF[int(line[4])] += 1
                        
                                         
                    if not int(line[5]) < -4000:
                        
                        #Unfiltered
                        if not ieJunc.has_key(int(line[5])):
                            ieJunc[int(line[5])] = 1
                        else:
                            ieJunc[int(line[5])] += 1
                            
                        if not iiJunc.has_key(int(line[5])):
                            iiJunc[int(line[5])] = 1
                        else:
                            iiJunc[int(line[5])] += 1  
                            
                        #Filtered
                        if int(line[6]) == 3 or int(line[6]) == 1:
                            if not ieJuncF.has_key(int(line[5])):
                                ieJuncF[int(line[5])] = 1
                            else:
                                ieJuncF[int(line[5])] += 1
                            
                            if not iiJuncF.has_key(int(line[5])):
                                iiJuncF[int(line[5])] = 1
                            else:
                                iiJuncF[int(line[5])] += 1  
                              
                    if line[8] == "protein_coding":
                        line[8] = line[8]+"_intron"
                        
                    if not regions.has_key("intron"):
                        regions["intron"] = 1
                    else:
                        regions["intron"] += 1
                else:
                    if not regions.has_key("intergenic"):
                        regions["intergenic"] = 1
                    else:
                        regions["intergenic"] += 1
                
                if not types.has_key(line[8]):
                    types[line[8]] = 1
                else:
                    types[line[8]] += 1
        
        #Filling up dicts with 0 counts if for some distances no values are available        
        for i in range(-4000, 4000):
            if not eiJunc.has_key(i):
                eiJunc[i] = 0
            if not eiJuncF.has_key(i):
                eiJuncF[i] = 0   
            if not ieJunc.has_key(i):
                ieJunc[i] = 0
            if not ieJuncF.has_key(i):
                ieJuncF[i] = 0   
            if not eeJunc.has_key(i):
                eeJunc[i] = 0
            if not eeJuncF.has_key(i):
                eeJuncF[i] = 0   
            if not iiJunc.has_key(i):
                iiJunc[i] = 0
            if not iiJuncF.has_key(i):
                iiJuncF[i] = 0     
        
        #Sort for plotting against   
        eiJunc = OrderedDict(sorted(eiJunc.items(), key=lambda x: x[0]))
        ieJunc = OrderedDict(sorted(ieJunc.items(), key=lambda x: x[0]))
        eeJunc = OrderedDict(sorted(eeJunc.items(), key=lambda x: x[0]))
        iiJunc = OrderedDict(sorted(iiJunc.items(), key=lambda x: x[0]))
        
        eiJuncF = OrderedDict(sorted(eiJuncF.items(), key=lambda x: x[0]))
        ieJuncF = OrderedDict(sorted(ieJuncF.items(), key=lambda x: x[0]))
        eeJuncF = OrderedDict(sorted(eeJuncF.items(), key=lambda x: x[0]))
        iiJuncF = OrderedDict(sorted(iiJuncF.items(), key=lambda x: x[0]))
        
        regions = OrderedDict(sorted(regions.items(), key=lambda x: (-x[1], x[0])))
        types   = OrderedDict(sorted(types.items(), key=lambda x: (-x[1], x[0])))
            
        #Preparing data for plotting
        for k in eiJunc:
            eiJuncDist.append(k)
            eiJuncCount.append(eiJunc[k])
         
        for k in ieJunc:
            ieJuncDist.append(k)
            ieJuncCount.append(ieJunc[k])
            
        for k in eeJunc:
            eeJuncDist.append(k)
            eeJuncCount.append(eeJunc[k])
         
        for k in iiJunc:
            iiJuncDist.append(k)
            iiJuncCount.append(iiJunc[k])
            
        for k in eiJuncF:
            eiJuncDistF.append(k)
            eiJuncCountF.append(eiJuncF[k])
         
        for k in ieJuncF:
            ieJuncDistF.append(k)
            ieJuncCountF.append(ieJuncF[k])
            
        for k in eeJuncF:
            eeJuncDistF.append(k)
            eeJuncCountF.append(eeJuncF[k])
         
        for k in iiJuncF:
            iiJuncDistF.append(k)
            iiJuncCountF.append(iiJuncF[k])
            
        for k in regions:
            region.append(k)
            regionCount.append(regions[k])
            
        for k in types:
            type.append(k)
            typeCount.append(types[k])
            
        #Data dicts              
        ieData = {}  
        eiData = {}
        eeData = {}
        iiData = {}
        
        ieDataF = {}  
        eiDataF = {}
        eeDataF = {}
        iiDataF = {}
        
        typeData   = {}
        regionData = {}
        
        ieData = {
            'ieJuncDist' : ieJuncDist,
            'ieJuncCounts' : ieJuncCount        
        }
        
        eiData = {
            'eiJuncDist' : eiJuncDist,
            'eiJuncCounts' : eiJuncCount
        }
        
        eeData = {
            'eeJuncDist' : eeJuncDist,
            'eeJuncCounts' : eeJuncCount
        }
        
        iiData = {
            'iiJuncDist' : iiJuncDist,
            'iiJuncCounts' : iiJuncCount
        }
        
        ieDataF = {
            'ieJuncDist' : ieJuncDistF,
            'ieJuncCounts' : ieJuncCountF        
        }
        
        eiDataF = {
            'eiJuncDist' : eiJuncDistF,
            'eiJuncCounts' : eiJuncCountF
        }
        
        eeDataF = {
            'eeJuncDist' : eeJuncDistF,
            'eeJuncCounts' : eeJuncCountF
        }
        
        iiDataF = {
            'iiJuncDist' : iiJuncDistF,
            'iiJuncCounts' : iiJuncCountF
        }
        
        typeData = {
            'type' : type,
            'typeCounts' : typeCount       
        }
        
        regionData = {
            'region' : region,
            'regionCounts' : regionCount            
        }
             
                  
        #Toolbar for each plot       
        TOOLS="resize,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,save"
         
         
        #Junction
        output_file(path+".html")
        
        ##########################################################################################################
        #Unfiltered
        #Intron-Exon Junction
        #---------------------------------------------------------------------------------------------------------
        ieP = figure(title=heading+"_ieJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        ieP.line(ieJuncDist, ieJuncCount)
        
        ieData_table = DataTable(source=ColumnDataSource(ieData), columns=[TableColumn(field="ieJuncDist", title="Distances"), TableColumn(field="ieJuncCounts", title="Count")], width=600, height=600)
        
        iePlot = vplot(
                hplot(ieP, ieData_table)
        )
        tab1 = Panel(child=iePlot, title="Intron-Exon Junction")
        
        ieDataOutput = open(path+"_ieData.txt", "w")
        ieDataOutput.write("Distances\tCount\n")
        for i in range(0, len(ieJuncDist)):
            ieDataOutput.write(str(ieJuncDist[i]) + "\t" + str(ieJuncCount[i]) + "\n")
        
        ieDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        
        #Exon-Intron Junction
        #---------------------------------------------------------------------------------------------------------
        eiP = figure(title=heading+"_eiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        eiP.line(eiJuncDist, eiJuncCount)
        
        eiData_table = DataTable(source=ColumnDataSource(eiData), columns=[TableColumn(field="eiJuncDist", title="Distances"), TableColumn(field="eiJuncCounts", title="Count")], width=600, height=600)
        
        eiPlot = vplot(
                hplot(eiP, eiData_table)
        )
        tab2 = Panel(child=eiPlot, title="Exon-Intron Junction")
        
        eiDataOutput = open(path+"_eiData.txt", "w")
        eiDataOutput.write("Distances\tCount\n")
        for i in range(0, len(eiJuncDist)):
            eiDataOutput.write(str(eiJuncDist[i]) + "\t" + str(eiJuncCount[i]) + "\n")
        
        eiDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        
        #Exon-Exon Junction
        #---------------------------------------------------------------------------------------------------------
        eeP = figure(title=heading+"_eeJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        eeP.line(eeJuncDist, eeJuncCount)
        
        eeData_table = DataTable(source=ColumnDataSource(eeData), columns=[TableColumn(field="eeJuncDist", title="Distances"), TableColumn(field="eeJuncCounts", title="Count")], width=600, height=600)
        
        eePlot = vplot(
                hplot(eeP, eeData_table)
        )
        tab3 = Panel(child=eePlot, title="Exon-Exon Junction")
        
        eeDataOutput = open(path+"_eeData.txt", "w")
        eeDataOutput.write("Distances\tCount\n")
        for i in range(0, len(eeJuncDist)):
            eeDataOutput.write(str(eeJuncDist[i]) + "\t" + str(eeJuncCount[i]) + "\n")
        
        eeDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        
        #Intron-Intron Junction
        #---------------------------------------------------------------------------------------------------------
        iiP = figure(title=heading+"_iiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        iiP.line(iiJuncDist, iiJuncCount)
        
        iiData_table = DataTable(source=ColumnDataSource(iiData), columns=[TableColumn(field="iiJuncDist", title="Distances"), TableColumn(field="iiJuncCounts", title="Count")], width=600, height=600)
        
        iiPlot = vplot(
                hplot(iiP, iiData_table)
        )
        tab4 = Panel(child=iiPlot, title="Intron-Intron Junction")
        
        iiDataOutput = open(path+"_iiData.txt", "w")
        iiDataOutput.write("Distances\tCount\n")
        for i in range(0, len(iiJuncDist)):
            iiDataOutput.write(str(iiJuncDist[i]) + "\t" + str(iiJuncCount[i]) + "\n")
        
        iiDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        ##########################################################################################################
        
        ##########################################################################################################
        #Filtered
        #Intron-Exon Junction
        #---------------------------------------------------------------------------------------------------------
        iePF = figure(title=heading+"_ieJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        iePF.line(ieJuncDistF, ieJuncCountF)
        
        ieData_tableF = DataTable(source=ColumnDataSource(ieDataF), columns=[TableColumn(field="ieJuncDist", title="Distances"), TableColumn(field="ieJuncCounts", title="Count")], width=600, height=600)
        
        iePlotF = vplot(
                hplot(iePF, ieData_tableF)
        )
        tab5 = Panel(child=iePlotF, title="Intron-Exon Junction")
        
        ieFDataOutput = open(path+"_ieFData.txt", "w")
        ieFDataOutput.write("Distances\tCount\n")
        for i in range(0, len(ieJuncDistF)):
            ieFDataOutput.write(str(ieJuncDistF[i]) + "\t" + str(ieJuncCountF[i]) + "\n")
        
        ieFDataOutput.close()
            
        #---------------------------------------------------------------------------------------------------------
        
        #Exon-Intron Junction
        #---------------------------------------------------------------------------------------------------------
        eiPF = figure(title=heading+"_eiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        eiPF.line(eiJuncDistF, eiJuncCountF)
        
        eiData_tableF = DataTable(source=ColumnDataSource(eiDataF), columns=[TableColumn(field="eiJuncDist", title="Distances"), TableColumn(field="eiJuncCounts", title="Count")], width=600, height=600)
        
        eiPlotF = vplot(
                hplot(eiPF, eiData_tableF)
        )
        tab6 = Panel(child=eiPlotF, title="Exon-Intron Junction")
        
        eiFDataOutput = open(path+"_eiFData.txt", "w")
        eiFDataOutput.write("Distances\tCount\n")
        for i in range(0, len(eiJuncDistF)):
            eiFDataOutput.write(str(eiJuncDistF[i]) + "\t" + str(eiJuncCountF[i]) + "\n")
        
        eiFDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        
        #Exon-Exon Junction
        #---------------------------------------------------------------------------------------------------------
        eePF = figure(title=heading+"_eeJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        eePF.line(eeJuncDistF, eeJuncCountF)
        
        eeData_tableF = DataTable(source=ColumnDataSource(eeDataF), columns=[TableColumn(field="eeJuncDist", title="Distances"), TableColumn(field="eeJuncCounts", title="Count")], width=600, height=600)
        
        eePlotF = vplot(
                hplot(eePF, eeData_tableF)
        )
        tab7 = Panel(child=eePlotF, title="Exon-Exon Junction")
        
        eeFDataOutput = open(path+"_eeFData.txt", "w")
        eeFDataOutput.write("Distances\tCount\n")
        for i in range(0, len(eeJuncDistF)):
            eeFDataOutput.write(str(eeJuncDistF[i]) + "\t" + str(eeJuncCountF[i]) + "\n")
        
        eeFDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        
        #Intron-Intron Junction
        #---------------------------------------------------------------------------------------------------------
        iiPF = figure(title=heading+"_iiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Counts')
        iiPF.line(iiJuncDistF, iiJuncCountF)
        
        iiData_tableF = DataTable(source=ColumnDataSource(iiDataF), columns=[TableColumn(field="iiJuncDist", title="Distances"), TableColumn(field="iiJuncCounts", title="Count")], width=600, height=600)
        
        iiPlotF = vplot(
                hplot(iiPF, iiData_tableF)
        )
        tab8 = Panel(child=iiPlotF, title="Intron-Intron Junction")
        
        iiFDataOutput = open(path+"_iiFData.txt", "w")
        iiFDataOutput.write("Distances\tCount\n")
        for i in range(0, len(iiJuncDistF)):
            iiFDataOutput.write(str(iiJuncDistF[i]) + "\t" + str(iiJuncCountF[i]) + "\n")
        
        iiFDataOutput.close()
        #---------------------------------------------------------------------------------------------------------
        ##########################################################################################################    
        
        #Types and Regions
        #---------------------------------------------------------------------------------------------------------
#         typeP = Bar(typeData, title=heading+"_types", values='typeCounts', label=CatAttr(columns=['type'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Counts")
#         typeData_table = DataTable(source=ColumnDataSource(typeData), columns=[TableColumn(field="type", title="Types"),TableColumn(field="typeCounts", title="Counts")], width=600, height=600)
#          
#         tPlot = vplot(
#                 hplot(typeP, typeData_table),
#         )
#         
#         typeTab = Panel(child=tPlot, title="Types")  
#         
#         regionP = Bar(regionData, title=heading+"_regions", values='regionCounts', label=CatAttr(columns=['region'], sort=False),tools=TOOLS, xlabel="Regions", ylabel="Counts")
#         regionData_table = DataTable(source=ColumnDataSource(regionData), columns=[TableColumn(field="region", title="Regions"),TableColumn(field="regionCounts", title="Counts")], width=600, height=600)
#         
#         rPlot = vplot(
#                 hplot(regionP, regionData_table)
#         )
#         
#         regionTab = Panel(child=rPlot, title="Regions")
        
            
        #---------------------------------------------------------------------------------------------------------
        
        unfTabs = Tabs(tabs=[tab1, tab2, tab3, tab4])
        unfiltered = Panel(child=unfTabs, title="Junction unfiltered")
        
        fTabs = Tabs(tabs=[tab5, tab6, tab7, tab8])
        filtered = Panel(child=fTabs, title="Junction filtered")
        jt = Tabs(tabs=[unfiltered, filtered])
         
        save(jt)

        
        #===================================================================================