# --------------------------------------------------
# bokehCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
# Institution: EMBL Heidelberg
# Date: December 2015
# --------------------------------------------------

import gzip
from math import log10
from collections import OrderedDict
from random import randint
from numpy import log2

from bokeh.plotting import figure
from bokeh.charts import Scatter, Histogram, output_file, save, vplot, hplot, Bar
from bokeh.models.widgets import DataTable, TableColumn, Panel, Tabs
from bokeh.models import ColumnDataSource
from bokeh.charts.attributes import CatAttr 
    
class bokehCLIP:
    
    fInput = ''
    fCompare = ''
    fOutput = ''
    choice = ''
    
    def __init__(self, options):
        
        if hasattr(options, 'input'):
            self.fInput = options.input[0] 
            
        if hasattr(options, 'compare'):
            self.fCompare = options.compare
            
        if hasattr(options, 'output'):
            self.fOutput = options.output[0]
            
        if hasattr(options, 'choice'):
            self.choice = options.choice
        
    #===================================================================================
    #===================================================================================
    def plot_read(self):
        
        if self.fInput.endswith(".gz"):
            almnt_file = gzip.open(self.fInput, 'r') 
        else:        
            almnt_file = open(self.fInput, 'r')
               
        heading = self.fOutput
        
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
        
        
        size = len(almnt_file.read())
        almnt_file.seek(0)
        
        if size == 0:
            TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save"
            output_file(path+".html")
            hP = Histogram([-0.1,0.1], title="No values for this type of cross-link site available!", xlabel="No values", ylabel="No values", tools=TOOLS)
            save(hP)
        else:
               
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
            
            rlDataOutput = open(path+"_read-length-data.txt", "w")
            rlDataOutput.write("Read Length\tNumber of Reads\n")
            for i in range(0, len(rlKeys)):
                rlDataOutput.write(str(rlKeys[i]) + "\t" + str(rlCounts[i]) + "\n")
            
            rlDataOutput.close()
             
            #Duplicates         
            dpData_table = DataTable(source=ColumnDataSource(dpData), columns=[TableColumn(field="dpKeys", title="Duplication ratio"),TableColumn(field="dpCounts", title="Number of reads")], width=600, height=600)
            
            tab2 = Panel(child=dpData_table, title="Duplicates")
            
            dpDataOutput = open(path+"_duplicates-data.txt", "w")
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
            
            chromDataOutput = open(path+"_chromosomes-data.txt", "w")
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
            
            strandDataOutput = open(path+"_strands-data.txt", "w")
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
        
        if self.fInput.endswith(".gz"):
            almnt_file = gzip.open(self.fInput, 'r') 
        else:        
            almnt_file = open(self.fInput, 'r')
        
        heading = self.fOutput
        
        heading = heading.replace('\\','/')
        heading = heading.split(".html")
        path    = heading[0]
        heading = heading[0].split("/")
        heading = heading[-1]
        
        empty = True
        
        for line in almnt_file:
            if not line.startswith("#track"):
                empty = False
        almnt_file.seek(0)
        
        if empty == True:
            TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save"
            output_file(path+".html")
            hP = Histogram([-0.1,0.1], title="No values for this type of cross-link site available!", xlabel="No values", ylabel="No values", tools=TOOLS)
            save(hP)
        else:
                  
            chrNorm     = {}
            featureNorm = {}
            
            counts      = {}
            countKeys   = []
            countCounts = []
            total_counts = 0.0
    
            
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
            vType      = []
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
        
            totalField = 0 # Number of fields in the input file
            for line in almnt_file:
                
                if line.startswith("#track"):
                    l = line.split("\n")
                    l = l[0].split(" ")

                    if len(l)==2:
                        totalField = 18
                    elif len(l)==1:
                        totalField = 17

                    elif l[1] == "chr":
                        chrNorm[l[2]] = int(l[3])
                    elif l[1] == "type":
                        featureNorm[l[2]] = int(l[3])
                    else:
                        error = "Unknown track annotation: "+l[1]+". Check your data!!"
                        raise ValueError(error)

                elif line.startswith('Chromosome') or line  == '\n':
                    continue
                elif "intergenic" in line:
                    line = line.split("\n")
                    line = line[0].split("\t")
                    
                    if not countsPerType.has_key("intergenic"):
                        countsPerType["intergenic"] = int(line[totalField-6])
                    else:
                        countsPerType["intergenic"] += int(line[totalField-6])



                else:

                    line = line.split("\n")
                    line = line[0].split("\t")

                    if not chromosomes.has_key(line[0]):
                        chromosomes[line[0]] = 1
                    else:
                        chromosomes[line[0]] += 1
                        
                    if not counts.has_key(int(line[totalField-6])):
                        counts[int(line[totalField-6])] = 1
                        total_counts += int(line[totalField-6])
                    else:
                        counts[int(line[totalField-6])] += 1
                        total_counts += int(line[totalField-6])

                    if not density.has_key(float(line[totalField-3])):
                        density[float(line[totalField-6])] = 1
                    else:
                        density[float(line[totalField-6])] += 1
                        
                    if not duplicates.has_key(int(line[totalField-2])):
                        duplicates[int(line[totalField-2])] = 1
                    else:
                        duplicates[int(line[totalField-2])] += 1

                    if not countsPergene.has_key(line[3]):
                        countsPergene[line[3]] = int(line[totalField-6])
                    else:
                        countsPergene[line[3]] += int(line[totalField-6])
                    if line[totalField-11] == "exon" and line[totalField-8] == "protein_coding":
                        line[totalField-8] = "protein_coding_exon"
                        if not types.has_key(line[totalField-8]):
                            types[line[totalField-8]] = 1
                        else:
                            types[line[totalField-8]] += 1
                            
                        if not dcQC.has_key(line[totalField-8]):
                            dcQC[line[totalField-8]] = [int(line[totalField-2]), int(line[totalField-6])]
                        else:
                            dcQC[line[totalField-8]][0] += int(line[totalField-2])
                            dcQC[line[totalField-8]][1] += int(line[totalField-6])
                            
                            
                    elif line[totalField-11] == "intron" and line[totalField-8] == "protein_coding":
                        line[totalField-8] = "protein_coding_intron"
                        if not types.has_key(line[totalField-8]):
                            types[line[totalField-8]] = 1
                        else:
                            types[line[totalField-8]] += 1
                            
                        if not dcQC.has_key(line[totalField-8]):
                            dcQC[line[totalField-8]] = [int(line[totalField-2]), int(line[totalField-6])]
                        else:
                            dcQC[line[totalField-8]][0] += int(line[totalField-2])
                            dcQC[line[totalField-8]][1] += int(line[totalField-6])

                    elif line[totalField-11] == "5UTR" and line[totalField-8] == "protein_coding":
                        line[totalField-8] = "protein_coding_5UTR"
                        #line[10] = "5UTR"

                        if not types.has_key(line[totalField-8]):
                            types[line[totalField-8]] = 1
                        else:
                            types[line[totalField-8]] += 1

                        if not dcQC.has_key(line[totalField-8]):
                            dcQC[line[totalField-8]] = [int(line[totalField-2]), int(line[totalField-6])]
                        else:
                            dcQC[line[totalField-8]][0] += int(line[totalField-2])
                            dcQC[line[totalField-8]][1] += int(line[totalField-6])

                    elif line[totalField-11] == "3UTR" and line[totalField-8] == "protein_coding":
                        line[totalField-8] = "protein_coding_3UTR"
                        #line[10] = "3UTR"
                        if not types.has_key(line[totalField-8]):
                            types[line[totalField-8]] = 1
                        else:
                            types[line[totalField-8]] += 1

                        if not dcQC.has_key(line[totalField-8]):
                            dcQC[line[totalField-8]] = [int(line[totalField-2]), int(line[totalField-6])]
                        else:
                            dcQC[line[totalField-8]][0] += int(line[totalField-2])
                            dcQC[line[totalField-8]][1] += int(line[totalField-6])

                    elif line[totalField-11] == "CDS" and line[totalField-8] == "protein_coding":
                        line[totalField-8] = "protein_coding_CDS"
                        #line[10] = "CDS"
                        if not types.has_key(line[totalField-8]):
                            types[line[totalField-8]] = 1
                        else:
                            types[line[totalField-8]] += 1

                        if not dcQC.has_key(line[totalField-8]):
                            dcQC[line[totalField-8]] = [int(line[totalField-2]), int(line[totalField-6])]
                        else:
                            dcQC[line[totalField-8]][0] += int(line[totalField-2])
                            dcQC[line[totalField-8]][1] += int(line[totalField-6])

                    else:
                        if not types.has_key(line[totalField-8]):
                            types[line[totalField-8]] = 1
                        else:
                            types[line[totalField-8]] += 1
                            
                        if not dcQC.has_key(line[totalField-8]):
                            dcQC[line[totalField-8]] = [int(line[totalField-2]), int(line[totalField-6])]
                        else:
                            dcQC[line[totalField-8]][0] += int(line[totalField-2])
                            dcQC[line[totalField-8]][1] += int(line[totalField-6])
                            
                    if not regions.has_key(line[totalField-11]):
                        regions[line[totalField-11]] = 1
                    else:
                        regions[line[totalField-11]] += 1
                        
                    if not strands.has_key(line[totalField-12]):
                        strands[line[totalField-12]] = 1
                    else:
                        strands[line[totalField-12]] += 1
                    
                    
            for k in dcQC:
                dupPerType[k] = dcQC[k][0]
                countsPerType[k] = dcQC[k][1]
                dcQC[k] = dcQC[k][0]/float(dcQC[k][1])
                
            chrN = {}
            cn = {}
            dn = {}
            total_length = 0.0
            for k in chromosomes:
                chrN[k] = (chromosomes[k]/ float(chrNorm[k]))

            for k in countsPerType:
                if not k == "intergenic":
                    cn[k] = (countsPerType[k] / float(total_counts))
                
            for k in dupPerType:
                if not k == "intergenic":
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
                vType.append(k)
                typeCount.append(types[k])
            
            typeData = {
                'type' : vType,
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
            
            cptDataOutput = open(path+"_counts-per-type-data.txt", "w")
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
            
            cptNDataOutput = open(path+"_counts-per-typeNorm.txt", "w")
            cptNDataOutput.write("Types\tNormalized count\n")
            for i in range(0, len(cptNKeys)):
                cptNDataOutput.write(str(cptNKeys[i]) + "\t" + str(cptNCount[i]) + "\n")
            
            cptNDataOutput.close()
            
            dptP = Bar(dptData, title=heading+"_Duplicates_per_Type", values='dptCount', label=CatAttr(columns=['dptKeys'], sort=False),tools=TOOLS, xlabel="Types", ylabel="Count")
            dptData_table = DataTable(source=ColumnDataSource(dptData), columns=[TableColumn(field="dptKeys", title="Types"),TableColumn(field="dptCount", title="Number of types")], width=600, height=600)
            
            dptPlot = vplot(
                    hplot(dptP, dptData_table)
            )
            
            tab3 = Panel(child=dptPlot, title="Duplicates per type")
            
            dptDataOutput = open(path+"_duplicates-per-type.txt", "w")
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
            
            dptNDataOutput = open(path+"_duplicates-per-typeNorm.txt", "w")
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
            
            dcQCDataOutput = open(path+"_dupicates-to-count-ratio.txt", "w")
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
            
            typeDataOutput = open(path+"_types-per-feature.txt", "w")
            typeDataOutput.write("Types\tNumber of types\n")
            for i in range(0, len(vType)):
                typeDataOutput.write(str(vType[i]) + "\t" + str(typeCount[i]) + "\n")
            
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
            
            chromDataOutput = open(path+"_chromosomes.txt", "w")
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
            
            chromDataOutput = open(path+"_chromosomes-norm.txt", "w")
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
            
            strandDataOutput = open(path+"_strands.txt", "w")
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
            
            regionDataOutput = open(path+"_total-regions.txt", "w")
            regionDataOutput.write("Regions\tNumber of Regions\n")
            for i in range(0, len(region)):
                regionDataOutput.write(str(region[i]) + "\t" + str(regionCount[i]) + "\n")
            
            regionDataOutput.close()
            #--------------------------------------------------
            
            #Counts
            cData_table = DataTable(source=ColumnDataSource(countData), columns=[TableColumn(field="cKeys", title="Counts"),TableColumn(field="cCounts", title="Number of types")], width=600, height=600)
             
            tab10 = Panel(child=cData_table, title="Counts")
            
            countDataOutput = open(path+"_total-counts.txt", "w")
            countDataOutput.write("Counts\tNumber of counts\n")
            for i in range(0, len(countKeys)):
                countDataOutput.write(str(countKeys[i]) + "\t" + str(countCounts[i]) + "\n")
            
            countDataOutput.close()
               
            #Duplicates         
            dpData_table = DataTable(source=ColumnDataSource(dpData), columns=[TableColumn(field="dpKeys", title="Duplication ratio"),TableColumn(field="dpCounts", title="Number of types")], width=600, height=600)
             
            tab11 = Panel(child=dpData_table, title="Duplicates")
            
            dupDataOutput = open(path+"_total-duplicates.txt", "w")
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
        
        if self.fInput.endswith(".gz"):
            replicate1 = gzip.open(self.fInput, 'r') 
        else:        
            replicate1 = open(self.fInput, 'r')
            
        if self.fCompare.endswith(".gz"):
            replicate2 = gzip.open(self.fCompare, 'r') 
        else:        
            replicate2 = open(self.fCompare, 'r')
            
        headingR1 = self.fInput
        headingR1 = headingR1.replace('\\','/')
        headingR1 = headingR1.split('/')
        headingR1 = headingR1[-1]
             
        headingR2 = self.fCompare
        headingR2 = headingR2.replace('\\','/')
        headingR2 = headingR2.split('/')
        headingR2 = headingR2[-1]
         
        heading = self.fOutput
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
            
            sd = [randint(0,len(r1ScatterCount)-1) for p in range(0,50000)] # @UnusedVariable
        
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
        TOOLS="resize,save" 
        
        #Scatterplot
        scPlot = Scatter(ScatterDataCount, x='r1', y='r2', title="Counts vs Counts (N = 50.000) ", marker='circle', tools=TOOLS, xlabel=headingR1+" [log10(counts+1)]", ylabel=headingR2+" [log10(counts+1)]") 
        tab1 = Panel(child=scPlot, title="Counts vs Counts")
        
        sdPlot = Scatter(ScatterDataDensity, x='r1', y='r2', title="Density vs Density (N = 50.000)", marker='circle', tools=TOOLS, xlabel=headingR1+" [log10(counts+1)]", ylabel=headingR2+" [log10(counts+1)]") 
        tab2 = Panel(child=sdPlot, title="Density vs Density")
        
        tabs = Tabs(tabs=[ tab1, tab2 ])
        save(tabs)
       
       
    #===================================================================================
    #===================================================================================  
    def plot_junction(self):
        
        if self.fInput.endswith(".gz"):
            almnt_file = gzip.open(self.fInput, 'r') 
        else:        
            almnt_file = open(self.fInput, 'r')
         
        heading = self.fOutput
        
        heading = heading.replace('\\','/')
        heading = heading.split(".html")
        path = heading[0]
        heading = heading[0].split("/")
        heading = heading[-1]
        
        empty = True
        
        for line in almnt_file:
            if not line.startswith("track"):
                empty = False
        almnt_file.seek(0)
        
        if empty == True:
            TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save"
            output_file(path+".html")
            hP = Histogram([-0.1,0.1], title="No values for this type of cross-link site available!", xlabel="No values", ylabel="No values", tools=TOOLS)
            save(hP)
        else:
            
            ieJunc = {}
            eiJunc = {}
            eeJunc = {}
            iiJunc = {}
            igeJunc = {}
            eigJunc = {}
            es11Junc = {}
            ee11Junc = {}
            
            eiJuncDist  = []
            eiJuncCount = []
            
            ieJuncDist  = []
            ieJuncCount = []
            
            eeJuncDist  = []
            eeJuncCount = []
            
            iiJuncDist  = []
            iiJuncCount = []
            
            igeJuncDist = []
            igeJuncCount = []
            
            eigJuncDist = []
            eigJuncCount = []
            
            ee11JuncDist = []
            ee11JuncCount = []
            
            es11JuncDist = []
            es11JuncCount = []
            
            
            ieJuncF = {}
            eiJuncF = {}
            eeJuncF = {}
            iiJuncF = {}
            igeJuncF = {}
            eigJuncF = {}
            es11JuncF = {}
            ee11JuncF = {}
            
            eiJuncDistF  = []
            eiJuncCountF = []
            
            ieJuncDistF  = []
            ieJuncCountF = []
            
            eeJuncDistF  = []
            eeJuncCountF = []
            
            iiJuncDistF  = []
            iiJuncCountF = []
            
            igeJuncDistF = []
            igeJuncCountF = []
            
            eigJuncDistF = []
            eigJuncCountF = []
            
            ee11JuncDistF = []
            ee11JuncCountF = []
            
            es11JuncDistF = []
            es11JuncCountF = []
            
            types     = {}
            vType      = []
            typeCount = []
            
            regions     = {}
            region      = []
            regionCount = []
            
            #Normalization
            dtl = {}
             
            for line in almnt_file:
                line = line.split("\n")
                line = line[0].split("\t") 
                 
                if not line[9] == "intergenic":
                    totalLength = abs(int(line[4])) + abs(int(line[5]))
                    expr = line[1]+line[2]+line[3]
                    if not dtl.has_key(line[9]):
                        dtl[line[9]] = {totalLength : [expr]}
                    else:
                        if not dtl[line[9]].has_key(totalLength):
                            dtl[line[9]][totalLength] = [expr]
                        else:
                            dtl[line[9]][totalLength].append(expr)
            
            
            #Setting iterator back to 0
            almnt_file.seek(0)
            
            efl = {}
            ifl = {}
            
            for f in dtl:
                for l in dtl[f]:
                    unique_words = set(dtl[f][l])
                    unique_word_count = len(unique_words)
                    if f == "exon":
                        efl[l] = unique_word_count
                    else:
                        ifl[l] = unique_word_count

            #Adding information to the arrays for the plots
            for line in almnt_file:
                 
                line = line.split("\n")
                line = line[0].split("\t")      
                        
                if line[9] == "exon":
                     
                    totalLength = abs(int(line[4])) + abs(int(line[5]))
                            
                    dts = int(line[4])
                    dte = int(line[5])
                    
                    line[10] = line[10].split("/")
                                       
                    if not dts > 200: 
                        
                        if line[10][0] == '1' and line[10][1] == '1':
                            
                            #Unfiltered
                            if not es11Junc.has_key(dts):
                                es11Junc[dts] = 1/float(efl[totalLength])
                            else:
                                es11Junc[dts] += 1/float(efl[totalLength])
                                        
                            #Filtered
                            if int(line[6]) == 3 or int(line[6]) == 2:
                                if not es11JuncF.has_key(dts):
                                    es11JuncF[dts] = 1/float(efl[totalLength])
                                else:
                                    es11JuncF[dts] += 1/float(efl[totalLength])
                                     
                        elif (line[10][0] == '1' and not line[10][1] == '1') or (line[10][0] == line[10][1]):
                            
                            #Unfiltered
                            if not igeJunc.has_key(dts):
                                igeJunc[dts] = 1/float(efl[totalLength])
                            else:
                                igeJunc[dts] += 1/float(efl[totalLength])
                                        
                            #Filtered
                            if int(line[6]) == 3 or int(line[6]) == 2:
                                if not igeJuncF.has_key(dts):
                                    igeJuncF[dts] = 1/float(efl[totalLength])
                                else:
                                    igeJuncF[dts] += 1/float(efl[totalLength])
                                
                        else:
                                        
                            #Unfiltered
                            if not ieJunc.has_key(dts):
                                ieJunc[dts] = 1/float(efl[totalLength])
                            else:
                                ieJunc[dts] += 1/float(efl[totalLength])
                                        
                            if not eeJunc.has_key(dts):
                                eeJunc[dts] = 1/float(efl[totalLength])
                            else:
                                eeJunc[dts] += 1/float(efl[totalLength])
                                        
                            #Filtered
                            if int(line[6]) == 3 or int(line[6]) == 2:
                                if not ieJuncF.has_key(dts):
                                    ieJuncF[dts] = 1/float(efl[totalLength])
                                else:
                                    ieJuncF[dts] += 1/float(efl[totalLength])
                                        
                                if not eeJuncF.has_key(dts):
                                    eeJuncF[dts] = 1/float(efl[totalLength])
                                else:
                                    eeJuncF[dts] += 1/float(efl[totalLength])
                                                           
                    if not dte < -200:
                        
                        if line[10][0] == '1' and line[10][1] == '1':
                            
                            #Unfiltered
                            if not ee11Junc.has_key(dte):
                                ee11Junc[dte] = 1/float(efl[totalLength])
                            else:
                                ee11Junc[dte] += 1/float(efl[totalLength])
                                        
                            #Filtered
                            if int(line[6]) == 3 or int(line[6]) == 2:
                                if not ee11JuncF.has_key(dte):
                                    ee11JuncF[dte] = 1/float(efl[totalLength])
                                else:
                                    ee11JuncF[dte] += 1/float(efl[totalLength])
                                     
                        elif (line[10][0] == '1' and not line[10][1] == '1') or (line[10][0] == line[10][1]):
                            
                            #Unfiltered
                            if not eigJunc.has_key(dte):
                                eigJunc[dte] = 1/float(efl[totalLength])
                            else:
                                eigJunc[dte] += 1/float(efl[totalLength])
                                        
                            #Filtered
                            if int(line[6]) == 3 or int(line[6]) == 2:
                                if not eigJuncF.has_key(dte):
                                    eigJuncF[dte] = 1/float(efl[totalLength])
                                else:
                                    eigJuncF[dte] += 1/float(efl[totalLength])
                        
                        else:      
                              
                            #Unfiltered
                            if not eiJunc.has_key(dte):
                                eiJunc[dte] = 1/float(efl[totalLength])
                            else:
                                eiJunc[dte] += 1/float(efl[totalLength])
                                        
                            if not eeJunc.has_key(dte):
                                eeJunc[dte] = 1/float(efl[totalLength])
                            else:
                                eeJunc[dte] += 1/float(efl[totalLength])
                                    
                            #Filtered
                            if int(line[6]) == 3 or int(line[6]) == 1:
                                if not eiJuncF.has_key(dte):
                                    eiJuncF[dte] = 1/float(efl[totalLength])
                                else:
                                    eiJuncF[dte] += 1/float(efl[totalLength])
                                        
                                if not eeJuncF.has_key(dte):
                                    eeJuncF[dte] = 1/float(efl[totalLength])
                                else:
                                    eeJuncF[dte] += 1/float(efl[totalLength])
                      
                    if line[8] == "protein_coding":
                        line[8] = line[8]+"_exon"
                                
                    if not regions.has_key("exon"):
                        regions["exon"] = 1
                    else:
                        regions["exon"] += 1
                            
                elif line[9] == "intron":
                     
                    totalLength = abs(int(line[4])) + abs(int(line[5]))
                            
                    dts = int(line[4])
                    dte = int(line[5])
                            
                    if not dts > 200: 
                                
                        #Unfiltered
                        if not eiJunc.has_key(dts):
                            eiJunc[dts] = 1/float(ifl[totalLength])
                        else:
                            eiJunc[dts] += 1/float(ifl[totalLength])
                                    
                        if not iiJunc.has_key(dts):
                            iiJunc[dts] = 1/float(ifl[totalLength])
                        else:
                            iiJunc[dts] += 1/float(ifl[totalLength])
                                    
                        #Filtered
                        if int(line[6]) == 3 or int(line[6]) == 2:
                            if not eiJuncF.has_key(dts):
                                eiJuncF[dts] = 1/float(ifl[totalLength])
                            else:
                                eiJuncF[dts] += 1/float(ifl[totalLength])
                                    
                            if not iiJuncF.has_key(dts):
                                iiJuncF[dts] = 1/float(ifl[totalLength])
                            else:
                                iiJuncF[dts] += 1/float(ifl[totalLength])
                                
                                                 
                    if not int(dte) < -200:
                                
                        #Unfiltered
                        if not ieJunc.has_key(dte):
                            ieJunc[dte] = 1/float(ifl[totalLength])
                        else:
                            ieJunc[dte] += 1/float(ifl[totalLength])
                                    
                        if not iiJunc.has_key(dte):
                            iiJunc[dte] = 1/float(ifl[totalLength])
                        else:
                            iiJunc[dte] += 1/float(ifl[totalLength])
                                    
                        #Filtered
                        if int(line[6]) == 3 or int(line[6]) == 1:
                            if not ieJuncF.has_key(dte):
                                ieJuncF[dte] = 1/float(ifl[totalLength])
                            else:
                                ieJuncF[dte] += 1/float(ifl[totalLength])
                                    
                            if not iiJuncF.has_key(dte):
                                iiJuncF[dte] = 1/float(ifl[totalLength])
                            else:
                                iiJuncF[dte] += 1/float(ifl[totalLength])
                                      
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
            for i in range(-200, 200):
                
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
                        
                if not igeJunc.has_key(i):
                    igeJunc[i] = 0         
                if not igeJuncF.has_key(i):
                    igeJuncF[i] = 0 
                              
                if not eigJunc.has_key(i):
                    eigJunc[i] = 0       
                if not eigJuncF.has_key(i):
                    eigJuncF[i] = 0
                    
                if not es11Junc.has_key(i):
                    es11Junc[i] = 0     
                if not es11JuncF.has_key(i):
                    es11JuncF[i] = 0 
                              
                if not ee11Junc.has_key(i):
                    ee11Junc[i] = 0
                if not ee11JuncF.has_key(i):
                    ee11JuncF[i] = 0
                
            #Sort for plotting against   
            eiJunc   = OrderedDict(sorted(eiJunc.items(), key=lambda x: x[0]))
            ieJunc   = OrderedDict(sorted(ieJunc.items(), key=lambda x: x[0]))
            eeJunc   = OrderedDict(sorted(eeJunc.items(), key=lambda x: x[0]))
            iiJunc   = OrderedDict(sorted(iiJunc.items(), key=lambda x: x[0]))
            igeJunc  = OrderedDict(sorted(igeJunc.items(), key=lambda x: x[0]))
            eigJunc  = OrderedDict(sorted(eigJunc.items(), key=lambda x: x[0]))
            es11Junc = OrderedDict(sorted(es11Junc.items(), key=lambda x: x[0]))
            ee11Junc = OrderedDict(sorted(ee11Junc.items(), key=lambda x: x[0]))
                
            eiJuncF   = OrderedDict(sorted(eiJuncF.items(), key=lambda x: x[0]))
            ieJuncF   = OrderedDict(sorted(ieJuncF.items(), key=lambda x: x[0]))
            eeJuncF   = OrderedDict(sorted(eeJuncF.items(), key=lambda x: x[0]))
            iiJuncF   = OrderedDict(sorted(iiJuncF.items(), key=lambda x: x[0]))
            igeJuncF  = OrderedDict(sorted(igeJuncF.items(), key=lambda x: x[0]))
            eigJuncF  = OrderedDict(sorted(eigJuncF.items(), key=lambda x: x[0]))
            es11JuncF = OrderedDict(sorted(es11JuncF.items(), key=lambda x: x[0]))
            ee11JuncF = OrderedDict(sorted(ee11JuncF.items(), key=lambda x: x[0]))
                
            regions   = OrderedDict(sorted(regions.items(), key=lambda x: (-x[1], x[0])))
            types     = OrderedDict(sorted(types.items(), key=lambda x: (-x[1], x[0])))
                          
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
                
            for k in igeJunc:
                igeJuncDist.append(k)
                igeJuncCount.append(igeJunc[k])
                 
            for k in eigJunc:
                eigJuncDist.append(k)
                eigJuncCount.append(eigJunc[k])
                
            for k in es11Junc:
                es11JuncDist.append(k)
                es11JuncCount.append(es11Junc[k])
                 
            for k in ee11Junc:
                ee11JuncDist.append(k)
                ee11JuncCount.append(ee11Junc[k])    
                
                            
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
                
            for k in igeJuncF:
                igeJuncDistF.append(k)
                igeJuncCountF.append(igeJuncF[k])
                 
            for k in eigJuncF:
                eigJuncDistF.append(k)
                eigJuncCountF.append(eigJuncF[k])
                
            for k in es11JuncF:
                es11JuncDistF.append(k)
                es11JuncCountF.append(es11JuncF[k])
                 
            for k in ee11JuncF:
                ee11JuncDistF.append(k)
                ee11JuncCountF.append(ee11JuncF[k])  
                    
            for k in regions:
                region.append(k)
                regionCount.append(regions[k])
                    
            for k in types:
                vType.append(k)
                typeCount.append(types[k])
                     
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
            
            eigData = {
                'eigJuncDist'   : eigJuncDist,
                'eigJuncCounts' : eigJuncCount      
            }
            
            igeData = {
                'igeJuncDist'   : igeJuncDist,
                'igeJuncCounts' : igeJuncCount       
            }
            
            eigDataF = {
                'eigJuncDist'   : eigJuncDistF,
                'eigJuncCounts' : eigJuncCountF      
            }
            
            igeDataF = {
                'igeJuncDist'   : igeJuncDistF,
                'igeJuncCounts' : igeJuncCountF       
            }
            
            es11Data = {
                'es11JuncDist'   : es11JuncDist,
                'es11JuncCounts' : es11JuncCount      
            }
            
            ee11Data = {
                'ee11JuncDist'   : ee11JuncDist,
                'ee11JuncCounts' : ee11JuncCount       
            }
            
            es11DataF = {
                'es11JuncDist'   : es11JuncDistF,
                'es11JuncCounts' : es11JuncCountF      
            }
            
            ee11DataF = {
                'ee11JuncDist'   : ee11JuncDistF,
                'ee11JuncCounts' : ee11JuncCountF       
            }
            
            
                
    #         typeData = {
    #             'type' : vType,
    #             'typeCounts' : typeCount       
    #         }
    #         
    #         regionData = {
    #             'region' : region,
    #             'regionCounts' : regionCount            
    #         }
                     
                          
            #Toolbar for each plot       
            TOOLS="resize,pan,wheel_zoom,box_zoom,reset,save"
                 
                 
            #Junction
            output_file(path+".html")
                
            ##########################################################################################################
            #Unfiltered
            #Intron-Exon Junction
            #---------------------------------------------------------------------------------------------------------
            ieP = figure(title=heading+"_ieJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
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
            #--------------------------------------------------------VM-------------------------------------------------
                 
            #Exon-Intron Junction
            #---------------------------------------------------------------------------------------------------------
            eiP = figure(title=heading+"_eiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
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
            eeP = figure(title=heading+"_eeJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
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
            iiP = figure(title=heading+"_iiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
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
            
            #Intergenic-Exon Junction
            #---------------------------------------------------------------------------------------------------------
            igeP = figure(title=heading+"_igeJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            igeP.line(igeJuncDist, igeJuncCount)
                 
            igeData_table = DataTable(source=ColumnDataSource(igeData), columns=[TableColumn(field="igeJuncDist", title="Distances"), TableColumn(field="igeJuncCounts", title="Count")], width=600, height=600)
                 
            igePlot = vplot(
                    hplot(igeP, igeData_table)
            )
            tab5 = Panel(child=igePlot, title="Intergenic-Exon Junction")
                 
            igeDataOutput = open(path+"_igeData.txt", "w")
            igeDataOutput.write("Distances\tCount\n")
            for i in range(0, len(igeJuncDist)):
                igeDataOutput.write(str(igeJuncDist[i]) + "\t" + str(igeJuncCount[i]) + "\n")
                 
            igeDataOutput.close()
            
            
            #Exon-Intergenic Junction
            #---------------------------------------------------------------------------------------------------------
            eigP = figure(title=heading+"_eigJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            eigP.line(eigJuncDist, eigJuncCount)
                 
            eigData_table = DataTable(source=ColumnDataSource(eigData), columns=[TableColumn(field="eigJuncDist", title="Distances"), TableColumn(field="eigJuncCounts", title="Count")], width=600, height=600)
                 
            eigPlot = vplot(
                    hplot(eigP, eigData_table)
            )
            tab6 = Panel(child=eigPlot, title="Exon-Intergenic Junction")
                 
            eigDataOutput = open(path+"_eigData.txt", "w")
            eigDataOutput.write("Distances\tCount\n")
            for i in range(0, len(eigJuncDist)):
                eigDataOutput.write(str(eigJuncDist[i]) + "\t" + str(eigJuncCount[i]) + "\n")
                 
            eigDataOutput.close()
            
            
            #Exon-Exon 1/1 Start Junction
            #---------------------------------------------------------------------------------------------------------
            es11P = figure(title=heading+"_e-1/1-s-Junction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            es11P.line(es11JuncDist, es11JuncCount)
                 
            es11Data_table = DataTable(source=ColumnDataSource(es11Data), columns=[TableColumn(field="es11JuncDist", title="Distances"), TableColumn(field="es11JuncCounts", title="Count")], width=600, height=600)
                 
            es11Plot = vplot(
                    hplot(es11P, es11Data_table)
            )
            
            tab7 = Panel(child=es11Plot, title="Exon-1/1 Start") 
            
                 
            es11DataOutput = open(path+"_es11Data.txt", "w")
            es11DataOutput.write("Distances\tCount\n")
            for i in range(0, len(es11JuncDist)):
                es11DataOutput.write(str(es11JuncDist[i]) + "\t" + str(es11JuncCount[i]) + "\n")
                 
            es11DataOutput.close()
            
            #Exon-Exon 1/1 End Junction
            #---------------------------------------------------------------------------------------------------------
            ee11P = figure(title=heading+"_e-1/1-e-Junction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            ee11P.line(ee11JuncDist, ee11JuncCount)
                 
            ee11Data_table = DataTable(source=ColumnDataSource(ee11Data), columns=[TableColumn(field="ee11JuncDist", title="Distances"), TableColumn(field="ee11JuncCounts", title="Count")], width=600, height=600)
                 
            ee11Plot = vplot(
                    hplot(ee11P, ee11Data_table)
            )
            
            tab8 = Panel(child=ee11Plot, title="Exon-1/1 End")
                 
            ee11DataOutput = open(path+"_ee11Data.txt", "w")
            ee11DataOutput.write("Distances\tCount\n")
            for i in range(0, len(ee11JuncDist)):
                ee11DataOutput.write(str(ee11JuncDist[i]) + "\t" + str(ee11JuncCount[i]) + "\n")
                 
            ee11DataOutput.close()
            
            
            #---------------------------------------------------------------------------------------------------------
            ##########################################################################################################
                 
            ##########################################################################################################
            #Filtered
            #Intron-Exon Junction
            #---------------------------------------------------------------------------------------------------------
            iePF = figure(title=heading+"_ieJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            iePF.line(ieJuncDistF, ieJuncCountF)
                 
            ieData_tableF = DataTable(source=ColumnDataSource(ieDataF), columns=[TableColumn(field="ieJuncDist", title="Distances"), TableColumn(field="ieJuncCounts", title="Count")], width=600, height=600)
                 
            iePlotF = vplot(
                    hplot(iePF, ieData_tableF)
            )
            tab9 = Panel(child=iePlotF, title="Intron-Exon Junction")
                 
            ieFDataOutput = open(path+"_ieFData.txt", "w")
            ieFDataOutput.write("Distances\tCount\n")
            for i in range(0, len(ieJuncDistF)):
                ieFDataOutput.write(str(ieJuncDistF[i]) + "\t" + str(ieJuncCountF[i]) + "\n")
                 
            ieFDataOutput.close()
                     
            #---------------------------------------------------------------------------------------------------------
                 
            #Exon-Intron Junction
            #---------------------------------------------------------------------------------------------------------
            eiPF = figure(title=heading+"_eiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            eiPF.line(eiJuncDistF, eiJuncCountF)
                 
            eiData_tableF = DataTable(source=ColumnDataSource(eiDataF), columns=[TableColumn(field="eiJuncDist", title="Distances"), TableColumn(field="eiJuncCounts", title="Count")], width=600, height=600)
                 
            eiPlotF = vplot(
                    hplot(eiPF, eiData_tableF)
            )
            tab10 = Panel(child=eiPlotF, title="Exon-Intron Junction")
                 
            eiFDataOutput = open(path+"_eiFData.txt", "w")
            eiFDataOutput.write("Distances\tCount\n")
            for i in range(0, len(eiJuncDistF)):
                eiFDataOutput.write(str(eiJuncDistF[i]) + "\t" + str(eiJuncCountF[i]) + "\n")
                 
            eiFDataOutput.close()
            #---------------------------------------------------------------------------------------------------------
                 
            #Exon-Exon Junction
            #---------------------------------------------------------------------------------------------------------
            eePF = figure(title=heading+"_eeJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            eePF.line(eeJuncDistF, eeJuncCountF)
                 
            eeData_tableF = DataTable(source=ColumnDataSource(eeDataF), columns=[TableColumn(field="eeJuncDist", title="Distances"), TableColumn(field="eeJuncCounts", title="Count")], width=600, height=600)
                 
            eePlotF = vplot(
                    hplot(eePF, eeData_tableF)
            )
            tab11 = Panel(child=eePlotF, title="Exon-Exon Junction")
                 
            eeFDataOutput = open(path+"_eeFData.txt", "w")
            eeFDataOutput.write("Distances\tCount\n")
            for i in range(0, len(eeJuncDistF)):
                eeFDataOutput.write(str(eeJuncDistF[i]) + "\t" + str(eeJuncCountF[i]) + "\n")
                 
            eeFDataOutput.close()
            #---------------------------------------------------------------------------------------------------------
                 
            #Intron-Intron Junction
            #---------------------------------------------------------------------------------------------------------
            iiPF = figure(title=heading+"_iiJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            iiPF.line(iiJuncDistF, iiJuncCountF)
                 
            iiData_tableF = DataTable(source=ColumnDataSource(iiDataF), columns=[TableColumn(field="iiJuncDist", title="Distances"), TableColumn(field="iiJuncCounts", title="Count")], width=600, height=600)
                 
            iiPlotF = vplot(
                    hplot(iiPF, iiData_tableF)
            )
            tab12 = Panel(child=iiPlotF, title="Intron-Intron Junction")
                 
            iiFDataOutput = open(path+"_iiFData.txt", "w")
            iiFDataOutput.write("Distances\tCount\n")
            for i in range(0, len(iiJuncDistF)):
                iiFDataOutput.write(str(iiJuncDistF[i]) + "\t" + str(iiJuncCountF[i]) + "\n")
                 
            iiFDataOutput.close()
                 
            #Intergenic-Exon Junction
            #---------------------------------------------------------------------------------------------------------
            igePF = figure(title=heading+"_igeJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            igePF.line(igeJuncDistF, igeJuncCountF)
                 
            igeData_tableF = DataTable(source=ColumnDataSource(igeDataF), columns=[TableColumn(field="igeJuncDist", title="Distances"), TableColumn(field="igeJuncCounts", title="Count")], width=600, height=600)
                 
            igePlotF = vplot(
                    hplot(igePF, igeData_tableF)
            )
            tab13 = Panel(child=igePlotF, title="Intergenic-Exon Junction")
                 
            igeFDataOutput = open(path+"_igeFData.txt", "w")
            igeFDataOutput.write("Distances\tCount\n")
            for i in range(0, len(igeJuncDistF)):
                igeFDataOutput.write(str(igeJuncDistF[i]) + "\t" + str(igeJuncCountF[i]) + "\n")
                 
            igeFDataOutput.close()
            
            
            #Exon-Intergenic Junction
            #---------------------------------------------------------------------------------------------------------
            eigPF = figure(title=heading+"_eigJunction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            eigPF.line(eigJuncDistF, eigJuncCountF)
                 
            eigData_tableF = DataTable(source=ColumnDataSource(eigDataF), columns=[TableColumn(field="eigJuncDist", title="Distances"), TableColumn(field="eigJuncCounts", title="Count")], width=600, height=600)
                 
            eigPlotF = vplot(
                    hplot(eigPF, eigData_tableF)
            )
            tab14 = Panel(child=eigPlotF, title="Exon-Intergenic Junction")
                 
            eigFDataOutput = open(path+"_eigData.txt", "w")
            eigFDataOutput.write("Distances\tCount\n")
            for i in range(0, len(eigJuncDistF)):
                eigFDataOutput.write(str(eigJuncDistF[i]) + "\t" + str(eigJuncCountF[i]) + "\n")
                 
            eigFDataOutput.close()
                   
            #Exon-Exon 1/1 Start Junction
            #---------------------------------------------------------------------------------------------------------
            es11PF = figure(title=heading+"_e-1/1-s-Junction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            es11PF.line(es11JuncDistF, es11JuncCountF)
                 
            es11Data_tableF = DataTable(source=ColumnDataSource(es11DataF), columns=[TableColumn(field="es11JuncDist", title="Distances"), TableColumn(field="es11JuncCounts", title="Count")], width=600, height=600)
                 
            es11PlotF = vplot(
                    hplot(es11PF, es11Data_tableF)
            )
            tab15 = Panel(child=es11PlotF, title="Exon-1/1 Start")
                 
            es11FDataOutput = open(path+"_es11Data.txt", "w")
            es11FDataOutput.write("Distances\tCount\n")
            for i in range(0, len(es11JuncDistF)):
                es11FDataOutput.write(str(es11JuncDistF[i]) + "\t" + str(es11JuncCountF[i]) + "\n")
                 
            es11FDataOutput.close()
            
            #Exon-Exon 1/1 End Junction
            #---------------------------------------------------------------------------------------------------------
            ee11PF = figure(title=heading+"_e-1/1-e-Junction", tools=TOOLS, x_axis_label='Distances', y_axis_label='Normalized Counts')
            ee11PF.line(ee11JuncDistF, ee11JuncCountF)
                 
            ee11Data_tableF = DataTable(source=ColumnDataSource(ee11DataF), columns=[TableColumn(field="ee11JuncDist", title="Distances"), TableColumn(field="ee11JuncCounts", title="Count")], width=600, height=600)
                 
            ee11PlotF = vplot(
                    hplot(ee11PF, ee11Data_tableF)
            )
            
            tab16 = Panel(child=ee11PlotF, title="Exon-1/1 End")
                 
            ee11FDataOutput = open(path+"_ee11Data.txt", "w")
            ee11FDataOutput.write("Distances\tCount\n")
            for i in range(0, len(ee11JuncDistF)):
                ee11FDataOutput.write(str(ee11JuncDistF[i]) + "\t" + str(ee11JuncCountF[i]) + "\n")
                 
            ee11FDataOutput.close()
            
            #########################################################################################################    
                
    #         Types and Regions
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
                
            unfTabs = Tabs(tabs=[tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8])
            unfiltered = Panel(child=unfTabs, title="Junction unfiltered")
                 
            fTabs = Tabs(tabs=[tab9, tab10, tab11, tab12, tab13, tab14, tab15, tab16])
            filtered = Panel(child=fTabs, title="Junction filtered")
            jt = Tabs(tabs=[unfiltered, filtered])
                 
            save(jt)      
        #===================================================================================