__author__ = 'Nadia'

# --------------------------------------------------
# featureCLIP class
# Authors: Nadia Ashraf, nadia.ashraf@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: Aug 2016
# --------------------------------------------------

from bokeh.charts import HeatMap, output_file, show
from bokeh.palettes import YlOrRd3 as palette
from density import extract
import subprocess


class HeatMap:

    fInput = ''
    sInput = ''
    fOutput = ''
    sOutput = ''
    width = 5
    height = 5
    feat = ''

    def __init__(self,options):

        if hasattr(options, 'input'):
            if len(options.input) == 1:
                self.fInput = options.input[0]
            else:
                self.fInput = options.input[0]
                self.sInput = options.input[1]

        if hasattr(options, 'output'):
            self.fOutput = options.output[0]
            self.sOutput = options.output[1]

        if hasattr(options,'height'):
            self.height = options.height

        if hasattr(options,'width'):
            self.width = options.width

        if hasattr(options,'element'):
            self.feat = options.element


    def heatmap(self):
        extract(self.fOutput,self.fInput,self.sInput,self.feat)
        command = 'Rscript'
        args = [self.fOutput,self.sOutput,str(self.width),str(self.height)]
        cmd = [command,'Heatmap.R']+args
        subprocess.call(cmd)





