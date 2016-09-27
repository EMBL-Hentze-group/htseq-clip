__author__ = 'Nadia'

# --------------------------------------------------
# featureCLIP class
# Authors: Nadia Ashraf, nadia.ashraf@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: Aug 2016
# --------------------------------------------------

from collections import OrderedDict


def extract(fOut,fIn,sIn,feat):

    output = open(fOut, 'w')
    region_file = open(fIn, 'r')
    updown_file = open(sIn, 'r')

    InDensity = {}
    InDensityKeys = []
    InDensityVal = []

    UpDensity = {}
    UpDensityKeys = []
    UpDensityVal = []


    for line in region_file:


        if line.startswith('#track') or line.startswith('track') or line.startswith('Chromosome'):
            continue


        else:

            line = line.split("\n")
            line = line[0].split("\t")

            if feat in line[3]:

                if InDensity.has_key(line[3]):
                    InDensity[line[3]] += float(line[10])
                else:
                    InDensity[line[3]] = float(line[10])

    for line in updown_file:
        if line.startswith('#track') or line.startswith('track') or line.startswith('Chromosome'):
            continue

        else:

            line = line.split("\n")
            line = line[0].split("\t")
            if feat in line[3]:

                if not UpDensity.has_key(line[3]):
                    UpDensity[line[3]] = float(line[12])

                else:
                    UpDensity[line[3]] += float(line[12])




    InDensity = OrderedDict(sorted(InDensity.items(), key=lambda x: (-x[1], x[0])))
    UpDensity = OrderedDict(sorted(UpDensity.items(), key=lambda x: (-x[1], x[0])))


    for k in InDensity:
        InDensityKeys.append(k)
        InDensityVal.append(InDensity[k])

    for k in UpDensity:
        UpDensityKeys.append(k)
        UpDensityVal.append(UpDensity[k])



    seq = ('Repeat','in_region','upstream_region','downstream_region')
    output.write(str("\t").join(seq)+'\n')

    den = {}

    for k in InDensity:
        up = 'upstream_'+k
        down = 'downstream_'+k
        name = str(k)
        den[name] = [0.0,0.0,0.0]

        if InDensity.has_key(name):
            if den.has_key(name):
                den[name][1]+=InDensity[k]
            else:
                den[name][1] = InDensity[k]
        if UpDensity.has_key(up):


            if den.has_key(name):
                den[name][0] += UpDensity[up]

            else:

                den[name][0] = UpDensity[up]

        if UpDensity.has_key(down):

            if not den.has_key(name):

                den[name][2] = UpDensity[down]


            else:
                den[name][2] += UpDensity[down]


    for k in den:
        seq = (k,str(den[k][1]),str(den[k][0]),str(den[k][2]))
        output.write(str("\t").join(seq)+'\n')

