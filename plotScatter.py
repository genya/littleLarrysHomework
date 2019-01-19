# Allows one to plot many data points with many different customizations (marker size, color, legend label etc),
# without writing code. Instead the customizations are specified by the input spreadsheets, in which each row is
# a data point (run code without inputs for instructions). Documentation and variable names assume plotting of
# CRISPR score data but can be used generically.
#
# MIT License
#
# Copyright (c) 2019 Evgeni (Genya) Frenkel
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import matplotlib.pyplot as plt
import os.path

outputFileBaseName = 'scatter'

#plot spec file column names
colNames = ['Gene name', 'Plot or not?', 'Color', 'Marker', 'Size', 'Alpha', 'Layer', 'Label or not?',
            'Label font size','Legend group']  # note gene name needs to be first entry in this list or code breaks

if len(sys.argv)<3:
    print('\n' + sys.argv[0] + ' requires at least two inputs. These are gene-level log-fold changes or CRISPR scores (CS) and are provided as follows:'
          + '\npython ' +  sys.argv[0] + ' [CS treatment A vs initial] [CS treatment B vs initial] [plot spec file]'
          + '\n\nThe output is two image files called scatter.svg, scatter.png (dpi=500). To modify these output formats, edit the code at the very bottom of the script.'
          + '\n\nFormat of the first two input arguments (data files) should be tab or comma separated columns consisting of :'
          + '\nGene name \t CRISPR score '
            '\n\nAny gene name that is not present in both files will be ignored. '
            'The header of the second column will be the axis label for those data.'
          + '\n\nThe third (optional) argument provides specifications for how all or some data points should be plotted.'
          + ' The plot spec file is also tab or comma delimited, ' + str(len(colNames)) + ' columns in total:'
          + '\n' + '\t'.join(colNames)
          + '\n\nThe plot spec file needs to have column headers (in any order) exactly matching these column names, '
            'but can have additional columns (e.g. notes about genes, other data), which will be ignored. '
            '\nLikewise in the CRISPR score files, any columns beyond the first two will be ignored.'
            '\n\nIf value in "Plot or not?" column = 1, then data for that gene will be plotted. '
            'Any value other than 1 will be treated as false. '
            'Likewise value in "Label or not?" = 1 means text of gene name will be overlayed on the data point.'
            '\n\nLayer should be a number and points with higher value layer are plotted on top. If no layer specified, default is bottom layer.'
            '\n\nThe permitted values and meanings for columns Color, Marker, Size, and Alpha can be found in the matplotlib/pyplot documentation:'
            '\n https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html'
          + '\n\nThis code was written and tested for python 2.7, might not work with other versions.\n'
          )
    sys.exit()


fileAvI, fileBvI = sys.argv[1], sys.argv[2]
fileGOI = ''
if len(sys.argv)>3:
    fileGOI = sys.argv[3]

def getVals(fileXvY):
    CS_XvY = {}
    with open(fileXvY,'r') as f:
        line = f.readline()
        delim = '\t'
        if len(line.strip().split('\t'))==1:
            delim = ','
        label = line.split(delim)[1]

        for line in f:
            line = line.strip().split(delim)
            CS_XvY[line[0]] = float(line[1])

    return CS_XvY, label

#load score values
CS_A, xlabel = getVals(fileAvI)
CS_B, ylabel = getVals(fileBvI)
geneList = [g for g in CS_A if g in CS_B]

#load plot specs
GOIs = {}
layers = [-float('inf')]
layerLists = {} #layers[x] = [list of genes to be plotted as layer x]
layerSpecified = [] #list of genes with layer specified
if len(fileGOI)>0:
    with open(fileGOI,'r') as f:

        #tab or comma delimeter?
        header = f.readline().strip()
        delim = '\t'
        if len(header.strip().split('\t'))==1:
            delim = ','
        header = header.split(delim)

        #find index of relevant columns
        colInds = {x:i for i,x in enumerate(header) if x.strip() in colNames}
        for x in colNames:
            error = False
            if x not in colInds:
                print('Error: cannot find column `' + x  + '` in ' + fileGOI)
                error = True
            if error : sys.exit()

        for line in f:
            line = [x.strip() for x in line.split(delim)]

            GOIs[line[colInds['Gene name']]] = {x:line[colInds[x]] for x in colNames[1:]}

            try:
                if int(line[colInds['Layer']]) not in layers:
                    layers.append(int(line[colInds['Layer']]))
                    layerLists[int(line[colInds['Layer']])] = []
                layerLists[int(line[colInds['Layer']])].append(line[colInds['Gene name']])
                layerSpecified.append(line[colInds['Gene name']])
            except ValueError:
                print('Error: Layer column contains non-integer value in ' + fileGOI + ' for gene ' + line[colInds['Gene name']])
                sys.exit()

layers = sorted(layers)
layerLists[-float('inf')] = [g for g in geneList if g not in layerSpecified]


###plot
fig=plt.figure()
ax = plt.subplot()
#determine axes bounds
marginFactor = 1.05
xlim = [marginFactor*min([CS_A[g] for g in geneList] + [0]), marginFactor*max([CS_A[g] for g in geneList])]
ylim = [marginFactor*min([CS_B[g] for g in geneList] + [0]), marginFactor*max([CS_B[g] for g in geneList])]

###MANUALLY SET AXIS BOUNDS HERE
#xlim = [-3, 3]
#ylime = [-3, 3]

ax.plot(xlim,[0, 0],'--',linewidth=1,color='silver',zorder=0)
ax.plot([0, 0],ylim,'--',linewidth=1,color='silver',zorder=0)

legendHandles = []
legendSymbols = []
legendLabels = []

numPtsPlotted = 0
for layerInd, layerKey in enumerate(layers):
        for g in layerLists[layerKey]:

            #coordinates
            x, y = CS_A[g], CS_B[g]

            #get plot specs
            if g in GOIs: #custom specs
                plotOrNot = GOIs[g]['Plot or not?'] == '1'
                if plotOrNot:
                    alpha = float(GOIs[g]['Alpha'])
                    markerColor = GOIs[g]['Color']
                    markerShape = GOIs[g]['Marker']
                    if markerShape == '0': #LibreOffice Calc converts periods into zeros, which aren't valid plot shape
                        markerShape = '.'
                    markerSize = float(GOIs[g]['Size'])
                    legendGroup = GOIs[g]['Legend group']
                    labelOrNot = GOIs[g]['Label or not?'] == '1'
                    labelFont = float(GOIs[g]['Label font size'])
            else: #default specs
                plotOrNot = True
                alpha = 1
                markerColor = 'b'
                markerShape = '.'
                markerSize = 10
                legendGroup = ''
                labelOrNot = False
                labelFont = 0

            #add point to figure
            if plotOrNot:
                ax.scatter(x, y,
                             color=markerColor,
                             marker=markerShape,
                             alpha=alpha,
                             s=markerSize,
                             zorder=layerInd)
                numPtsPlotted += 1

                if numPtsPlotted%100==0:
                    print(str(numPtsPlotted) + ' data points plotted')

                # assign to legend group?
                if legendGroup != '' and legendGroup not in legendLabels:
                    legendSymbols.append(markerShape + markerColor)
                    legendHandles.append(plt.scatter(-1000, -1000, color=markerColor, marker=markerShape, alpha=alpha,
                                     s=markerSize*2, zorder=0))
                    legendLabels.append(legendGroup)

                #overlay gene name?
                if labelOrNot:
                    ax.text(x, y, g, fontsize=labelFont, zorder=max(layers) + 1)
#add legend
if len(legendHandles) > 0:
    ax.legend(tuple(legendHandles), tuple(legendLabels), fontsize=6)  # ,location=outside)

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim(xlim)
plt.ylim(ylim)
ax.set_aspect(aspect='equal')
plt.tight_layout()

#save plot to png, svg files
if os.path.isfile(outputFileBaseName + '.png') or os.path.isfile(outputFileBaseName + '.svg'):
    fileInd = 0
    while True:
        newOutputFileBaseName = outputFileBaseName + '_' + str(fileInd)
        if os.path.isfile(newOutputFileBaseName + '.png') or os.path.isfile(newOutputFileBaseName + '.svg'):
            fileInd += 1
            newOutputFileBaseName = outputFileBaseName + '_' + str(fileInd)
        else:
            plt.savefig(newOutputFileBaseName + '.svg')
            plt.savefig(newOutputFileBaseName + '.png',dpi=500)
            break
else:
    plt.savefig(outputFileBaseName + '.svg')
    plt.savefig(outputFileBaseName+ '.png', dpi=500)