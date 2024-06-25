from matplotlib import pyplot
from matplotlib_venn import venn3_unweighted, venn2_unweighted, venn3_circles, venn2_circles
import numpy, os
from math import log2

class indiv_venn:
    def __init__(this, sample):
        this.data = {
            '1':0,      '2':0,      '3':0,      '4':0,      '123':0,
            '12':0,     '13':0,     '14':0,     '23':0,     '24':0,     '34':0,
            } #where 1 is arg, 2 is arg-mge, 3 is mge, 4 is pacbio
        this.groupCount = {}
        this.sample = sample

    def addToCount(this, filepath, index):
        sampleFile = open(filepath, "r")
        lineNum = 0
        for line in sampleFile:
            lineNum += 1
            if lineNum < 20:
                continue
            line_list = line.split(',')
            gene_header = line_list[0].split('|')
            if gene_header[4] not in this.groupCount:
                this.groupCount[gene_header[4]]=[False,False,False,False]
            this.groupCount[gene_header[4]][index-1] = True

    def findFinalCount(this):
        for geneCount in this.groupCount.values():
            trueIndex = []
            for i in range(0,4):
                if geneCount[i]:
                    trueIndex.append(i+1)
            if len(trueIndex) == 4:
                for intersection in this.data:
                    this.data[intersection] += 1
            elif len(trueIndex) > 0:
                this.data[str(trueIndex[0])] += 1
                if len(trueIndex) > 1:
                    this.data[str(trueIndex[1])] += 1
                    this.data[str(trueIndex[0])+str(trueIndex[1])] += 1
                    if len(trueIndex) > 2:
                        this.data[str(trueIndex[2])] += 1
                        this.data[str(trueIndex[0])+str(trueIndex[2])] += 1
                        this.data[str(trueIndex[1])+str(trueIndex[2])] += 1
                        if trueIndex[2] != 4:
                            this.data[str(trueIndex[0])+str(trueIndex[1])+str(trueIndex[2])] += 1

    def makeFigure(this, outputFolder, VENN):
        def makeVenn3():
            dataForVenn = {
                '1':0,      '2':0,      '3':0,      '123':0,
                '12':0,     '13':0,     '23':0
                }
            dataForVenn['123'] = this.data['123']
            dataForVenn['12'] = this.data['12'] - dataForVenn['123']
            dataForVenn['13'] = this.data['13'] - dataForVenn['123']
            dataForVenn['23'] = this.data['23'] - dataForVenn['123']
            dataForVenn['1'] = this.data['1'] - dataForVenn['12'] - dataForVenn['13'] - dataForVenn['123']
            dataForVenn['2'] = this.data['2'] - dataForVenn['12'] - dataForVenn['23'] - dataForVenn['123']
            dataForVenn['3'] = this.data['3'] - dataForVenn['13'] - dataForVenn['23'] - dataForVenn['123']
            return dataForVenn
        def makeVenn2(index):
            dataForVenn = {
                index:0,    '4':0,       index+'4':0
                }
            dataForVenn[index+'4'] = this.data[index+'4']
            dataForVenn[index] = this.data[index] - dataForVenn[index+'4']
            dataForVenn['4'] = this.data['4'] - dataForVenn[index+'4']
            return dataForVenn
        gs_kw = dict(width_ratios=[3,1],height_ratios=[1,1,1])
        fig, axs = pyplot.subplot_mosaic([['l','ur'],['l','cr'],['l','lr']], gridspec_kw=gs_kw, constrained_layout=True)
        fig.suptitle(this.sample, fontsize=20)

        pyplot.sca(axs['l'])
        data = makeVenn3()
        v = venn3_unweighted(subsets=(data['1'], data['2'], data['12'], data['3'], data['13'], data['23'], data['123']), 
                             subset_areas=(round(log2(data['1']+1)), round(log2(data['2']+1)), round(log2(data['12']+1)), round(log2(data['3']+1)), 
                                           round(log2(data['13']+1)), round(log2(data['23']+1)), round(log2(data['123']+1))),
                             set_labels=('','',''),
                             set_colors=('red', 'forestgreen', 'royalblue'), 
                             alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
                if text.get_text() == '0':
                    text.set_text('')
        c = venn3_circles(subsets=(round(log2(data['1']+1)), round(log2(data['2']+1)), round(log2(data['12']+1)), round(log2(data['3']+1)), 
                                   round(log2(data['13']+1)), round(log2(data['23']+1)), round(log2(data['123']+1))), 
                          linestyle="solid", linewidth=.5, color='black')

        if this.sample == "Mock+V2":
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(-0.275, -0.015), ha='center', fontsize=15)
        elif this.sample == "Mock+XT":
            v.get_label_by_id('110').set_text('')
            pyplot.annotate(data['12'], xy=v.get_label_by_id('110').get_position(), xytext=(-0.05, 0.35), ha='center', fontsize=15)
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(-0.05, -0.1), ha='center', fontsize=15)


        pyplot.sca(axs['ur'])
        data = makeVenn2('1')
        v = venn2_unweighted(subsets=(data['1'], data['4'], data['14']), 
                             subset_areas=(round(log2(data['1']+1)), round(log2(data['4']+1)), round(log2(data['14']+1))),
                             set_labels=('',''), 
                             set_colors=('red', 'gold'), 
                             alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
                if text.get_text() == '0':
                    text.set_text('')
        c = venn2_circles(subsets=(round(log2(data['1']+1)), round(log2(data['4']+1)), round(log2(data['14']+1))), 
                          linestyle="solid", linewidth=.5, color='black')
        
        pyplot.sca(axs['cr'])
        data = makeVenn2('2')
        v = venn2_unweighted(subsets=(data['2'], data['4'], data['24']), 
                             subset_areas=(round(log2(data['2']+1)), round(log2(data['4']+1)), round(log2(data['24']+1))),
                             set_labels=('',''), 
                             set_colors=('forestgreen', 'gold'), 
                             alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
                if text.get_text() == '0':
                    text.set_text('')
        c = venn2_circles(subsets=(round(log2(data['2']+1)), round(log2(data['4']+1)), round(log2(data['24']+1))), 
                          linestyle="solid", linewidth=.5, color='black')
        
        pyplot.sca(axs['lr'])
        data = makeVenn2('3')
        v = venn2_unweighted(subsets=(data['3'], data['4'], data['34']), 
                             subset_areas=(round(log2(data['3']+1)), round(log2(data['4']+1)), round(log2(data['34']+1))),
                             set_labels=('',''), 
                             set_colors=('royalblue', 'gold'), 
                             alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
                if text.get_text() == '0':
                    text.set_text('')
        c = venn2_circles(subsets=(round(log2(data['3']+1)), round(log2(data['4']+1)), round(log2(data['34']+1))), 
                          linestyle="solid", linewidth=.5, color='black')
        
        pyplot.gcf()
        if not(os.path.exists(outputFolder)):
            os.makedirs(outputFolder)
        pyplot.savefig(outputFolder + "/arg_" + this.sample + VENN)
        pyplot.close()
