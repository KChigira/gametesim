#!/usr/bin/env python3
import math
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from gametesim import rwhap
from gametesim.params import Params
from gametesim.utils import time_stamp

pm = Params('genovisual')
args = pm.set_options()

colors =['lightgray', 'black', 'red', 'blue', 'yellow', 
         'darkgreen', 'orange', 'pink', 'aqua', 'lightgreen',
         'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray',
         'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray',
         'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray',
         'gray', 'gray']

class GenoVisual(object):
    def __init__(self, args):
        self.args = args
        self.dir = args.directory

        #Check args
        pm.genovisual_check_args(args)

        self.char_1 = args.character.split(',')

        #Make output directory
        self.out = '{}_visual'.format(self.dir)
        os.mkdir(self.out)

    def run(self):
        files = os.listdir(self.dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        parent_list = []
        sample_num = 0
        for f in files_hap_1:
            f_stem = ''.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                parent_list.append('{}/{}'.format(self.dir, f_stem))
                sample_num = sample_num + 1
        if sample_num < 1:
            print(time_stamp(),
                  'Disgnated directory does not cantain valid haploid files.',
                  flush=True)
            sys.exit(1) 

        #Get length of chromosome
        chr = []
        chr_len = []
        hap = [rwhap.readhap('{}_1.hap'.format(parent_list[0])),
               rwhap.readhap('{}_2.hap'.format(parent_list[0]))]
        format = hap[0][1]
        for i in range(len(hap[0]) - 2):
            chr.append(hap[0][i+2][0])
            seq = hap[0][i+2][1]
            chr_len.append(len(seq) * format)

        #number of digits in the length of the longest chromosome
        max_len = max(chr_len)
        digits = math.floor(math.log10(max_len))
        standard = 10**(digits)
        #if the longest chr length is 23098790, standard = 10000000
        if(max_len / standard < 2):
            standard = standard / 5
        elif(max_len / standard < 5):
            standard = int(standard / 2)
        #if the longest chr length is 23098790, standard = 5000000
        y_axis_at = range(0, standard*11, standard)
        y_axis_lab = []
        if(standard >= 100000):
            st_lab = standard/1000000
            sign = 'M'
        elif(standard >= 100):
            st_lab = standard/1000
            sign = 'K'
        else:
            st_lab = standard
            sign = 'bp'
        for i in range(11):
            y_axis_lab.append('{}{}'.format(round(st_lab * i, 1), sign))

        #Make figure
        for h in range(len(parent_list)):
            hap = [rwhap.readhap('{}_1.hap'.format(parent_list[h])),
                   rwhap.readhap('{}_2.hap'.format(parent_list[h]))]
            chr_num = len(hap[0]) - 2
            # Create a figure
            fig = plt.figure(figsize=(5, 5), dpi=144)
            ax = fig.add_subplot(111,
                                xlim=[-1, len(chr)],
                                xticks=range(len(chr)),
                                xticklabels=chr,
                                xlabel="Chromosome",
                                ylim=[max_len*1.05, -max_len*0.05],
                                yticks=y_axis_at,
                                yticklabels=y_axis_lab,
                                ylabel="Position")
            plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
            plt.xticks(rotation=45)
            plt.xlim(-1, len(chr))
            plt.ylim(max_len*1.05, -max_len*0.05)

            legends = []
            for i in range(len(self.char_1)):
                legends.append(patches.Patch(color=colors[i], 
                                             label=self.char_1[i]))
            ax.legend(handles=legends, loc="lower right")

            #each chromosome
            for i in range(chr_num):
                str_1 = hap[0][i + 2][1]
                str_2 = hap[1][i + 2][1]
                for j in range(len(str_1)):
                    x1 = -1
                    x2 = -1
                    for k in range(len(self.char_1)):
                        if str_1[j] == self.char_1[k]:
                            x1 = k
                            break
                    for k in range(len(self.char_1)):
                        if str_2[j] == self.char_1[k]:
                            x2 = k
                            break
                    if x1 == -1 or x2 == -1:
                        print(time_stamp(),
                              'Haploid files have unrecognized characterr.',
                              flush=True)
                        sys.exit(1) 
                    r = patches.Rectangle(xy=(i-0.4, j * format), width=0.4, 
                                          height=format, ec=None, 
                                          fc=colors[x1], fill=True)
                    ax.add_patch(r)
                    r = patches.Rectangle(xy=(i, j * format), width=0.4, 
                                          height=format, ec=None, 
                                          fc=colors[x2], fill=True)
                    ax.add_patch(r)

                #Draw rectangle of chromosome
                r = patches.Rectangle(xy=(i-0.4, 0), width=0.4,
                    height=chr_len[i], ec='black', fill=False)
                ax.add_patch(r)
                r = patches.Rectangle(xy=(i, 0), width=0.4,
                    height=chr_len[i], ec='black', fill=False)
                ax.add_patch(r)

            # Save figure
            file = '{}/{}.png'.format(self.out, hap[0][0])
            fig.savefig(file, dpi=144)

            # Release memory
            plt.clf()
            plt.close()


def main():
    print(time_stamp(), 'Genotype visualize started.', flush=True)

    prog = GenoVisual(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
