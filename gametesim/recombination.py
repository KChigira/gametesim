import hashlib
import math
import random
import sys
from gametesim import rwhap
from gametesim.utils import time_stamp

def recombination(list):
    #Map file
    map = list[0]
    #Read haploid files
    hap = [[rwhap.readhap('{}_1.hap'.format(list[1])),
            rwhap.readhap('{}_2.hap'.format(list[1]))],
           [rwhap.readhap('{}_1.hap'.format(list[2])),
            rwhap.readhap('{}_2.hap'.format(list[2]))]]
    #check format
    if not(hap[0][0][1] == hap[0][1][1] == 
           hap[1][0][1] == hap[1][1][1]):
        print(time_stamp(),
            'Format of haploid files are not same among parents.',
            flush=True)
        sys.exit(1) 
    else:
        format = int(hap[0][0][1]) 
    
    #Output file name
    stem = list[3]
    #Seed
    seed = list[4]
    #Output directory
    outdir = list[5]
            
    for h in range(1,3):
        #make seed of random numbers
        seed_txt = '{}__{}__{}__{}'.format(outdir, h, stem, seed)
        hash = hashlib.md5(seed_txt.encode("utf-8")).hexdigest()
        #6 digits are used for seed
        hash = int(hash[:6], 16)  
        random.seed(hash)

        #Simulate recombinanation events
        events = []
        current_chr = ''
        current_hap = 0
        for j in range(1, len(map)): #map[0] is header
            if current_chr != map[j][0]:
                current_chr = map[j][0]
                current_hap = int(round(random.random()))
                events.append([current_chr, int(map[j][1]), current_hap])
            else:
                pos_start = int(map[j-1][1])
                pos_end =  int(map[j][1])
                cM_dif = float(map[j][2]) - float(map[j-1][2])
                #Check difference of posittion and cM
                if pos_end - pos_start <= 0 or cM_dif <= 0:
                    continue
                
                #Recombination events occurs based on Poisson distribution
                rand = random.random()
                c = 0
                prob = math.e ** (-cM_dif/100)
                while(rand > prob and c < 10): #limit of recombination: 10 
                    c = c + 1
                    prob = prob + ((cM_dif/100) ** c) * (math.e ** (-cM_dif/100)) / math.factorial(c)

                if c > 0:
                    focus_range = range(pos_start+1, pos_end-1)
                    cross_pos = random.sample(focus_range, c)
                    cross_pos.sort()
                    for pos in cross_pos:
                        if current_hap == 0:
                            current_hap = 1
                        else:
                            current_hap = 0
                        events.append([current_chr, pos, current_hap])

        #Apply recombination events
        new_gamete = []
        new_gamete.append('{}_{}'.format(stem, h)) #This is the name
        new_gamete.append(format)
        chrom = []
        current_chr = ''
        current_seq = ['', '']
        hap_row = 1
        state = 0
        char_cnt = 0
        for row in events:
            if row[0] != current_chr:
                if len(chrom) == 2:
                    if char_cnt < len(current_seq[0]):
                        chrom[1] = chrom[1] + current_seq[state][char_cnt:len(current_seq[0])]
                    new_gamete.append(chrom)

                current_chr = row[0]
                hap_row = hap_row + 1
                if hap[h-1][0][hap_row][0] != row[0]:
                    print(time_stamp(),
                            'Format of haploid files does not match the linkage map file.',
                            flush=True)
                    sys.exit(1) 
                current_seq = [hap[h-1][0][hap_row][1], hap[h-1][1][hap_row][1]]
                chrom = [current_chr, '']
                char_cnt = 0

            state = row[2]
            recombi_point = math.floor(row[1] / format)
            if recombi_point != char_cnt:
                chrom[1] = chrom[1] + current_seq[state][char_cnt:recombi_point]
                char_cnt = recombi_point

        if char_cnt < len(current_seq[0]):
            chrom[1] = chrom[1] + current_seq[state][char_cnt:len(current_seq[0])]
        new_gamete.append(chrom)

        rwhap.writehap(new_gamete, '{}/{}_{}.hap'.format(outdir, stem, h))
