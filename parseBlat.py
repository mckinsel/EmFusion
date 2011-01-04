#!/usr/bin/python
import re
import sys

TRANSCRIPT_TAG = 'ucsc'

class tx_regions(object):
    
    def __init__(self, txfile, fuzz = 1500):
        
        self._chrom_dict = {}
        self._fuzz = int(fuzz)
        for line in file(txfile):
            ls = line.strip('\n').split('\t')
            self._chrom_dict.setdefault(ls[0], []).append( (int(ls[1]), int(ls[2])) )
        
        for chrom in self._chrom_dict:
            self._chrom_dict[chrom].sort(key=lambda x:x[0])
            
    def get_tx_range(self, chrom, coordinate):
        if chrom not in self._chrom_dict:
            return [coordinate - self._fuzz, coordinate + self._fuzz]
        
#        print 'chrom', chrom, 'coordinate:', coordinate
        for region in self._chrom_dict[chrom]:
#            print "testing region", region
            if region[1] < coordinate:
#                print "too early"
                continue
            elif region[0] <= coordinate:
#                print "hit!", [region[0] - self._fuzz, region[1] + self._fuzz]
                return [region[0] - self._fuzz, region[1] + self._fuzz]
            else:
#                print "failed"
                return [coordinate - self._fuzz, coordinate + self._fuzz]
            
        return [coordinate - self._fuzz, coordinate + self._fuzz]


class psl_iterator(object):
    
    def __init__(self, pslfile):
        self._file = open(pslfile, 'r')
        
    def __iter__(self):
        return self
    
    def next(self):
        nl = self._file.readline()
        
        if len(nl) == 0:
            return None
        else:
            nl_list = nl.strip('\n').split('\t')
            while len(nl_list) < 15 or nl_list[0] == 'match' or nl_list[1] == 'match':
                nl = self._file.readline()
                nl_list = nl.strip('\n').split('\t')

        match_count = int(nl_list[0])
        mismatch_count = int(nl_list[1])
        read_id = nl_list[9]
        mapping = nl_list[13]
        mapping_start = int(nl_list[15])
        mapping_end = int(nl_list[16])
        read_length = int(nl_list[10])
        
        if TRANSCRIPT_TAG in mapping:
            transcript_coordinate = int(mapping.split('_')[-1])
            mapping_start += transcript_coordinate
            mapping_end += transcript_coordinate
            mapping = '_'.join(mapping.split('_')[2:-1])

        return {'match_count':match_count, 'mismatch_count':mismatch_count,
                'read_id':read_id, 'mapping':mapping, 'mapping_start':mapping_start,
                'mapping_end':mapping_end, 'read_length':read_length}
        
        

def main(blat_results_file, tx_regions_file, fuzz = '1500'):

    bi = psl_iterator(blat_results_file)
    #print "created psl_iterator"
    txf = tx_regions(tx_regions_file, fuzz)
    #print "created tx_regions"
    

    read_alns_1 = {} ##Dictionary from read_id to list of hits
    read_alns_2 = {} ##Same for other paired end
    
    counter = 0
    while 1:

        rec = bi.next()

        if rec is None:
            break
        
        counter += 1

        read_id = rec['read_id']

        ##This is kind of specific to the format of your read ids
        ##end shoudl be '1' or '2'
        ##read_id should end with '/1' or '/2'

        end = read_id[-1]
        base_read_id = read_id[:-2]
        read_length = rec['read_length']

        if base_read_id not in eval('read_alns_' + end):
            eval('read_alns_' + end)[base_read_id] = []


        mapping = rec['mapping']

        if rec['match_count'] >= .8 * read_length:
            eval('read_alns_' + end)[base_read_id].append((mapping, rec['mapping_start']))
    
    #print "read in blat mappings"
    
    for read in read_alns_1:
        maps1 = read_alns_1[read]
        
        try:
            maps2 = read_alns_2[read]
        except KeyError:
            continue
        
        concordant = False

        for map1 in maps1:
            for map2 in maps2:
#                print map1, map2
                if map1[0] == map2[0]:
#                    print "checking", map1, map2
#                    print "is", map1[1], "in", txf.get_tx_range(map2[0], map2[1])
#                    print "or", map2[1], "in", txf.get_tx_range(map1[0], map1[1])
                    map1range = txf.get_tx_range(map1[0], map1[1])
                    map2range = txf.get_tx_range(map2[0], map2[1])
                    if map2range[0] <= map1[1] <= map2range[1]  or \
                       map1range[0] <= map2[1] <= map1range[1]:
#                           print "concordant is True!"
                           concordant = True
                           break
        if concordant:
           print read

if __name__ == '__main__':
    main(*sys.argv[1:])

