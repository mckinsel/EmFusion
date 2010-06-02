#!/usr/bin/python

import sys

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
        
        return {'match_count':match_count, 'mismatch_count':mismatch_count,
                'read_id':read_id, 'mapping':mapping, 'mapping_start':mapping_start,
                'mapping_end':mapping_end, 'read_length':read_length}
        
        

def main(blat_results_file, max_distance = '15000'):

    bi = psl_iterator(blat_results_file)
    
    max_dist = int(max_distance)

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
    
    for read in read_alns_1:
        maps1 = read_alns_1[read]
        
        try:
            maps2 = read_alns_2[read]
        except KeyError:
            continue
        
        concordant = False

        for map1 in maps1:
            for map2 in maps2:
                if map1[0] == map2[0] and abs(map1[1] - map2[1]) < max_dist:
                    concordant = True
        if concordant:
           print read

if __name__ == '__main__':
    main(*sys.argv[1:])

