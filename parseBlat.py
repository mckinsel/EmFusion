from Bio.Blast import NCBIStandalone
import sys


def main(blat_results_file, max_distance = '15000'):
    bf = open(blat_results_file, 'r')
    bp = NCBIStandalone.BlastParser()
    bi = NCBIStandalone.Iterator(bf, bp)

    max_dist = int(max_distance)

    read_alns_1 = {}
    read_alns_2 = {}
    
    counter = 0
    while 1:

        rec = bi.next()

        if rec is None:
            break
        
        counter += 1

        read_id = rec.query
        end = read_id[-1]
        base_read_id = read_id[:-2]
        read_length = rec.query_letters

        if base_read_id not in eval('read_alns_' + end):
            eval('read_alns_' + end)[base_read_id] = []

        for aln in rec.alignments:
            mapping = aln.title

            for hit in aln.hsps:
                if hit.identities[0] >= .9 * read_length:
                    eval('read_alns_' + end)[base_read_id].append((mapping, hit.sbjct_start))
    
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

