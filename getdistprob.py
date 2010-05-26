from __future__ import division
import sys

def main(mapped_distance_file, ofile):
    
    curr_count = 0
    total_count = 0
    curr_val = None
    prob_dict = {}
    
    
    for line in file(mapped_distance_file):
        val = int(line.strip())
        if val != curr_val:
            if curr_val:
                prob_dict[curr_val] = curr_count
            curr_count = 1
            curr_val = val
        else:
            curr_count += 1
        total_count +=1
        
    for val in prob_dict:
        freq = prob_dict[val]
        prob = freq/total_count
        prob_dict[val] = prob
        
    of = open(ofile, 'w')
    
    for val in sorted(prob_dict.keys()):
        of.write('\t'.join([str(val), str(prob_dict[val])]))
        of.write('\n')
    
if __name__ == '__main__':
    main(*sys.argv[1:])