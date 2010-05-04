import sys

class exon_seq_cont(object):
    
    def __init__(self, filename):
        
        
        self._seq_dict = {}
        
        current_id = ""
        current_seq = ""
        
        for line in file(filename):
            if line[0] == '>':
                if current_id:
                    self._seq_dict[current_id] = current_seq
                current_id = line.strip()[1:]
                current_seq = ""
            else:
                current_seq += line.strip()
                
        self._seq_dict[current_id] = current_seq
    
    def getseq(self, sid):
        return self._seq_dict[sid]
        
        
def main(exon_seq_file, fusion_exons_file):
    
    esq = exon_seq_cont(exon_seq_file)
    
    
    outf = open(fusion_exons_file + '.seqs', 'w')
    
    
    count = 0
    for line in file(fusion_exons_file):
        count += 1
        seq = ""
        ls = line.strip().split('\t')
        for exon in ls[2:]:
            seq += esq.getseq(exon)
        outf.write('>F_' + ls[0] + '_' + ls[1] + '_' + str(count) + '|F_' + str(count) +'|fusion|fusions|0|100|1\n')
        outf.write(seq)
        outf.write('\n')
        
        
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])