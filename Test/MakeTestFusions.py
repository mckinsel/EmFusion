import random
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
        

def main(transcript_structure_file, exon_seq_file, num_of_fusions):
    """
    transcript_structure_file contains information about the order of exons
    in known transcripts.
    
    exon_seq_file contains the actual sequence for the exons.
    
    num_of_fusions is the number of fusion transcripts to include in the 
    test data.
    """
    
    tsf = open(transcript_structure_file)
    transcript_to_exons = {}
    transcript_to_gene = {}
    curr_trans = ""
    exon_list = []
    
    ##First read in the transcript structures
    for line in tsf:
        ls = line.split('\t')
        transcriptid = ls[0]
        geneid= ls[1]
        exonid = ls[4]
        
        if transcriptid != curr_trans:
            if curr_trans:
                transcript_to_exons[curr_trans] = exon_list
                
            curr_trans = transcriptid
            exon_list = [exonid]
            transcript_to_gene[transcriptid] = geneid
        else:
            exon_list.append(exonid)
    transcript_to_exons[curr_trans] = exon_list
    
   
    assert set(transcript_to_exons.keys()) == set(transcript_to_gene.keys())
    
    ##Now make fusions randomly
    made_fusions = 0
    fusions = []
    
    while made_fusions < int(num_of_fusions):
        trans1 = random.choice(transcript_to_exons.keys())
        trans2 = random.choice(transcript_to_exons.keys())
        
        gene1 = transcript_to_gene[trans1]
        gene2 = transcript_to_gene[trans2]
        if gene1 == gene2:
            continue
        
        trans1_bound = random.choice(range(len(transcript_to_exons[trans1])))
        trans2_bound = random.choice(range(len(transcript_to_exons[trans2])))
        
        fusion_exons = transcript_to_exons[trans1][:trans1_bound + 1] +\
                       transcript_to_exons[trans2][trans2_bound:]
                       
        fusion_id = '|'.join([trans1+trans2,gene1+gene2,'fusion', 'fusion', 
                              str(trans1_bound + 1), str(trans2_bound + 1), '1'])
        fusions.append((fusion_id, fusion_exons))
        made_fusions += 1
        

        
    esc = exon_seq_cont(exon_seq_file)
    
    for fusion in fusions:
        sys.stdout.write('>' + fusion[0] + '\n')
        for exon in fusion[1]:
            sys.stdout.write(esc.getseq(exon))
        sys.stdout.write('\n')
        
if __name__ == '__main__':
    if len(sys.argv[1:]) != 3:
        print "usage MakeTestFusions.py <transcript_structure_file> <exon_seq_file> <num_of_fusions>"
        sys.exit(1)
    main(*sys.argv[1:])