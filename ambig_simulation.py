from __future__ import division
import random
import string
import sys

bowtie_dir =  '~/bowtie-0.10.1/'
bowtie_bin = bowtie_dir + 'bowtie'
reference_root = '/home/mckinsel/FusionData/RisingSun/refseq_index/refseq_transcripts'

SR = random.SystemRandom()

QUALITY_LIST = [25.80, 25.78, 25.77, 25.74, 25.73, 25.67, 25.75, 25.71, 25.69, 
                25.58, 25.61, 25.61, 25.65, 25.62, 25.60, 25.53, 25.54, 25.50, 
                25.46, 25.41, 25.37, 25.27, 25.21, 25.19, 25.11, 25.04, 24.96, 
                24.86, 24.75, 24.64, 20.33, 20.22, 20.13, 20.07, 20.00, 15.41, 
                15.31, 15.22, 15.13, 15.04, 14.93, 14.85, 14.77, 14.69, 14.62, 
                14.55, 14.54, 14.60, 14.82, 15.82]

QUALITY_LIST.extend([15]*100)


ERROR_PROBABILITIES = [10**(-k/10) for k in QUALITY_LIST]
QUALITY_SEQUENCE = ''.join([chr(int(k + 33)) for k in QUALITY_LIST])

print "ERROR_PROBABILITIES"
print ERROR_PROBABILITIES

print "QUALITY_SEQUENCE"
print QUALITY_SEQUENCE


def add_errors(sequence):
    newseq = []
    for i in range(len(sequence)):
        ep = ERROR_PROBABILITIES[i]
        testval = SR.random()
        if testval < ep:
            newchar = SR.choice([k for k in ['A','C','T','G'] if k != sequence[i]])
        else:
            newchar = sequence[i]
        newseq.append(newchar)
    return ''.join(newseq)
        

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

def reverse_complement(seq):
    complement = string.maketrans("AGCTagctnNRYWSMKHBVD", "TCGAtcganNYRWSKMDVBH")
    return seq.translate(complement)[::-1]
        
def parse_line(efline):
    
    out = {}
    
    lseq = efline.strip().split('\t')
    out['transcriptid'] = lseq[0]
    out['geneid'] = lseq[1]
    out['exonid'] = lseq[4]
    out['start'] = int(lseq[5])
    out['end'] = int(lseq[6])
    out['length'] = abs(out['start'] - out['end']) ##Add + 1 for Ensembl data!
    
    exon_position = int(lseq[7])
    out['position'] = exon_position

    return out

def make_fusion_seqs(transcript1, transcript2, exon_list1, exon_list2, 
                     seq_cont, readlength, insertlength, pairs):
    

    first_exon_index = random.choice(range(len(exon_list1) - 1)) ##We don't really want the last one
    second_exon_index = random.choice(range(1, len(exon_list2))) ##Nor the first one here
    
        
    upstream_exons = exon_list1[:first_exon_index + 1]
    downstream_exons = exon_list2[second_exon_index:]

    upstream_sequence = ''
    downstream_sequence = ''
    
    for ux in upstream_exons:
        upstream_sequence += seq_cont.getseq(ux)
    for dx in downstream_exons:
        downstream_sequence += seq_cont.getseq(dx)    

    ##If the read isn't long enough for this, raise LookupError
    if len(upstream_sequence) < readlength + insertlength or \
       len(downstream_sequence) < readlength + insertlength:
        raise LookupError
        
    reads = []
    
    if pairs:
        for k in range(insertlength):
            if k == 0:
                read1seq = upstream_sequence[-(readlength + k):]
                read2seq = reverse_complement(downstream_sequence[insertlength - k: insertlength - k + readlength])
            else:
                read1seq = upstream_sequence[-(readlength + k):-k]
                read2seq = reverse_complement(downstream_sequence[insertlength - k: insertlength - k + readlength])
            reads.append((add_errors(read1seq), add_errors(read2seq)))
            
    return reads
            
        
    
    
    
def main(exoninfofile, exonseqfile, readlength, insertlength, n, paralogfile = '', tag = ''):
    
    
    of1 = open('simulated_reads.fastq.1' + tag, 'w')
    of2 = open('simulated_reads.fastq.2' + tag, 'w')
    
    paraloglist = []
    if paralogfile:
        for line in file(paralogfile):
            paraloglist.append(line.strip())
    
    
    esq = exon_seq_cont(exonseqfile)
    print "Loaded exon seqs"
    
    current_transcriptid = ""
    
    transcript_dict = {}
    gene_dict = {}
    exon_list = []

    for line in file(exoninfofile):
        exon_d = parse_line(line)
        gene_dict[exon_d['transcriptid']] = exon_d['geneid']
        
        if exon_d['transcriptid'] != current_transcriptid:
            if exon_list:
                transcript_dict[current_transcriptid] = exon_list
            
            exon_list = []
            current_transcriptid = exon_d['transcriptid']
            
        
        exon_list.append(exon_d['exonid'])
        
   
    transcript_dict[current_transcriptid] = exon_list
    
    transcript_list = transcript_dict.keys()
    
    print "Loaded exon structures."
    
    i = 0
    while i < int(n):
        if i%1000 == 0: print "On pair", i
        transcript1 = random.choice(transcript_list)
        transcript2 = random.choice(transcript_list)
        
        if gene_dict[transcript1] == gene_dict[transcript2]: 
            print "Same Genes", transcript1, transcript2, gene_dict[transcript1], gene_dict[transcript2]
            continue
        
        if paraloglist:    
            if not (transcript1 in paraloglist and transcript2 in paraloglist):
                continue
            
        exon_list1 = transcript_dict[transcript1]
        exon_list2 = transcript_dict[transcript2]

        try:
            sim_reads = make_fusion_seqs(transcript1, transcript2, exon_list1,
                             exon_list2, esq, int(readlength), int(insertlength),
                             pairs = True)
        except LookupError:
            continue
            
        k = 0
        for read in sim_reads:
            of1.write('@' + transcript1 + '_' + transcript2 + '_' + tag + '-' + str(i) + '_' + str(k) +  '/1\n')
            of2.write('@' + transcript1 + '_' + transcript2 + '_' + tag + '-' + str(i) + '_' + str(k) + '/2\n')
            
            of1.write(read[0] + '\n')
            of2.write(read[1] + '\n')
            
            of1.write('+\n' + QUALITY_SEQUENCE[:len(read[0])] + '\n')
            of2.write('+\n' + QUALITY_SEQUENCE[:len(read[1])] + '\n')
            k+=1
        i += 1

if __name__ == '__main__':
    if len(sys.argv[1:]) in (5,6,7):
        main(*sys.argv[1:])
    else:
        print 'main(exoninfofile, exonseqfile, readlength, insertlength, n)'
