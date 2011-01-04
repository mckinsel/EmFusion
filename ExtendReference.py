from copy import deepcopy
import networkx as nx
#import matplotlib.pyplot as plt
from pprint import pprint
import sys

MAX_INSERT_SIZE = 1000
MIN_MAPPING_COUNT = 2 
MIN_GOOD_OVERLAP = 6
PREPENDED_SEQUENCE_LENGTH = 50
POLY_A_LENGTH = 50

def good_read_length(length, mismatches):
    
    out_len = int(length)
    good_count = 0
    while out_len > 0:
        if (out_len - 1) not in mismatches:
            good_count += 1
        out_len -= 1
        if good_count >= MIN_GOOD_OVERLAP:
            break
    return out_len

def good_read_start(start, length, mismatches):
    
    base = int(start)
    out_start = int(start)
    good_count = 0
    while out_start < base + length - 1:
        if (out_start - base) not in mismatches:
            good_count += 1
        out_start += 1
        if good_count >= MIN_GOOD_OVERLAP:
            break
    return out_start

def dfs(v, g, firstgenename, secondgenename, mapped_reads_1, mapped_reads_2, 
        exon_bounds_1, exon_bounds_2, outfileh, read_status = [], 
        exon_list = [], length_in_1 = 0):
            
    fail = False
    #mapped_reads are lists of tuples:
    #(read_id, start_position, end_position)
    
    #read_status is a list of lists:
    #[read_id, current_implied_insert_size, closed?]
    
    #exon_bounds are dictionaries
    #eb[exon_id] = (exon_start, exon_end)
    
    #v is the name of the exon we're at
#    print '*'*80
#    print "Arrived at node " + v
#    print "Read_status "
#    pprint(read_status)
#    print "Exon list " 
#    pprint(exon_list)
    
    if sorted([k[0] for k in mapped_reads_1]) != sorted([k[0] for k in mapped_reads_2]):
        pprint(sorted(mapped_reads_1))
        pprint(sorted(mapped_reads_2))
    assert sorted([k[0] for k in mapped_reads_1]) == \
           sorted([k[0] for k in mapped_reads_2])
    
    #First we're going to check for some fail conditions:
    
    seen_reads = sorted([k[0] for k in read_status]) #List of reads we've seen 
                                                     #so far
    all_reads = sorted([k[0] for k in mapped_reads_1]) #List of all reads
    
    open_reads = sorted([k[0] for k in read_status if k[2] == 0])
    
    current_gene = g.node[v]['graph'] #The gene of the exon we're at
    
    working_read_status = deepcopy(read_status)
    
#    print 'exon_bounds_1', exon_bounds_1
#    print 'exon_bounds_2', exon_bounds_2
    if 'source' in v or 'sink' in v:
        exon_start, exon_end = (-1, -1)
    elif current_gene == 1:
        exon_start, exon_end = exon_bounds_1[v[2:]]
    elif current_gene == 2:
        exon_start, exon_end = exon_bounds_2[v[2:]]

    at_end_of_first_gene = False
    if current_gene == 1:
        if exon_end == max([k[1] for k in exon_bounds_1.values()]):
            at_end_of_first_gene = True
    
    at_start_of_second_gene = current_gene == 2 and exon_start == 0

#    print "seen reads"
#    pprint(seen_reads)
#    print "open reads"
#    pprint(open_reads)
#    print "current gene", current_gene
    
    #If we've gone to the second graph without hitting all the reads mapping to
    #the first, stop
    if current_gene == 2 and seen_reads != all_reads:
#        print "FAILING BECAUSE DIDN'T SEE ALL READS IN GENE1"
        fail = True
        
    #Now see if the implied insert size for any of the reads has become too
    #large:
    if working_read_status and max(k[1] for k in working_read_status) > \
                                                            MAX_INSERT_SIZE:
#        print "FAILING BECAUSE INSERT SIZE", max(k[1] for k in working_read_status)
#        working_read_status.sort(key=lambda x: x[1], reverse=True)
#        print "FOR GENE", working_read_status[:10]
        fail = True
    
    
    #So we didn't fail! Now, if we're in gene 1, see if we should open up any
    #reads and then add the length of the exon to the implied insert size of the
    #already open reads.
    if not fail and current_gene == 1:
        for read_tuple in mapped_reads_1:
            if read_tuple[2] < exon_end and read_tuple[0] not in seen_reads:
                if at_end_of_first_gene: ##If we're at the end of the first gene, don't add to the 
                                         ##gap size. This increases read-through sensitivity.
                    working_read_status.append([read_tuple[0], 0, 0])
                else:
                    working_read_status.append(
                                [read_tuple[0], exon_end - read_tuple[2], 0]
                                )
#                print "opening read " + read_tuple[0]
            elif read_tuple[0] in seen_reads:
                status_index = [k[0] for k in working_read_status].index(read_tuple[0])
                if at_end_of_first_gene:
                    pass
                else:
                    working_read_status[status_index][1] += (exon_end - exon_start)
    #If we're in gene 2, close the open reads that we hit and add the length of
    #the exon to those we don't
    elif not fail and current_gene == 2:
        for read_tuple in mapped_reads_2:
            if exon_start <= read_tuple[1] < exon_end:
                status_index = [k[0] for k in working_read_status].index(read_tuple[0])
                if at_start_of_second_gene:
                    pass
                else:
                    working_read_status[status_index][1] += read_tuple[1] - exon_start
                working_read_status[status_index][2] = 1
#                print "closing read " + read_tuple[0]
            elif read_tuple[0] in open_reads:
                status_index = [k[0] for k in working_read_status].index(read_tuple[0])
                if at_start_of_second_gene: ##Again with not penalizing read-through ends
                    pass
                else:
                    working_read_status[status_index][1] += (exon_end - exon_start)
    
    #Now, we see if we've come to the end of the line in gene 2
    if not fail and v == '2-' + secondgenename + 'sink':
#        print "at sink"
#        print "open_reads", open_reads
        #Make sure no reads are left open
#        g1_exon_count = len([k for k in g.nodes() if k[:2] == '1-']) - 1 #Won't include sink
#        g2_exon_count = len([k for k in g.nodes() if k[:2] == '2-']) - 2 #Won't include source or sink
        
#        exons_from_g1 = len([k for k in exon_list if k[:2] == '1-'])
#        exons_from_g2 = len([k for k in exon_list if k[:2] == '2-'])
#        print "all nodes", g.nodes()
#        print 'g1_ec', [k for k in g.nodes() if k[:2] == '1-']
#        print 'g2_ec', [k for k in g.nodes() if k[:2] == '2-']
#        print 'ec1', [k for k in exon_list if k[:2] == '1-']
#        print 'ec2', [k for k in exon_list if k[:2] == '2-']

        if len(open_reads) == 0: # and g1_exon_count != exons_from_g1 and g2_exon_count != exons_from_g2:
            outfileh.write(firstgenename + '.' + secondgenename + '\t' + \
                            str(length_in_1 - PREPENDED_SEQUENCE_LENGTH) + '\t')
            outfileh.write('\t'.join([k[2:] for k in exon_list[1:]]))
            outfileh.write('\n')
            
    
    
    #Otherwise, just keep moving
    if not fail:
        if current_gene == 1:
            new_length_in_1 = length_in_1 + (exon_end - exon_start)
        elif current_gene == 2:
            new_length_in_1 = int(length_in_1)
        
        new_exon_list = deepcopy(exon_list)
        new_exon_list.append(v)
        for neighbor in g.neighbors(v):
            dfs(neighbor, g, firstgenename, secondgenename, mapped_reads_1, 
            mapped_reads_2,exon_bounds_1, exon_bounds_2, outfileh, 
            working_read_status, new_exon_list, new_length_in_1)



def write_exon_orders(eg1, eg2, gene1, gene2, outfileh, mapped_reads_1, 
                      mapped_reads_2, exon_bounds_1, exon_bounds_2):
    
#    print gene1, gene2
#    print set(eg1.nodes()).intersection(eg2.nodes())
    g1 = nx.union(eg1, eg2, rename = ('1-', '2-'))
    
    for node1 in eg1.nodes():
        g1.node['1-' + node1] = eg1.node[node1]
    for node2 in eg2.nodes():
        g1.node['2-' + node2] = eg2.node[node2]   
        
    for node1 in eg1.nodes():
        g1.node['1-' + node1]['graph'] = 1
        
    for node2 in eg2.nodes():
        g1.node['2-' + node2]['graph'] = 2 
        
    for node1 in eg1.nodes():
        for node2 in eg2.nodes():
            if node1 not in (gene1 + "source", gene1 + "sink") and \
               node2 not in (gene2 + "source", gene2 + "sink"):
                g1.add_edge('1-' + node1, '2-' + node2)
                
#    exon_bounds_1 = exon_dict[gene1]
#    exon_bounds_2 = exon_dict[gene2]
        
    dfs('1-' + gene1 + "source", g1, gene1, gene2, mapped_reads_1, mapped_reads_2,
        exon_bounds_1, exon_bounds_2, outfileh)

    
#def draw_graphs(exon_graph_d):
#    
#    for geneid in exon_graph_d:
#        nx.draw_circular(exon_graph_d[geneid])
#        plt.savefig(geneid +'.png')
#        plt.close()
        
def parse_tpdm_line(tpdmline):
    
    out = {}
    ts = tpdmline.strip('\n').split('\t')
    
    out['read_id'] = ts[0]
    
    out['gene1'] = ts[3]
    out['gene2'] = ts[4]
    out['transcript1'] = ts[1].split('|')[0]
    out['chrom1'] = ts[1].split('|')[1]
    out['tstart1'] = int(ts[1].split('|')[2])
    out['tend1'] = int(ts[1].split('|')[3])
    out['transcript2'] = ts[2].split('|')[0]
    out['chrom2'] = ts[2].split('|')[1]
    out['tstart2'] = int(ts[2].split('|')[2])
    out['tend2'] = int(ts[2].split('|')[3])
    out['pos1'] = int(ts[5])
    out['pos2'] = int(ts[6])
    out['length1'] = int(ts[7])
    out['length2'] = int(ts[8])
    out['mismatches1'] = [int(k) for k in ts[9].split(',')[:-1]]
                                            ##There's an extra ',' at the end
    out['mismatches2'] = [int(k) for k in ts[10].split(',')[:-1]]
    
    return out
    
    
def compare_transcripts(ttuple1, ttuple2):
    return ttuple1[0] == ttuple2[0] and abs(ttuple1[2] - ttuple2[2]) < 125 and\
           abs(ttuple1[3] - ttuple2[3]) < 125 and ttuple1[1] == ttuple2[1]
    
    
def parse_line(efline):
    
    out = {}
    
    lseq = efline.strip().split('\t')
    out['transcriptid'] = lseq[0]
    out['geneid'] = lseq[1]
    out['exonid'] = lseq[4]
    out['chrom'] = '_'.join(lseq[4].split('_')[:-2])
    out['start'] = int(lseq[5])
    out['end'] = int(lseq[6])
    out['transcriptstart'] = int(lseq[2])
    out['transcriptend'] = int(lseq[3])
    
    out['length'] = abs(out['start'] - out['end']) ##Add + 1 for Ensembl data!

    exon_position = int(lseq[7])
    
    if exon_position == 1:
        out['length'] += PREPENDED_SEQUENCE_LENGTH
    return out
    
class exon_graph_lookup(object):
    
    def __init__(self):
        self._transcriptid_dict = {}
        
    def add_graph(self, transcriptid, chrom, start, end, egraph, ebounds):
        
        self._transcriptid_dict.setdefault(transcriptid, []).append(
                    (chrom, start, end, egraph, ebounds)
                    )
                    
    def __distance(self, start1, start2, end1, end2, chrom1, chrom2):
        pad = 0
        if chrom1 != chrom2:
            pad += 1000000000
        return abs(start1 - start2) + abs(end1 - end2) + pad
        
    def get_graph(self, transcriptid, chrom, start, end):
        
        egraphs = self._transcriptid_dict[transcriptid]
        
        if len(egraphs) == 1:
            return egraphs[0][3:]
            
        egraphs.sort(key=lambda k:self.__distance(start,k[1],end,k[2],chrom,k[0]))
        
        print "Matching", transcriptid, chrom, start, end, "with", egraphs[0][:3] 
        return egraphs[0][3:]



def main(exonfilename, tpdmfilename, concordantfilename):
    
    
    concordantset = set()
    for line in file(concordantfilename):
        concordantset.add(line.strip())

    
#==============================================================================
#     First, build the exon_graph.
#==============================================================================
    exon_graph_d = exon_graph_lookup()
    exon_graph = nx.DiGraph()
    
#    current_geneid = ""
    current_transcriptid = ""
    current_chrom = ""
    current_tstart = None
    current_tend = None
#    current_exonid = ""    
    
    current_exon_boundaries = None
    
    current_node = ""
    
    position_in_transcript = 0
    
    for line in file(exonfilename):
        exon_d = parse_line(line)
                                    
        if (exon_d['transcriptid'], exon_d['chrom'], exon_d['transcriptstart'], exon_d['transcriptend']) != \
            (current_transcriptid, current_chrom, current_tstart, current_tend):
            
            if current_transcriptid:
                exon_graph.node[current_node]['length'] += POLY_A_LENGTH
                exon_graph.add_edge(current_node, current_transcriptid + "sink")
                exon_graph_d.add_graph(current_transcriptid, current_chrom, 
                                       current_tstart, current_tend, exon_graph,
                                       current_exon_boundaries)
                
            exon_graph = nx.DiGraph() ##Make a new exon_graph
            exon_graph.add_node(exon_d['transcriptid'] + "source", length = 0) ##Add its source and sink
            exon_graph.add_node(exon_d['transcriptid'] + "sink", length = 0)
            
            current_node = exon_d['transcriptid'] + "source"
            current_transcriptid = exon_d['transcriptid']
            current_chrom = exon_d["chrom"]
            current_tstart = exon_d['transcriptstart']
            current_tend = exon_d['transcriptend']
            position_in_transcript = 0
            
            current_exon_boundaries = {}
        
        exon_graph.add_node(exon_d['exonid'], length = exon_d['length'])

        exon_graph.add_edge(current_node, exon_d['exonid'])
        current_node = exon_d['exonid']
        
        current_exon_boundaries[exon_d['exonid']] = \
            (position_in_transcript, position_in_transcript + exon_d['length'])
        position_in_transcript += exon_d['length']
        
    exon_graph.node[current_node]['length'] += POLY_A_LENGTH
    exon_graph.add_edge(current_node, current_transcriptid + "sink")
    exon_graph_d.add_graph(current_transcriptid, current_chrom, 
            current_tstart, current_tend, exon_graph, current_exon_boundaries)
    
    print "Exon structures loaded."
#===============================================================================
#     Now, iterate through the tdpm file and mark mapped exons.
#===============================================================================
    
    current_transcript_pair = (None, None)
    current_chrom_pair = (None, None)
    current_tstart_pair = (None, None)
    current_tend_pair = (None, None)
#    current_gene_pair = (None, None)
    
    exon_graph1 = None
    exon_graph2 = None
    exon_bounds1 = None
    exon_bounds2 = None
    
    mapped_reads_1 = []
    mapped_reads_2 = []
    
    outef = open(tpdmfilename + '.exons', 'w')
    read_set = set()
    linecount = 0
    concordantcount = 0
    
    for line in file(tpdmfilename):
        
        linecount += 1
        if linecount % 500000 == 0: 
            print linecount, concordantcount
        
        tpdm_d = parse_tpdm_line(line)

        if tpdm_d['read_id'] in concordantset:
            concordantcount += 1
            continue
        
        if (tpdm_d['transcript1'],tpdm_d['transcript2']) != current_transcript_pair or\
           (tpdm_d['chrom1'],tpdm_d['chrom2']) != current_chrom_pair or\
           (tpdm_d['tstart1'],tpdm_d['tstart2']) != current_tstart_pair or\
           (tpdm_d['tend1'],tpdm_d['tend2']) != current_tend_pair:
            
            if current_transcript_pair[0] and len(read_set) >= MIN_MAPPING_COUNT:
                write_exon_orders(exon_graph1, exon_graph2, 
                    current_transcript_pair[0],current_transcript_pair[1], 
                    outef, mapped_reads_1, mapped_reads_2, exon_bounds1,
                    exon_bounds2)
                
            current_transcript_pair = (tpdm_d['transcript1'],tpdm_d['transcript2'])
            current_chrom_pair = (tpdm_d['chrom1'],tpdm_d['chrom2'])
            current_tstart_pair = (tpdm_d['tstart1'],tpdm_d['tstart2'])
            current_tend_pair = (tpdm_d['tend1'],tpdm_d['tend2'])
            
            mapped_reads_1 = []
            mapped_reads_2 = []
            
            exon_graph1, exon_bounds1 = exon_graph_d.get_graph(current_transcript_pair[0],
                                                 current_chrom_pair[0],
                                                 current_tstart_pair[0],
                                                 current_tend_pair[0])
            exon_graph2, exon_bounds2 = exon_graph_d.get_graph(current_transcript_pair[1],
                                                 current_chrom_pair[1],
                                                 current_tstart_pair[1],
                                                 current_tend_pair[1])
                    
            read_set = set()
#===============================================================================
#         Recall that t1_bounds looks like [0, 123, 323, 542] where each number is
#         the start site of the kth exon.
#===============================================================================
              
        read_set.add((tpdm_d['pos1'], tpdm_d['pos2']))
        
        
        read1_start = good_read_start(tpdm_d['pos1'], tpdm_d['length1'], 
                                                        tpdm_d['mismatches1'])
        read1_end = tpdm_d['pos1'] + good_read_length(tpdm_d['length1'], 
                                                      tpdm_d['mismatches1'])
        mapped_reads_1.append( (tpdm_d['read_id'], read1_start, read1_end) )
        
        read2_start = good_read_start(tpdm_d['pos2'], tpdm_d['length2'], 
                                                        tpdm_d['mismatches2'])
        read2_end = tpdm_d['pos2'] + good_read_length(tpdm_d['length2'], 
                                                      tpdm_d['mismatches2'])
        mapped_reads_2.append( (tpdm_d['read_id'], read2_start, read2_end) )
        

        
        

    if len(read_set) > MIN_MAPPING_COUNT:
        write_exon_orders(exon_graph1, exon_graph2, current_transcript_pair[0], 
                current_transcript_pair[1], outef, mapped_reads_1, 
                mapped_reads_2, exon_bounds1, exon_bounds2)   
        
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
