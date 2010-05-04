import networkx as nx
#import matplotlib.pyplot as plt
import sys

MAX_MARK_LENGTH = 750
READ_LENGTH = 40
MIN_MAPPING_COUNT = 4
MIN_GOOD_OVERLAP = 5


def dfs(v, g, sinkname, marked1_set, marked2_set, outfileh, firstgene, 
        secondgene, length_1 = 0, exon_list = [], 
        seenmarked1 = set(), seenmarked2 = set(), marked1 = False, 
        marked2 = False, marked_distance = 0):
            
    fail = False
#    print '*'*80
#    print "Visiting node", v
#===============================================================================
#     Check if we've passed marked nodes and if so, how long between marked1 and
#     marked2.
#===============================================================================
    if (not marked1) and g.node[v]['marked'] == True and g.node[v]['graph'] == 1:
        ismarked1 = True
    else:
        ismarked1 = marked1
    
    if (not marked2) and g.node[v]['marked'] == True and g.node[v]['graph'] == 2:
        ismarked2 = True
    else:
        ismarked2 = marked2
        
    if marked1 and not ismarked2:
        m_d = marked_distance + g.node[v]['length']
    else:
        m_d = int(marked_distance)
    
    if m_d > MAX_MARK_LENGTH:
#        print "Failing because we've gone too far after a marked node."
        fail = True
    
#    print "graph", g.node[v]['graph']
#    print "marked?", g.node[v]['marked']
#    print "ismarked1?", ismarked1, "ismarked2?", ismarked2
#    print "m_d", m_d
#    print "exon_list", exon_list
#    print "seenmarked1", seenmarked1
#    print "seenmarked2", seenmarked2
#===============================================================================
#     See if we've entered into the second graph without hitting all the marked
#     nodes in the first.
#===============================================================================
    if g.node[v]['graph'] == 2 and seenmarked1 != marked1_set:
#        print "Failing because in graph 2 but haven't seen all marked."
        fail = True

    
#===============================================================================
#     See if we're at the sink of the graph and check if we've found something
#     good.
#===============================================================================
    
    if v == sinkname:
        if seenmarked1 == marked1_set and seenmarked2 == marked2_set:
#            print "WRITING"
            outfileh.write(firstgene + '.' + secondgene + '\t' + str(length_1) + '\t')
            outfileh.write('\t'.join(exon_list[1:]))
            outfileh.write('\n')
            
        fail = True
    
    
    
    if not fail:
        
        if g.node[v]['graph'] == 1:
            l_1 = int(length_1) + g.node[v]['length']
        elif g.node[v]['graph'] == 2:
            l_1 = int(length_1)
        
        if g.node[v]['marked'] == True and g.node[v]['graph'] == 1:
            seenmarked1.add(v)
        if g.node[v]['marked'] == True and g.node[v]['graph'] == 2:
            seenmarked2.add(v)
        e_l = list(exon_list)
        e_l.append(v)
        for neighbor in g.neighbors(v):
            dfs(neighbor, g, sinkname, marked1_set, marked2_set, outfileh,
                firstgene, secondgene, l_1, e_l, set(seenmarked1), 
                set(seenmarked2), ismarked1, ismarked2, m_d)


def write_exon_orders(eg1, eg2, gene1, gene2, outfileh, marked1, marked2):
        
    g1 = nx.union(eg1, eg2)
    
    for node1 in eg1.nodes():
        g1.node[node1] = eg1.node[node1]
    for node2 in eg2.nodes():
        g1.node[node2] = eg2.node[node2]   
        
    for node1 in eg1.nodes():
        g1.node[node1]['graph'] = 1
        
    for node2 in eg2.nodes():
        g1.node[node2]['graph'] = 2 
        
    for node1 in eg1.nodes():
        for node2 in eg2.nodes():
            if node1 not in (gene1 + "source", gene1 + "sink") and \
               node2 not in (gene2 + "source", gene2 + "sink"):
                g1.add_edge(node1, node2)
    
#    print gene1 + "source", gene2 + "sink"
#    print "marked1", marked1
#    print "marked2", marked2
    
    dfs(gene1 + "source", g1, gene2 + "sink", marked1, marked2, outfileh, gene1, gene2)
    
    ##############################################################
                
#    g2 = nx.union(eg1, eg2)
#    
#    for node1 in eg1.nodes():
#        g2.node[node1] = eg1.node[node1]
#    for node2 in eg2.nodes():
#        g2.node[node2] = eg2.node[node2]
#        
#    for node1 in eg1.nodes():
#        g2.node[node1]['graph'] = 2
#        
#    for node2 in eg2.nodes():
#        g2.node[node2]['graph'] = 1    
#    
#    
#    for node1 in eg1.nodes():
#        for node2 in eg2.nodes():
#            if node1 not in (gene1 + "source", gene1 + "sink") and \
#               node2 not in (gene2 + "source", gene2 + "sink"):
#                g2.add_edge(node2, node1)
#   
#        
#    dfs(gene2 + "source", g2, gene1 + "sink", marked2, marked1, outfileh, gene2, gene1)    
    
def draw_graphs(exon_graph_d):
    
    for geneid in exon_graph_d:
        nx.draw_circular(exon_graph_d[geneid])
        plt.savefig(geneid +'.png')
        plt.close()
        
def parse_tpdm_line(tpdmline):
    
    out = {}
    ts = tpdmline.strip('\n').split('\t')
    
    out['read_id'] = ts[0]
    
    out['gene1'] = ts[3]
    out['gene2'] = ts[4]
    out['transcript1'] = ts[1]
    out['transcript2'] = ts[2]
    out['pos1'] = int(ts[5])
    out['pos2'] = int(ts[6])
    out['mismatches1'] = [int(k) for k in ts[7].split(',')[:-1]] ##There's an extra ',' at the end
    out['mismatches2'] = [int(k) for k in ts[8].split(',')[:-1]]
    
    return out

def clear_exon_graph(eg):
    
    for node in eg.nodes():
        eg.node[node]['marked'] = False
    
def parse_line(efline):
    
    out = {}
    
    lseq = efline.strip().split('\t')
    out['transcriptid'] = lseq[0]
    out['geneid'] = lseq[1]
    out['exonid'] = lseq[4]
    out['start'] = int(lseq[5])
    out['end'] = int(lseq[6])
    out['length'] = abs(out['start'] - out['end']) + 1

    return out

def main(exonfilename, tpdmfilename):
    
    
#===============================================================================
#     First, build the exon_graph.
#===============================================================================
    exon_graph_d = {}
    exon_graph = nx.DiGraph()
    
    current_geneid = ""
    current_transcriptid = ""
    current_exonid = ""    
    
    current_node = ""
    
    position_in_transcript = 0
    transcript_to_exon = {}
    
    for line in file(exonfilename):
        exon_d = parse_line(line)
        
        if exon_d['geneid'] != current_geneid:  ##If new gene
            current_geneid = exon_d['geneid']
                            
        if exon_d['transcriptid'] != current_transcriptid:
            
            if current_transcriptid:
                exon_graph.add_edge(current_node, current_transcriptid + "sink")
                exon_graph_d[current_transcriptid] = exon_graph
                
            exon_graph = nx.DiGraph() ##Make a new exon_graph
            exon_graph.add_node(exon_d['transcriptid'] + "source", marked = False, length = 0) ##Add its source and sink
            exon_graph.add_node(exon_d['transcriptid'] + "sink", marked = False, length = 0)
            
            current_node = exon_d['transcriptid'] + "source"
            current_transcriptid = exon_d['transcriptid']
            position_in_transcript = 0
            transcript_to_exon[current_transcriptid] = [ [], [] ]
        
        exon_graph.add_node(exon_d['exonid'], marked = False, length = exon_d['length'])

        exon_graph.add_edge(current_node, exon_d['exonid'])
        current_node = exon_d['exonid']
        
        transcript_to_exon[current_transcriptid][0].append(position_in_transcript)
        transcript_to_exon[current_transcriptid][1].append(exon_d['exonid'])
        position_in_transcript += exon_d['length']
        
        
    exon_graph.add_edge(current_node, current_transcriptid + "sink")
    exon_graph_d[current_transcriptid] = exon_graph
    
    print "Exon structures loaded."
#===============================================================================
#     Now, iterate through the tdpm file and mark mapped exons.
#===============================================================================
    
    current_transcript_pair = (None, None)
    current_gene_pair = (None, None)
    
    exon_graph1 = None
    exon_graph2 = None
    
    marked1 = set()
    marked2 = set()
    
    outef = open(tpdmfilename + '.exons', 'w')
    read_set = set()
    linecount = 0
    
    for line in file(tpdmfilename):
        
        linecount += 1
        if linecount % 500000 == 0: 
            print linecount
        
        tpdm_d = parse_tpdm_line(line)
        
        if (tpdm_d['transcript1'],tpdm_d['transcript2']) != current_transcript_pair:
            
            if current_transcript_pair[0] and len(read_set) > MIN_MAPPING_COUNT:
                write_exon_orders(exon_graph1, exon_graph2, current_transcript_pair[0], 
                                  current_transcript_pair[1], outef, marked1, marked2)
                
            current_gene_pair = (tpdm_d['gene1'],tpdm_d['gene2'])
            current_transcript_pair = (tpdm_d['transcript1'],tpdm_d['transcript2'])
            
            exon_graph1 = exon_graph_d[current_transcript_pair[0]]
            exon_graph2 = exon_graph_d[current_transcript_pair[1]]

            clear_exon_graph(exon_graph1)
            clear_exon_graph(exon_graph2)
                    
            marked1 = set()
            marked2 = set()
            
            read_set = set()
#===============================================================================
#         Recall that t1_exons looks like [0, 123, 323, 542] where each number is
#         the start site of the kth exon.
#===============================================================================
              
        t1_bounds, t1_exons = transcript_to_exon[tpdm_d['transcript1']]
        assert len(t1_bounds) == len(t1_exons)

        read_set.add(tpdm_d['read_id'])
        
        read_bound1_start = tpdm_d['pos1']
        read_bound1_end = tpdm_d['pos1'] + READ_LENGTH - 1

        for i in range(len(t1_exons)):
            exon_start = t1_bounds[i]
            try:
                exon_end = t1_bounds[i+1] - 1
            except IndexError:
                exon_end = 999999
            
            if (exon_start <= read_bound1_start < exon_end) or \
               (exon_start <= read_bound1_end < exon_end) or \
               (read_bound1_start <= exon_start and exon_end < read_bound1_end):
                   
#                exon_positions = range(max(exon_start,read_bound1_start) , min(exon_end,read_bound1_end) + 1)
#                read_positions = range(read_bound1_start, read_bound1_end + 1)
#                overlap = list(set(read_positions).intersection(exon_positions))
                overlap = range(max(exon_start,read_bound1_start) , min(exon_end,read_bound1_end) + 1)
                if overlap:
                    overlap_read_indices = [x - read_bound1_start for x in overlap]
                    good_overlap = [x for x in overlap_read_indices \
                                    if x not in tpdm_d['mismatches1']]
                                    
                    if len(good_overlap) >= MIN_GOOD_OVERLAP:
                        exon_graph1.node[t1_exons[i]]['marked'] = True
                        marked1.add(t1_exons[i])
        
        
        t2_bounds, t2_exons = transcript_to_exon[tpdm_d['transcript2']]
        assert len(t2_bounds) == len(t2_exons)

        read_bound2_start = tpdm_d['pos2']
        read_bound2_end = tpdm_d['pos2'] + READ_LENGTH - 1
        
        for i in range(len(t2_exons)):
            exon_start = t2_bounds[i]
            try:
                exon_end = t2_bounds[i+1] - 1
            except IndexError:
                exon_end = 999999
            
            if (exon_start <= read_bound2_start < exon_end) or \
               (exon_start <= read_bound2_end < exon_end) or \
               (read_bound2_start <= exon_start and exon_end < read_bound2_end):
                   
#                exon_positions = range(max(exon_start,read_bound2_start) , min(exon_end,read_bound2_end) + 1)
#                read_positions = range(read_bound2_start, read_bound2_end + 1)
                overlap = range(max(exon_start,read_bound2_start) , min(exon_end,read_bound2_end) + 1)
                
                if overlap:
                    overlap_read_indices = [x - read_bound2_start for x in overlap]
                    good_overlap = [x for x in overlap_read_indices \
                                    if x not in tpdm_d['mismatches2']]
                    if len(good_overlap) >= MIN_GOOD_OVERLAP:
                        exon_graph2.node[t2_exons[i]]['marked'] = True
                        marked2.add(t2_exons[i])

    write_exon_orders(exon_graph1, exon_graph2, current_transcript_pair[0], 
                                current_transcript_pair[1], outef, marked1, marked2)    
        
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
