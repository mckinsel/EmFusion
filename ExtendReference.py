import networkx as nx
#import matplotlib.pyplot as plt
import sys

MAX_MARK_LENGTH = 500
READ_LENGTH = 40


def dfs(v, g, sinkname, marked1_set, marked2_set, outfileh, visited = None, 
        exon_list = [], seenmarked1 = set(), seenmarked2 = set(),
        marked1 = False, marked2 = False, marked_distance = 0):
            
    if visited is None: visited = set()
    
    
    fail = False
    
#    print "visiting node %s in graph %d" % (v, g.node[v]['graph'])
#    print "sinkname %s" % sinkname
#    print "marked1_set", marked1_set
#    print "marked2_set", marked2_set
#    print "exon_list", exon_list
    
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

#    print "marked1 %d ismarked1 %d" %(marked1, ismarked1)
#    print "marked2 %d ismarked2 %d" %(marked2, ismarked2)
#    print "m_d %d" % m_d
    
    if m_d > MAX_MARK_LENGTH:
#        print "FAILING due to m_d %d > MAX_MARK_LENGTH %d" % (m_d, MAX_MARK_LENGTH)
        fail = True
    
#===============================================================================
#     See if we've entered into the second graph without hitting all the marked
#     nodes in the first.
#===============================================================================
    if g.node[v]['graph'] == 2 and seenmarked1 != marked1_set:
#        print "FAILING due to getting to graph 2 before hitting marked1_set"
        fail = True

    
#===============================================================================
#     See if we're at the sink of the graph and check if we've found something
#     good.
#===============================================================================
    
    if v == sinkname:
#        print "reached sink!"
#        print exon_list
#        print seenmarked1
#        print seenmarked2
        if seenmarked1 == marked1_set and seenmarked2 == marked2_set:
            outfileh.write('\t'.join(exon_list[1:]))
            outfileh.write('\n')
            
        fail = True
    
    
    
    if not fail:
#        visited.add(v)
        if g.node[v]['marked'] == True and g.node[v]['graph'] == 1:
            seenmarked1.add(v)
        if g.node[v]['marked'] == True and g.node[v]['graph'] == 2:
            seenmarked2.add(v)
        e_l = list(exon_list)
        e_l.append(v)
        for neighbor in g.neighbors(v):
            if neighbor not in visited:
                dfs(neighbor, g, sinkname, marked1_set, marked2_set, outfileh,
                visited, e_l, set(seenmarked1), set(seenmarked2),
                ismarked1, ismarked2, m_d)


def write_exon_orders(eg1, eg2, gene1, gene2, outfileh, marked1, marked2):
    
#    print len(eg1), len(eg2)
#    print gene1, gene2
#    print set(eg1.nodes()).intersection(eg2.nodes())
#    for node in eg1.nodes():
#        print node, eg1.node[node]
#    for node in eg2.nodes():
#        print node, eg2.node[node]
#        
#    for edge in eg1.edges_iter():
#        print edge
#    for edge in eg2.edges_iter():
#        print edge
    
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
                

#    print 'eg1', eg1.node
#    print 'eg2', eg2.node
#    print 'g1', g1.node
    dfs(gene1 + "source", g1, gene2 + "sink", marked1, marked2, outfileh)
    
    ##############################################################
                
    g2 = nx.union(eg1, eg2)
    
    for node1 in eg1.nodes():
        g2.node[node1] = eg1.node[node1]
    for node2 in eg2.nodes():
        g2.node[node2] = eg2.node[node2]
        
    for node1 in eg1.nodes():
        g2.node[node1]['graph'] = 2
        
    for node2 in eg2.nodes():
        g2.node[node2]['graph'] = 1    
    
    
    for node1 in eg1.nodes():
        for node2 in eg2.nodes():
            if node1 not in (gene1 + "source", gene1 + "sink") and \
               node2 not in (gene2 + "source", gene2 + "sink"):
                g2.add_edge(node2, node1)
   
        
    dfs(gene2 + "source", g2, gene1 + "sink", marked2, marked1, outfileh)    
    
def draw_graphs(exon_graph_d):
    
    for geneid in exon_graph_d:
        nx.draw_circular(exon_graph_d[geneid])
        plt.savefig(geneid +'.png')
        plt.close()
        
def parse_tpdm_line(tpdmline):
    
    out = {}
    ts = tpdmline.strip().split('\t')
    
    out['gene1'] = ts[2]
    out['gene2'] = ts[3]
    out['transcript1'] = ts[0]
    out['transcript2'] = ts[1]
    out['pos1'] = int(ts[4])
    out['pos2'] = int(ts[5])
    
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
    out['length'] = abs(out['start'] - out['end'])

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
    
    for line in file(tpdmfilename):
        tpdm_d = parse_tpdm_line(line)
        
        if (tpdm_d['transcript1'],tpdm_d['transcript2']) != current_transcript_pair:
            
            if current_transcript_pair[0]:
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
#===============================================================================
#         Recall that t1_exons looks like [0, 123, 323, 542] where each number is
#         the start site of the kth exon.
#===============================================================================
        
        
        t1_bounds, t1_exons = transcript_to_exon[tpdm_d['transcript1']]
        assert len(t1_bounds) == len(t1_exons)

        read_bound1_start = tpdm_d['pos1']
        read_bound1_end = tpdm_d['pos1'] + READ_LENGTH
        
        for i in range(len(t1_exons)):
            exon_start = t1_bounds[i]
            try:
                exon_end = t1_bounds[i+1] - 1
            except IndexError:
                exon_end = 999999
        
            if (exon_start <= read_bound1_start < exon_end) or \
               (exon_start <= read_bound1_end < exon_end) or \
               (read_bound1_start <= exon_start and exon_end < read_bound1_end):
                   exon_graph1.node[t1_exons[i]]['marked'] = True
                   marked1.add(t1_exons[i])
        
        
        t2_bounds, t2_exons = transcript_to_exon[tpdm_d['transcript2']]
        assert len(t2_bounds) == len(t2_exons)

        read_bound2_start = tpdm_d['pos2']
        read_bound2_end = tpdm_d['pos2'] + READ_LENGTH
        
        for i in range(len(t2_exons)):
            exon_start = t2_bounds[i]
            try:
                exon_end = t2_bounds[i+1] - 1
            except IndexError:
                exon_end = 999999
        
            if (exon_start <= read_bound2_start < exon_end) or \
               (exon_start <= read_bound2_end < exon_end) or \
               (read_bound2_start <= exon_start and exon_end < read_bound2_end):
                   exon_graph2.node[t2_exons[i]]['marked'] = True
                   marked2.add(t2_exons[i])

    write_exon_orders(exon_graph1, exon_graph2, current_transcript_pair[0], 
                                current_transcript_pair[1], outef, marked1, marked2)    
        
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])