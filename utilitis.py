import streamlit.components.v1 as components
import networkx as nx
from pyvis.network import Network

def create_ggi(dpi_select, G,selected_drugs):
    # IF selected drugs>2: Create networkx graph object for gene-gene visualisation
    # Initiate PyVis network object
    got_net_ggi = Network(
                    height='750px',
                    width='100%',
                    bgcolor='#222222',
                    font_color='white'
                    )
    got_net_ggi.barnes_hut()
    got_net_ggi.repulsion(
                node_distance=270,
                central_gravity=0.33,
                spring_length=80,
                spring_strength=0.10,
                damping=0.95
                )


    # filter genes in two-hop
    # or construct a dataframe with 4 columns from lists
    paths = nx.all_simple_paths(G, source=selected_drugs[0], target=selected_drugs[1],cutoff=3)
    ps = [p for p in paths]
    #有时候没有two-hop.
    if len(ps) < 1:
        paths = nx.all_shortest_paths(G, source=selected_drugs[0], target=selected_drugs[1])
        ps = [p for p in paths]

    #从大到小sort list长度. 这样node的颜色不会丢失
    ps = sorted(ps, key=len, reverse=False)
    
    for ls in ps:
        got_net_ggi.add_node(ls[0], title=ls[0], color='#6BAEA9', shape='triangle',labelHighlightBold=True)
        got_net_ggi.add_node(ls[-1], title=ls[-1], color='#6BAEA9', shape='triangle',labelHighlightBold=True)

        rest = ls[1:-1]

        got_net_ggi.add_edge(ls[0], ls[-1], color='#CDCDCD')
        common_nodes = list()
        if len(rest) == 1:
            got_net_ggi.add_node(rest[0], title=rest[0], color='#0087bd', shape='star')
            got_net_ggi.add_edge(ls[0], rest[0], color='#CDCDCD')
            got_net_ggi.add_edge(ls[-1], rest[0], color='#CDCDCD')
            common_nodes.append(rest[0])
        elif len(rest) == 2:
            if rest[0] in common_nodes:
                #添加另一个
                got_net_ggi.add_node(rest[1], title=rest[1], color='#A0AA9B', shape='star')
            elif rest[1] in common_nodes:
                got_net_ggi.add_node(rest[0], title=rest[0], color='#A0AA9B', shape='star')

            else:
                got_net_ggi.add_node(rest[1], title=rest[1], color='#A0AA9B', shape='star')
                got_net_ggi.add_node(rest[0], title=rest[0], color='#A0AA9B', shape='star')

            got_net_ggi.add_edge(ls[0], rest[0], color='#CDCDCD')
            got_net_ggi.add_edge(ls[-1], rest[1], color='#CDCDCD')
            got_net_ggi.add_edge(rest[0], rest[1], color='#CDCDCD')
        
        elif len(rest) == 3:
            got_net_ggi.add_node(rest[0], title=rest[0], color='#A0AA9B', shape='star')
            got_net_ggi.add_node(rest[1], title=rest[1], color='#A0AA9B', shape='star')
            got_net_ggi.add_node(rest[2], title=rest[2], color='#A0AA9B', shape='star')

            got_net_ggi.add_edge(ls[0], rest[0], color='#CDCDCD')
            got_net_ggi.add_edge(ls[-1], rest[1], color='#CDCDCD')
            got_net_ggi.add_edge(rest[0], rest[1], color='#CDCDCD')
            got_net_ggi.add_edge(rest[1], rest[2], color='#CDCDCD')

    #0087bd

        #drug本身连接到的基因也要挤上
    dpi_edge_data = zip(dpi_select['source'], dpi_select['target'])

    for src, dst in dpi_edge_data: #gene_name, 
        #add nodes and edges to the graph
        #drug
        # got_net_ggi.add_node(src, title=src, color='#6BAEA9', shape='triangle',labelHighlightBold=True)
        #protein
        if dst not in common_nodes:
            got_net_ggi.add_node(dst, title=dst, color='#A0AA9B', shape='star')

        got_net_ggi.add_edge(src, dst, color='#CDCDCD')


    # Save and read graph as HTML file (on Streamlit Sharing)
    try:
        path = './tmp'
        got_net_ggi.save_graph(f'{path}/pyvis_graph_ggi.html')
        HtmlFile = open(f'{path}/pyvis_graph_ggi.html', 'r', encoding='utf-8')

    # Save and read graph as HTML file (locally)
    except:
        path = './html_files'
        got_net_ggi.save_graph(f'{path}/pyvis_graph_ggi.html')
        HtmlFile = open(f'{path}/pyvis_graph_ggi.html', 'r', encoding='utf-8')

    # Load HTML file in HTML component for display on Streamlit page
    components.html(HtmlFile.read(), height=435)

def add_drug_side(df_select,got_net):
     # create graph using pviz network 
    edge_data = zip(df_select['source'], df_select['target'], df_select['rel'])

    for src, dst, rel in edge_data:
        #add nodes and edges to the graph
        #drug
        got_net.add_node(src, src, title=src, color='#6BAEA9', shape='triangle',labelHighlightBold=True)
        #side effect
        got_net.add_node(dst, dst, title=dst, color='#F08327', shape='dot')

        got_net.add_edge(src, dst, color='#CDCDCD')

    
def add_drug_drugddi_select(ddi_select,got_net):
    # ddi
    ddi_edge_data = zip(ddi_select['DRUG_1_CONCEPT_NAME'], ddi_select['DRUG_2_CONCEPT_NAME'], ddi_select['EVENT_CONCEPT_NAME'],\
                        ddi_select['MICROMEDEX_SEV_LEVEL'])

    for src_1, src_2, rel, servere in ddi_edge_data:
        #drug
        got_net.add_node(src_1, src_1, title=src_1, color='#6BAEA9', shape='triangle',labelHighlightBold=True)
        #drug
        got_net.add_node(src_2, src_2, title=src_2, color='#6BAEA9', shape='triangle',labelHighlightBold=True)
        #side effect
        got_net.add_node(rel, rel, title=rel, color='#CDCDCD', shape='dot')
        
        #{'Contraindicated', 'Major', 'Minor', 'Moderate'}
        if servere == 'Major':
            val = 1.4
            clr = '#746AB0'
        elif servere == 'Moderate':
            val = 1.4
            clr = '#D7AA73'
        elif servere == 'Minor':
            val = 1
            clr = '#288BA8'
        else:
            val = 1
            clr = '#288BA8'
        got_net.add_edge(src_1, src_2, color=clr,value = val)

        got_net.add_edge(src_1, rel, color='#CDCDCD')
        got_net.add_edge(src_2, rel, color='#CDCDCD')