#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
from socket import CAN_BCM_CAN_FD_FRAME
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file)as f:
        lines = f.readlines()
        for i in range(1,len(lines),4):
            yield lines[i].strip()



def cut_kmer(read, kmer_size):
    flag_start = 0
    flag_end = kmer_size
    while flag_end <= len(read):
        yield(read[flag_start:flag_end])
        flag_start = flag_start+1
        flag_end = flag_end + 1


def build_kmer_dict(fastq_file, kmer_size):
    memoire = {}
    file = read_fastq(fastq_file)
    for read in file:
        kmer = cut_kmer(read, kmer_size)
        for k in kmer:
            if k not in memoire.keys():
                memoire[k] = 1
            else:
                memoire[k] += 1
    return memoire



def build_graph(kmer_dict):
    g = nx.DiGraph()
    for k,i in kmer_dict.items():
        g.add_edge(k[:-1], k[1:], weight = i)
    return g

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    entre, sortie = 1, None
    if delete_entry_node == True:
        entre = 0
    if delete_sink_node == False:
        sortie = -1
    for chemin in path_list:
        graph.remove_nodes_from(chemin[entre:sortie]) 
    return graph


def std(data):
    return statistics.stdev(data)   


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    max_weight = max(mean_weights)
    heaviest = [i for i, j in enumerate(mean_weights) if j == max_weight]
    if len(heaviest) > 1:
        max_len = max(path_lens)
        longest = [i for i in heaviest if path_lens[i] == max_len]
        if len(longest) > 1:
            Random.seed(9001)
            best = random.choice[longest]
        else:
            best = longest[0]
    else:
        best = heaviest[0]
    
    for p in paths:
        print(p)

    paths2 = [p for p in paths]
    paths2.pop(best)

    return remove_paths(graph, paths2, delete_entry_node, delete_sink_node)




def path_average_weight(graph, path):
    route = graph.subgraph(path).edges(data=True)
    poids = 0

    for chemin in route:
        poids += chemin[2]["weight"]

    poids = poids/(len(path)-1)
    return poids


def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    memoire= []
    for node in graph.nodes():
        if(len(list(graph.predecessors(node)))== 0):
            memoire.append(node)
    return(memoire)

def get_sink_nodes(graph):  
    memoire= []
    for node in graph.nodes():
        if(len(list(graph.successors(node)))== 0):
            memoire.append(node)
    return(memoire)

def get_contigs(graph, starting_nodes, ending_nodes):
    memoire = []
    for entre in starting_nodes:
        for sortie in ending_nodes:
            if(nx.has_path(graph, entre, sortie)):
                az = nx.all_simple_paths(graph, entre, sortie)
                for t in az:
                    contig = "".join(t[::len(t[0])])
                    memoire.append((contig, len(contig)))
    return memoire

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as f:
        for i,typl in enumerate(contigs_list):
            f.write(f">contig_{i} len={typl[1]}\n")
            f.write(fill(typl[0]))
            f.write("\n")



def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))
    
def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    g = build_graph(build_kmer_dict(args.fastq_file,args.kmer_size))
    start = get_starting_nodes(g)
    end = get_sink_nodes(g)
    contig = get_contigs(g, start, end)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
