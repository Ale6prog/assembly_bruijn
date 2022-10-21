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


import pickle
import argparse
import os
import random
from random import randint
import statistics
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
random.seed(9001)


__author__ = "Bel Alexis"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Bel Alexis"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Bel Alexis"
__email__ = "alexbel28@yahoo.fr"
__status__ = "Developpement"


def isfile(path_file):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path_file):
        if os.path.isdir(path_file):
            msg = f"{path_file} is a directory"
        else:
            msg = f"{path_file} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path_file


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="-h")
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
    """Read file fasta
      :Parameters:
        fastq_file: fasta
    """
    with open(fastq_file, encoding="UTF-8")as file:
        lines = file.readlines()
        for i in range(1, len(lines), 4):
            yield lines[i].strip()


def cut_kmer(read, kmer_size):
    """Cut a sequence in k-mers.
    """
    flag_start = 0
    flag_end = kmer_size
    while flag_end <= len(read):
        yield read[flag_start:flag_end]
        flag_start = flag_start+1
        flag_end = flag_end + 1


def build_kmer_dict(fastq_file, kmer_size):
    """ Give a dictionnary of k-mers.
    """
    memoire = {}
    for read in read_fastq(fastq_file):
        for k in cut_kmer(read, kmer_size):
            if k not in memoire:
                memoire[k] = 1
            else:
                memoire[k] += 1
    return memoire


def build_graph(kmer_dict):
    """ Build a graph.
    """
    graph = nx.DiGraph()
    for kmer, compte in kmer_dict.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=compte)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Remove a path from the graph.
    """
    entre, sortie = 1, None
    if delete_entry_node:
        entre = 0
    if not delete_sink_node:
        sortie = -1
    for chemin in path_list:
        graph.remove_nodes_from(chemin[entre:sortie])
    return graph


def std(data):
    """Return the standard deviation of a data
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """ Select the best path on the graph by using a path list.
    """
    std_weight = std(weight_avg_list)
    std_lenght = std(path_length)

    if std_weight > 0:
        best_path = path_list[weight_avg_list.index(max(weight_avg_list))]
    else:
        if std_lenght > 0:
            best_path = path_list[path_length.index(max(path_length))]
        else:
            best_path = path_list[randint(0, len(path_list)-1)]
    path_list.remove(best_path)
    return remove_paths(graph, path_list, delete_entry_node, delete_sink_node)


def path_average_weight(graph, pathway):
    """ Give the mean of weight of a path.
    """
    route = graph.subgraph(pathway).edges(data=True)
    poids = 0

    for chemin in route:
        poids += chemin[2]["weight"]

    poids = poids/(len(pathway)-1)
    return poids


def solve_bubble(graph, ancestor_node, descendant_node):
    """ Solve a bubble from the graph
    """
    p_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    graph = select_best_path(graph, p_list, [len(path) for path in p_list],
                             [(path_average_weight(graph, path)) for path
                             in p_list])
    return graph


def simplify_bubbles(graph):
    """ Solve all the bubble from the graph.
    """
    flag = False

    for node in graph:
        node_n = node

        if graph.in_degree(node) >= 2:
            node_pred = []

            for node in graph.predecessors(node):
                node_pred += [node]

            if len(node_pred) <= 2:
                flag = True
                break

        if flag:
            break

    if not flag:
        ancestor = nx.lowest_common_ancestor(graph, node_pred[0], node_pred[1])
        graph = solve_bubble(graph, ancestor, node_n)

    return graph


def solve_entry_tips(graph, starting_nodes):
    """ Solve all the entry tips of the graph.
    """
    n_list = []

    if len(starting_nodes) > 2:
        n_list = [(starting_nodes[n_i], starting_nodes[n_j]) for n_i in range(len(starting_nodes)) for n_j
                  in range(n_i+1, len(starting_nodes))]
    elif len(starting_nodes) == 2:
        n_list = [tuple(starting_nodes)]
    else:
        return  graph

    a_list = []

    for node in n_list:
        a_list.append([nx.lowest_common_ancestor(graph.reverse(),
                       node[0], node[1])])

    p_list = []
    w_list = []
    p_len = []

    for i, node in enumerate(n_list):
        path_1 = list(nx.all_simple_paths(graph, node[0], a_list[i]))[0]
        path_2 = list(nx.all_simple_paths(graph, node[1], a_list[i]))[0]

        p_list += [path_1, path_2]

        p_len += [len(path_1), len(path_2)]
        w_list += [path_average_weight(graph, path_1),
                   path_average_weight(graph, path_2)]

    graph = select_best_path(graph, p_list, p_len, w_list,
                             delete_entry_node=True)

    return graph


def solve_out_tips(graph, ending_nodes):
    """ Solve all the exit tips of the graph.
    """
    n_list = []

    if len(ending_nodes) > 2:
        n_list = [(n_i, n_j) for n_i in ending_nodes for n_j in ending_nodes]
    elif len(ending_nodes) == 2:
        n_list = [tuple(ending_nodes)]
    else:
        return graph

    a_list = []

    for node in n_list:
        a_list += [nx.lowest_common_ancestor(graph, node[0], node[1])]

    p_list = []
    w_list = []
    p_len = []

    for i, node in enumerate(n_list):
        path_1 = list(nx.all_simple_paths(graph, a_list[i], node[0]))[0]
        path_2 = list(nx.all_simple_paths(graph, a_list[i], node[1]))[0]

        p_list += [path_1, path_2]

        p_len += [len(path_1), len(path_2)]
        w_list += [path_average_weight(graph, path_1),
                   path_average_weight(graph, path_2)]

    graph = select_best_path(graph, p_list, p_len, w_list,
                             delete_sink_node=True)

    return graph


def get_starting_nodes(graph):
    """ Give all the starting nodes from the graph.
    """
    memoire = [node for node in graph.nodes() if len(
               list(graph.predecessors(node))) == 0]
    return memoire


def get_sink_nodes(graph):
    """ Give all the contigs.
    """
    memoire = [node for node in graph.nodes() if len(
               list(graph.successors(node))) == 0]
    return memoire


def get_contigs(graph, starting_nodes, ending_nodes):
    """ Give all the starting nodes from the graph.
    """
    memoire = []
    for nodes in [(entre, sortie) for entre in starting_nodes for sortie
                  in ending_nodes if nx.has_path(graph, entre, sortie)]:
        for chemin in nx.all_simple_paths(graph, nodes[0], nodes[1]):
            contig = "".join(chemin[::len(chemin[0])])
            memoire.append((contig, len(contig)))
    return memoire


def save_contigs(contigs_list, output_file):
    """ Save the contigs in a file.
    """
    with open(output_file, "w", encoding="UTF-8") as file:
        for i, typl in enumerate(contigs_list):
            file.write(f">contig_{i+1} len={typl[1]}\n")
            file.write(f"{fill(typl[0])}\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight']
              > 3]
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight']
              <= 3]
    # Draw the graph with networkx
    pos = nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt", encoding="UTF-8") as save:
        pickle.dump(graph, save)


# Main program

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # k-mer
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    # Graph building
    graph = build_graph(kmer_dict)
    graph = simplify_bubbles(graph)
    start = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, start)
    end = get_sink_nodes(graph)
    graph = solve_out_tips(graph, end)
    end = get_sink_nodes(graph)

    # Contigs
    contig = get_contigs(graph, start, end)
    if args.output_file:
        save_contigs(contig, args.output_file)

    # Save image
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
