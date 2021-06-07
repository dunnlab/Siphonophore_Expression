#!/usr/bin/env python

def read_gaf(gaf_file):
    d = {}
    for line in open(gaf_file):
        args = line.split("\t")
        if args[0][0] != "!":
            gene_name = args[0]
            #print args
            go_id = args[2]
            try:
                d[gene_name].append(go_id)
            except KeyError:
                d[gene_name] = [go_id]
    return d


def main(gaf_file):
    gaf_dict = read_gaf(gaf_file)
    for gene in gaf_dict.keys():
        go_concat = ";GO:".join(gaf_dict[gene])
        anno = "{0}\tGO:{1}".format(gene,go_concat)
        print anno

main("Frillagalma_GO.out")
