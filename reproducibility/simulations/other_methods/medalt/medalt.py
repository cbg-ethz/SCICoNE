from optparse import OptionParser
from Edmonds import *
import os,sys
import subprocess
import numpy as np

#get the absolute path of input file
def getPath(path):
    path1=path.split("/")
    if path1[0] == ".":
        if (len(path1)==1):
            newpath=os.getcwd()
        else:
            newpath=os.getcwd()
            for ele in path1:
                if ele !=".":
                    newpath=newpath+"/"+ele
    elif path1[0]=="..":
        i = 0
        for ele in path1:
            if ele == "..":
                i=i+1
        path2=os.getcwd()
        path2=path2.split("/")
        newpath="/"+path2[0]
        if len(path2)-i > 1:
            for j in range(1,len(path2)-i):
                newpath=newpath+"/"+path2[j]
        for j in range(i,len(path1)):
            newpath=newpath+"/"+path1[j]
    else:
        newpath=path
    return newpath

# Returns a dictionary mapping node names to list of list of integers representing list of copy number list
def readl(l, n_bins):
    nodes = {}
    charlist=[]
    CNV={}
    line = l[0].split("\t")
    i = 0
    for line in l[1:]:
        array = line.split()
        snip = []
        CNVvalue = []
        snip.append(map(int, array[0:n_bins]))
        CNVvalue.extend(map(int, array[0:n_bins]))
        nodes[str(i)] = snip
        CNV[str(i)]=list(set(CNVvalue))
        i=i+1

    root = 'NA'
    for ele in CNV.keys():
        if CNV[ele] == [2]:
            root=ele
    if root == "NA":
        snip=[]
        snip.append([2]*(n_bins-0))
        nodes['root']=snip
        root='root'
    return nodes,root



def main():
    usage = "usage: python %prog <-P path> <-I input> <-D datatype>"
    description = "Input integer copy number profile. Columns correspond to chromosomal position. Rows correspond to cells."
    op = OptionParser(version="%prog 1.0",description=description,usage=usage,add_help_option=False)
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-I","--Input",dest="Input",type="str",
                  help="Input file")
    op.add_option("-O","--Output",dest="Output",type="str",
                  help="Output path")
                  
    (options,args) = op.parse_args()
    # check input parameters. Package path, input file, data type and genome version are required.
    #if not options.Path or not options.Input:
    #    op.print_help()
    #    sys.exit(1)

    # get the input parameters
    currentPath=os.getcwd()
    filename=options.Input
    filename=getPath(filename)
    outputpath = options.Output
    if not options.Output:
    	outpath=currentPath
    
    # Cell by bin matrix
    mat = np.loadtxt(filename, delimiter=' ')
    
    l = []
    l.append('\t'.join(np.arange(mat.shape[1]).astype(int).astype(str)))
    for i in range(0, mat.shape[0]):
        l.append('\t'.join(mat[i,:].astype(int).astype(str)))
    
    print("Inferring MEDALT.")

    #Identifying root node from input file.
    #If a diploidy genome is not input, will add an extra diploidy node as root
    (nodes,root) = readl(l,  mat[:,:].shape[1])
    print(root)
    node_name_list = nodes.keys()

    #calculation of MED distance
    g = create_tree(nodes, node_name_list,root)

    #Inference of tree and output
    result = compute_rdmst(g, root)
    write=open(outputpath,'w')
    tree=result[0]
    out1="from"+"\t"+"to"+"\t"+"dist"
    print >> write, out1
    for ele in tree.keys():
        out=ele
        for value in tree[ele].keys():
            out1=out+"\t"+value+"\t"+str(tree[ele][value])
            print >> write,out1
    write.close()
    print("MEDALT inferrence finish.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
