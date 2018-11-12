#!/usr/bin/env python
"""
A simple script to change key values in psfmap file for various versions of the
ngmixit program

"""
import string
import shutil
import getopt
import os
import sys
import glob


class Remap():
    
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile
       
    def remap_psf_map(self):
        
        if os.path.exists(self.outfile):
            os.remove(self.outfile)
        
        outF = open(self.outfile,'w')

        for line in open(self.infile,'r'):
            tokens = line.split(' ')
            line = 'D'+tokens[0]+'-'+tokens[1]+' '+tokens[2]
            outF.write(line)
        outF.close()
        return     
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 3:
        "Usage: Remap.py <required inputs>"
        print "  Required inputs:"
        print "  -i <infile> - input psfmap"
        print " -o <outfile> - output psfmap"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["infile=","outfile"])
    except getopt.GetoptError:
        "Usage: Remap.py <required inputs>"
        print "  Required inputs:"
        print "  -i <infile> - input psfmap"
        print " -o <outfile> - output psfmap"
        sys.exit(2)
    i_flag = 0
    o_flag = 0

    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            "Usage: Remap.py <required inputs>"
            print "  Required inputs:"
            print "  -i <infile> - input psfmap"
            print " -o <outfile> - output psfmap"
            sys.exit(2)
           
        elif opt in ("-i","--infile"):
            i_flag = 1
            infile = arg 
        elif opt in ("-o","--outfile"):
            o_flag = 1
            outfile = arg


    sumF = i_flag + o_flag  
    if sumF != 2:
        print "Usage: Remap.py <required inputs> \n"
        print "  Required inputs:"
        print "  -i <infile> - input psfmap"
        print " -o <outfile> - output psfmap"
        sys.exit(-2)
    remap = Remap(infile,outfile)
    remap.remap_psf_map()
    