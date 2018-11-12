#!/usr/bin/env python
""" The script is used in grid production to modify previously created psfnap files
    to reflect new directory structure. The script should be run from the directory where psfmap files are.
    By N. Kuropatkin  08/30/2018
"""
import sys
import getopt
import glob


if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 3:
        print "Usage: MakePSFMaps.py <required inputs>"
        print "  Required inputs:"
        print " -c <medsconf> - like y3v02 "
        print " -b <basedir> - base directorey where psfs dir is"
        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:b:",["medsconf=","basedir"])
    except getopt.GetoptError:
        print "Usage: MakePSFMaps.py <required inputs>"
        print "  Required inputs:"
        print " -c <medsconf> - like y3v02 "
        print " -b <basedir> - base directorey where psfs dir is"
        sys.exit(2)
    c_flag = 0
    b_flag = 0

    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage: MakePSFMaps.py <required inputs>"
            print "  Required inputs:"
            print " -c <medsconf> - like y3v02 "
            print " -b <basedir> - base directorey where psfs dir is"
            sys.exit(2)
           
        elif opt in ("-c","--medsconf"):
            c_flag = 1
            medsconf = arg  
        elif opt in ("-b","--basedir"):
            b_flag = 1
            basedir = arg

    sumF = c_flag + b_flag
    if sumF != 2:
        print "Usage: MakePSFMaps.py <required inputs>"
        print "  Required inputs:"
        print " -c <medsconf> - like y3v02 "
        print " -b <basedir> - base directorey where psfs dir is"
        sys.exit(-2)
    " First make list of psfmap.dat files "
    listF = glob.glob("*psfmap*.dat")
    for fileN in listF:
        lines = []
        with open(fileN) as infile:
            for line in infile:
                tokens = line.split(' ')
                newline = tokens[0] +' '+ tokens[1]
                subline = tokens[2].split(medsconf)[1]
                subline = basedir+medsconf+subline
                lines.append(newline+' '+subline)                
        with open(fileN, 'w') as outfile:
            for line in lines:    
                outfile.write(line)
            
        