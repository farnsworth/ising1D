#!/opt/local/bin/python

import numpy
import sys

if (len(sys.argv) != 2):
    print "I need a file to elaborate"
    exit()

filename = sys.argv[1]

try:
    data = numpy.loadtxt(filename,dtype='str')
except IOError:
    print "error opening the file"
    exit()

try:
    input = raw_input("Enter columns to select: ")
    columns = map(int, input.split())
except:
    print "Non valid input"
    exit()

try:
    ncol = len(columns)
except:
    print "Non valid columns"
    exit()

if ( max(columns) >= len(data[0,:])):
    print "column ",max(columns)," doesn't exist"
    print "max column ",len(data[0,:])
    exit()

if (min(columns) < 0):
    print "column ",min(columns),"doesn't exist"
    exit()

numpy.savetxt(filename+"_out",data[:,columns],fmt='%s',delimiter='\t')
