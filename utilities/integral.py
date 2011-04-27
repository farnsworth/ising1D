#!/opt/local/bin/python

import numpy
import sys

if (len(sys.argv) != 2):
    print "I need one file to elaborate"
    exit()

filename = sys.argv[1]

try:
    data = numpy.loadtxt(filename)
except IOError:
    print "error opening the file"
    exit()

try:
    column = int(raw_input("Enter column: "))
except:
    print "Non valid input"
    exit()

if (column >= len(data[0,:])):
    print "column ",column,"doesn't exist"
    exit()

integral = numpy.zeros(len(data[:,0]))

try:
    initval = float(raw_input("Enter initial time ["+`data[0,0]`+"]: "))
except:
    initval = data[0,0]

ifirst = -1

for i in range(1,len(data[:,column])):
    if (data[i,0]>=initval):
        integral[i] = integral[i-1] + data[i,column]*(data[i,0]-data[i-1,0])
        if (ifirst<0):
            ifirst = i-1

for i in range(ifirst,len(data[:,column])):
    integral[i] = integral[i]/(data[i,0]-data[ifirst,0])


final = numpy.column_stack((data,integral))
numpy.savetxt(filename+"_out",final,fmt='%.18e')
