#!/opt/local/bin/python
import numpy
import sys

def select(data,cols):
    return data[:,cols]

if (__name__=="__main__"):
    if (len(sys.argv) != 2):
        print "I need a file to elaborate"
        exit(1)

    filename = sys.argv[1]

    try:
        data = numpy.loadtxt(filename,dtype='str')
    except IOError:
        print "error opening the file"
        exit(1)

    try:
        input = raw_input("Enter columns to select: ")
        columns = map(int, input.split())
    except IOError:
        print "Non valid input"
        exit(1)

    try:
        ncol = len(columns)
    except:
        print "Non valid columns"
        exit(1)

    if ( max(columns) >= len(data[0,:])):
        print "column ",max(columns)," doesn't exist"
        print "max column ",len(data[0,:])
        exit(1)

    if (min(columns) < 0):
        print "column ",min(columns)," doesn't exist"
        exit()

    numpy.savetxt(filename+"_out",select(data,columns),fmt='%s',delimiter='\t')
