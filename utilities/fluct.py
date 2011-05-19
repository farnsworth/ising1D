#!/opt/local/bin/python

import numpy
import sys

def getFluct(datax,datay,delta):
    initialtime = datax[0] + delta
    av = 0.0
    av2 = 0.0
    elabdata=[]
    ndata = 0
    interval=[]

    for i in range(len(datay)):
        for val in interval:
            if (val[0]<datax[i]-delta):
                interval.remove(val)
                av2 -= val[1]*val[1]
                av -= val[1]
                ndata = ndata - 1
            else:
                break

        ndata = ndata + 1
        av2 += datay[i]*datay[i]
        av += datay[i]
        interval.append( [ datax[i],datay[i] ] )

        if ( datax[i]>=initialtime ):
            tmp = av/float(ndata)
            tmp2 = av2/float(ndata)
            elabdata.append( [datax[i]-delta/2,tmp,tmp2-tmp*tmp] )

    return numpy.array(elabdata)


if (__name__=='__main__'):
    if (len(sys.argv) != 2):
        print "I need a file to elaborate"
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

    try:
        delta = float(raw_input("Enter time interval: "))
    except:
        print "Non valid input"
        exit()

    numpy.savetxt(filename+"_out",getFluct(data[:,0],data[:,column],delta),fmt='%.18e')
