#!/usr/bin/python

import numpy  as np
import argparse

parser = argparse.ArgumentParser(description='Compute the average of 2D Chebyshev Moments')
parser.add_argument('FileList', metavar='L', type=str, nargs='+',
                    help='The list of .chebmom2D files')
parser.add_argument('-o', '--output',  type=str ,
                    help='The name of the output file',default="out.chebmom2D")

filelist= parser.parse_args().FileList
outfile = parser.parse_args().output


i = 0;
data  = None;
header= None;
for inputfile in filelist:
    print("added file ",i)
    with open(inputfile) as f:
        current_header = f.readline()+f.readline()[:-1]#This remove the last newline char. 

        #Check if header had been read
        #and that its the same as the current one
        if header is not None:
            if header != current_header:
                print("The headers:" ,header, current_header, "do not match");
                print("Aborting Average");
                break;

        #If is the same or is had not been read
        #save it.
        header = current_header;

    #Creata a zero array for data if is not
    if data is None:
        data = 0*np.genfromtxt(inputfile, dtype=None, skip_header=2);

    #Then just summ it
    data+=np.genfromtxt(inputfile, dtype=None, skip_header=2); 
    i=i+1;
#Compute average
data/=len(filelist)    

#Save with header
np.savetxt(outfile, data, header=header,comments="")
