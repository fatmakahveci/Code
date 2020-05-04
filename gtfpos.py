#!/bin/python

def lineFormat(line):
    fields=line.replace(";","\t").replace(" ","\t").replace("\"","").replace(";","").strip().split('\t');
    return fields;

directory="/Users/fatma/Downloads/"
outputFile=open(directory+"exonsInAGene.txt","w")

file=open(directory+"f","r")

line=file.readline()
fields=line.split('\t')

prevLine=line
prevFields=fields

while len(line) != 0:
    start = fields[1]
    end = fields[2]
    if fields[3] == prevFields[3]:
        while prevFields[3] == fields[3]:
            if int(fields[1]) < int(prevFields[1]):
                start = fields[1]
            if int(fields[2]) > int(prevFields[2]):
                end = fields[2]
            prevFields=fields
            line=file.readline();
            if line != "":
                fields=line.split('\t')
            else:
                break;
    outputFile.write(prevFields[0] + '\t' + start + '\t' + end + '\t' + prevFields[3] + '\t' + prevFields[4])
    prevFields=fields
    prevLine=line
    line=file.readline()
    fields=line.split('\t')

outputFile.close()
