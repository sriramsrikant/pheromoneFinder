##Code written by Sri (Sriram) Srikant
##Member of Murray and Gaudet Labs of Dept. of MCB
##Updated on 20200420

import sys
import os
from os import listdir
from os.path import isfile, join
import subprocess
import shlex

import traceback
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from collections import defaultdict
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import csv
import math
import operator
import itertools
from operator import attrgetter
from pprint import pprint

##Class definitions##
#This is to define the FASTAStructure
class FASTAStructure():
    def __init__(self, id, seq, start, stop, Cys):
        self.id = id
        self.seq = seq
        self.startPosition = start
        self.stopPosition = stop
        self.CysPosition = Cys
        
class FASTAStructureMinimal():
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

#This is to define a dictionary to contain the aa-distribution probabilities for scoring predicted pheromones
class ScoringSet(dict):
    def __missing__(self, key):
        return 0
 
#This is to define the Pheromone scoring structure
class pheromoneScoreClass():
    def __init__(self, id, seq, blastXevalue, blastNevalue, hydrophobicity, similarity, compositeScore):
        self.id = id
        self.seq = seq
        self.blastXevalue = blastXevalue
        self.blastNevalue = blastNevalue
        self.hydrophobicity = hydrophobicity
        self.similarityScore = similarity
        self.compositeScore = compositeScore

#This is to define the Pheromone candidate structure
class pheromoneCandidateClass():
    def __init__(self, id, seq, start, end, promoterSeq, flankSeq, promoterStart, promoterEnd):
        self.id = id
        self.seq = seq
        self.start = start
        self.end = end
        self.promoterSeq = promoterSeq
        self.flankSeq = flankSeq
        self.promoterStart = promoterStart
        self.promoterEnd = promoterEnd
 
##Function definitions##
def seqIdProcess(fasta_folder, dictionary_file):
    #This is to add the filename to the ids of chromosomes so that the species name (as contained in the file name) will be included in the chromosome identifiers.
    if not(os.path.isdir(fasta_folder)):
        print "Input provided is not a folder. Please use the folder containing all sequence files, genomes or proteins."
    else:
        seqFileList = [f for f in listdir(fasta_folder) if (f.endswith(".fasta") or f.endswith(".fas") or f.endswith(".fa") or f.endswith(".max.pep"))]
        writeCount = 1
        
        ##This is to create the folder to store the processed genome files.
        if not(os.path.isdir(str(fasta_folder[:-1])+"_Processed")):
            os.mkdir(str(fasta_folder[:-1])+"_Processed")
         
        ##This is to make the fungal genome to species dictionary.
        with open(dictionary_file, mode='r') as infile:
            reader = csv.reader(infile)
            fungalDictionary = {rows[0].split("\t")[0]: rows[0].split("\t")[2] for rows in reader}
        
        for g in seqFileList:
            records = open(os.path.join(fasta_folder, g), "rU")
            sequences = list(SeqIO.parse(records, "fasta"))
            records.close()
            
            #print(sequences[0])
            #print(sequences[0].description.replace(" ", "_||_"))
            #sys.exit()
            
            if (g.find(".f") >= 0): 
                speciesName = g.split(".f")[0]
            elif (g.find(".max.pep") >= 0):
                speciesName = g.split(".max.pep")[0]
                
            writeIndex = os.path.join(fasta_folder[:-1]+"_Processed/", g)
            outputHandle = open(writeIndex,"wb")
            if (g.find(".f") >= 0):
                for sequence_record in sequences:
                    outputHandle.write(">"+str(speciesName)+"|/|"+str(sequence_record.id)+"\n"+str(sequence_record.seq)+"\n")
            elif (g.find(".max.pep") >= 0):
                for sequence_record in sequences:
                    outputHandle.write(">"+str(speciesName)+"|/|"+str(sequence_record.description.replace(" ", "_||_"))+"\n"+str(sequence_record.seq)+"\n")
            outputHandle.close()
            writeCount += 1
          
        #I think I'm almost done!
        print "Done writing ", str(writeCount-1), " genome OR protein files in the "+str(fasta_folder[:-1])+"_Processed folder. Good luck!"

def genomeSplitter(genomeFile):
    #This is to split the genome into individual FASTA files in order to parallelize my pheromone finding script
    ##Using the biopython fasta parse we can read our fasta input
    records = open(genomeFile, "rU")
    genome = list(SeqIO.parse(records, "fasta"))
    records.close()

    genomeName = genomeFile.split('/')[-1]
    writeCount = 1
    for chromosome in genome:
        ##This is to open the file to read from within the loop.
        writeIndex = "Temp/split"+str(writeCount)+"_"+str(genomeName)
        outputHandle = open(writeIndex,"wb")
        outputHandle.write(">"+str(chromosome.id)+"\n"+str(chromosome.seq)+"\n")
        outputHandle.close()
        writeCount = writeCount + 1

    #I think I'm almost done!
    print "Done writing ", str(writeCount-1), " chromosome files in the Temp folder. Good luck!"

def scaffoldShortlist(scaffoldFolder, sizeCutoff):
    #This is to add the filename to the ids of chromosomes so that the species name (as contained in the file name) will be included in the chromosome identifiers.
    if not(os.path.isdir(scaffoldFolder) ):
        print "Input provided is not a folder. Please use the folder containing all genome files."
    else:
        scaffoldFileList = [f for f in listdir(scaffoldFolder) if (f.endswith(".fasta") or f.endswith(".fas") or f.endswith(".fa"))]
        writeCount = 1

        writeIndex = str("scaffoldList_Descriptions.txt")
        outputHandle = open(writeIndex,"wb")

        outputHandle2 = open("scaffoldList_large.txt", "wb") 
        
        for s in scaffoldFileList:
            records = open(scaffoldFolder+"/"+s, "rU")
            genome = list(SeqIO.parse(records, "fasta"))
            records.close()
            
            scaffoldName = s.split(".f")[0]
            scaffoldLength = 0
            for chromosome in genome:
                scaffoldLength = scaffoldLength + len(chromosome.seq)
            
            if scaffoldLength > int(sizeCutoff):
                outputHandle2.write(">"+str(s)+"\n")
                writeCount += 1

            outputHandle.write(">"+str(s)+"\t"+int(scaffoldLength))

        outputHandle.close()           
        outputHandle2.close()           
 
        #I think I'm almost done!
        print "Done selecting ", str(writeCount-1), " of ", str(len(scaffoldFileList)), " scaffold files in the "+str(genomeFolder)+" folder to identify scaffolds with ge ", str(sizeCutoff), ". Good luck!"

def scaffoldRunCheck(scaffoldFolder, scaffoldList):
    #This is to check the results of the job arrays of CAAX-Locator before the concatPutativePheromone can be run.
    if not(os.path.isdir(scaffoldFolder) ):
        print "Input provided is not a folder. Please use the folder containing all scaffold and CAAX-Locator output files."
    else:
        scaffoldCAAXs = [f[:f.rfind(".caax")] for f in listdir(scaffoldFolder) if (f.endswith(".caax"))]
        writeCount = 1

        scaffolds = [line.rstrip('\n') for line in open(str(scaffoldList))]

        if (len(scaffolds) < len(scaffoldCAAXs)):
            print "There are more CAAX-Locator results than scaffolds used. Something is horribly wrong. Printed excess files."
            writeIndex = str("excess_CAAXLocator-results.txt")
            outputHandle = open(writeIndex,"wb")

            exclude = set(scaffolds)
            listDiff = [i for i in scaffoldCAAXs if i not in exclude]
            for element in listDiff:
                outputHandle.write(str(element)+"\n")
            outputHandle.close()

        elif (len(scaffolds) == len(scaffoldCAAXs)):
            print "All CAAX-Locator results are accounted for. Lucky you!"

        elif (len(scaffolds) > len(scaffoldCAAXs)):
            print "There are some missing CAAX-Locator results. ", str(len(scaffolds)-len(scaffoldCAAXs)), " files are being written to the missedScaffoldList.txt to run as a batch. Good luck!"

            writeIndex = str("missedScaffoldList.txt")
            outputHandle = open(writeIndex,"wb")

            exclude = set(scaffoldCAAXs)
            listDiff = [i for i in scaffolds if i not in exclude]
            for element in listDiff:
                outputHandle.write(str(element)+"\n")
            outputHandle.close()

        #I think I'm almost done!
        print "Done inspecting ", str(scaffoldFolder), " for ", str(len(scaffolds)), " CAAX-Locator results from scaffold files. Good luck!"


def concatPutativePheromone(genomeFile,outputFile):
    #This is to concatenate all the putative pheromones from the genome after the parallelized chromosome search on multiple cores of the cluster. The outputFile has to be the same as what was used to write the splits in "Temp/"
    ##Using the biopython fasta parse we can read our fasta input
    records = open(genomeFile, "rU")
    genome = list(SeqIO.parse(records, "fasta"))
    records.close()

    chromosomeNumber = int(len(genome))
    #Define the hash to add the putative pheromone list to.
    putativePheromone = []

    for writeCount in xrange(chromosomeNumber):
        ##This is to open the file to write to with the loop.
        pheromoneList = []
        readIndex = "Temp/split"+str(writeCount+1)+"_"+str(outputFile)
        #print(readIndex)
        #exit()
        inputHandle = open(readIndex, "rU")
        pheromoneList = list(SeqIO.parse(readIndex, "fasta"))
        inputHandle.close()

        putativePheromone.extend(pheromoneList)

    ##This is to write the final output file with the complete set of the putative pheromones. 
    #outputFile = genomeFile[genomeFile[:rfind('.')]+"_caax.fasta"
    if not(os.path.isdir("CAAX-Locator")):
        os.mkdir("CAAX-Locator")
 
    outputHandle = open("CAAX-Locator/"+outputFile,"wb")
    for sequence in putativePheromone:
        outputHandle.write(">"+str(sequence.id)+"\n"+str(sequence.seq)+"\n")
    outputHandle.close()

    #I think I'm almost done!
    print "Done merging ", str(chromosomeNumber), " putative pheromone per chromosome output files in the Temp folder. This is a total of ", str(len(putativePheromone)), " putative pheromones. Good luck!"

def CAAXLocator(genomeFile,outputFile):
    #Create hash table to store identified CAAX-Stop motifs.
    putativePheromone=[]
    ##Create the core parameters of the system
    #These are parameters that determine the size of the excluded edges (in bp) of a chromosome in the search. These must be >= maxPheromoneSize (as defined below) in order to avoid definition errors
    scanStart = 300
    scanStop = 300
    #These are the limits on the size of the pheromone I set (in aa).
    minPheromoneSize = 20
    maxPheromoneSize = 100
    #These are the lists of amin-acids to lookup for the CAAX motif
    aliphaticAA1 = ["A","G","I","L","V","S","T","D"]
    aliphaticAA2 = ["A","G","I","L","V","M","S","T","N","Q","H"]
    farnesylXAA3 = ["A","G","I","L","V","M","C","S","N","Q"]

    ##Using the biopython fasta parse we can read our fasta input
    records = open(genomeFile, "rU")
    genome = list(SeqIO.parse(records, "fasta"))
    records.close()

    #I am going to scan through a genome chromosome-by-chromosome.
    for chromosome in genome:
        ##This is for the Watson strand
        chromosome.seq = chromosome.seq.upper()
        remainSeq = chromosome.seq[scanStart:-scanStop]
        currentPosition = scanStart
        stopPosition = 0
#        debugCount = 0
        print "The chromosome is:\t", chromosome.id, " of length ", len(chromosome.seq), "bp (usable length=", len(remainSeq), "bp)\nThere are ", remainSeq.count("TGA"), " STOP:TGA codons; ", remainSeq.count("TAG"), "STOP:TAG codons; ", remainSeq.count("TAA"), "STOP:TAA codons!"
        while ((remainSeq.count("TAG") + remainSeq.count("TGA") + remainSeq.count("TAA")) != 0 and currentPosition < (len(chromosome.seq)-scanStop)):
            #Now start by finding the first STOP codon including all TGA/TAG/TAA, in remainSeq; only if the corresponding STOP codon even remains, else this leads to a constant increment every loop.
            listStop = (remainSeq.find("TGA"), remainSeq.find("TAG"), remainSeq.find("TAA"))
            stopPosition = min([x for x in listStop if x !=-1])
#            print "The position of the STOP codon is ", currentPosition+stopPosition+1, ":\t", chromosome.seq[currentPosition+stopPosition-12:currentPosition+stopPosition+3], "\n"
            testRegion = chromosome.seq[currentPosition+stopPosition-3*maxPheromoneSize:currentPosition+stopPosition+3]
            testRegionTranslate = testRegion.translate()
#            print "The translated region of the ID'ed region is ", testRegionTranslate
#            exit()

            ##This is the part of the loop to start adding filters to ID the CAAX-STOP motif.
            if (testRegionTranslate[-minPheromoneSize:-1].count("*") == 0):
                if (testRegionTranslate[-5] == "C"):
                    if (testRegionTranslate[-4] in aliphaticAA1 and testRegionTranslate[-3] in aliphaticAA2):
                        if (testRegionTranslate[-2] in farnesylXAA3):
                            ##This is the next filter to locate the upstream START, making sure I see one before encountering another STOP. Also I need to keep counting STARTs until I hit a maxORF size.
                            if (testRegionTranslate[:-minPheromoneSize].count("M") != 0):
                                findingStartRegion = testRegionTranslate[:-minPheromoneSize]
                                countORF = 0
                                while (findingStartRegion.rfind("M") > findingStartRegion.rfind("*")):
                                    startPosition = findingStartRegion.rfind("M")
                                    #This is to make sure I count all possible ORFs per STOP that I find.
                                    countORF = countORF + 1
                                    putativePheromone.append(FASTAStructure(chromosome.id+"||[forward]["+str(countORF)+"]",testRegion[3*startPosition:],((currentPosition+stopPosition+1)+3*(-maxPheromoneSize+startPosition)),(currentPosition+stopPosition+1),((currentPosition+stopPosition+1)-12)))
                                    findingStartRegion = findingStartRegion[:startPosition]
                                    ##This is for STD-OUT comment at any time. 
                                    print "=-=-=\nPutative pheromones:\t#", len(putativePheromone), "\n The locus is at position [", putativePheromone[-1].startPosition, ",", putativePheromone[-1].stopPosition, "]bp on chromosome identified by ", putativePheromone[-1].id, "\nTranslated Sequence:\t", putativePheromone[-1].seq.translate(), "\n-=-=-"

            #Update variables for the next loop.
            remainSeq = chromosome.seq[currentPosition+stopPosition+3:-scanStop]
            currentPosition = currentPosition+stopPosition+3
#            print "Start at ", currentPosition+1, "with the sequence:\n", remainSeq, "\n"
#            debugCount += 1
#            if (debugCount == 3):
#                exit()

        ##This is for the Crick strand
        remainSeq = chromosome.seq.reverse_complement()[scanStart:-scanStop]
        currentPosition = scanStart
        stopPosition = 0
#        debugCount = 0
        print "The chromosome is:\t", chromosome.id, "||reverse of length ", len(chromosome.seq.reverse_complement()), "bp (usable length=", len(remainSeq), "bp)\nThere are ", remainSeq.count("TGA"), " STOP:TGA codons; ", remainSeq.count("TAG"), "STOP:TAG codons; ", remainSeq.count("TAA"), "STOP:TAA codons!"
        while ((remainSeq.count("TAG") + remainSeq.count("TGA") + remainSeq.count("TAA")) != 0 and currentPosition < (len(chromosome.seq.reverse_complement())-scanStop)):
            #Now start by finding the first STOP codon including all TGA/TAG/TAA, in remainSeq; only if the corresponding STOP codon even remains, else this leads to a constant increment every loop.
            listStop = (remainSeq.find("TGA"), remainSeq.find("TAG"), remainSeq.find("TAA"))
            stopPosition = min([x for x in listStop if x !=-1])
#            print "The position of the STOP codon is ", currentPosition+stopPosition+1, ":\t", chromosome.seq[currentPosition+stopPosition-12:currentPosition+stopPosition+3], "\n"
            testRegion = chromosome.seq.reverse_complement()[currentPosition+stopPosition-3*maxPheromoneSize:currentPosition+stopPosition+3]
            testRegionTranslate = testRegion.translate()
#            print "The translated region of the ID'ed region is ", testRegionTranslate
#            exit()

            ##This is the part of the loop to start adding filters to ID the CAAX-STOP motif.
            if (testRegionTranslate[-minPheromoneSize:-1].count("*") == 0):
                if (testRegionTranslate[-5] == "C"):
                    if (testRegionTranslate[-4] in aliphaticAA1 and testRegionTranslate[-3] in aliphaticAA2):
                        if (testRegionTranslate[-2] in farnesylXAA3):
                            ##This is the next filter to locate the upstream START, making sure I see one before encountering another STOP. Also I need to keep counting STARTs until I hit a maxORF size.
                            if (testRegionTranslate[:-minPheromoneSize].count("M") != 0):
                                findingStartRegion = testRegionTranslate[:-minPheromoneSize]
                                countORF = 0
                                while (findingStartRegion.rfind("M") > findingStartRegion.rfind("*")):
                                    startPosition = findingStartRegion.rfind("M")
                                    #This is to make sure I count all possible ORFs per STOP that I find.
                                    countORF = countORF + 1
                                    putativePheromone.append(FASTAStructure(chromosome.id+"||[reverse]["+str(countORF)+"]",testRegion[3*startPosition:],(len(chromosome.seq.reverse_complement())-((currentPosition+stopPosition+1)+3*(-maxPheromoneSize+startPosition))+1),(len(chromosome.seq.reverse_complement())-(currentPosition+stopPosition+1)+1),(len(chromosome.seq.reverse_complement())-((currentPosition+stopPosition+1)-12)+1))) 
                                    findingStartRegion = findingStartRegion[:startPosition]
                                    ##This is for STD-OUT comment at any time. 
                                    print "=-=-=\nPutative pheromones:\t#", len(putativePheromone), "\n The STOP is at position ", putativePheromone[-1].stopPosition, "bp on chromosome identified by ", putativePheromone[-1].id, "\nTranslated Sequence:\t", putativePheromone[-1].seq.translate(), "\n-=-=-"

            #Update variables for the next loop.
            remainSeq = chromosome.seq.reverse_complement()[currentPosition+stopPosition+3:-scanStop]
            currentPosition = currentPosition+stopPosition+3
#            print "Start at ", currentPosition+1, "with the sequence:\n", remainSeq, "\n"
#            debugCount += 1
#            if (debugCount == 3):
#                exit()

    ##Writing all putativePheromones to an output FASTA file as specified in the arguments of the function
    output_handle = open(outputFile, "w+")
    count=0
    for sequence in putativePheromone:
        output_handle.write(">"+sequence.id+"["+str(sequence.startPosition)+"-"+str(sequence.CysPosition)+"-"+str(sequence.stopPosition)+"][Translate:"+str(sequence.seq.translate())+"]\n"+str(sequence.seq)+"\n")
    output_handle.close()
    print "\nDone writing the putative Pheromone list with ", len(putativePheromone), "sequences. Good Luck!"

##Quick function to identify if there are any Asn/N residues in the working region of the candidates identified from CAAXLocator
def AsnLocator(inputFASTA):
     ##I am going to start by importing the FASTA list of all the putative pheromones from the output of CAAXLocator. 
    records = open(inputFASTA, "rU")
    putativePheromone = list(SeqIO.parse(records, "fasta"))
    records.close()

    refinedPheromone = []
#    maturePheromoneStartAA = ["Y","W","F","A","I","L","V"]
    maturePheromoneStartAA = ["Y","W","F","A","I","L","V","S","T","G","C","M"]
    asnPosition = 0
    for candidate in putativePheromone:
#        print candidate.seq.translate()
        if (len(candidate.seq.translate()) < 25):
            if (candidate.seq.translate()[:-10].count("N") != 0):
                asnPosition = candidate.seq.translate()[:-10].rfind("N")
                #if (candidate.seq.translate()[asnPosition+1] in maturePheromoneStartAA):
                refinedPheromone.append(candidate)
                asnPosition = 0
#                    print candidate.id, "\n", candidate.seq
        else:
            if (candidate.seq.translate()[-25:-10].count("N") != 0):
                asnPosition = candidate.seq.translate()[:-10].rfind("N")
                #if (candidate.seq.translate()[asnPosition+1] in maturePheromoneStartAA):
                refinedPheromone.append(candidate)
                asnPosition = 0
#                    print candidate.id, "\n", candidate.seq

    ##Writing all putativePheromones to an output FASTA file as specified in the arguments of the function
    if not(os.path.isdir("Asn-CAAX-Locator")):
        os.mkdir("Asn-CAAX-Locator")

    outputFile = str(inputFASTA).split("/")[-1]+".asn"
    output_handle = open("Asn-CAAX-Locator/"+outputFile,"wb")
    for sequence in refinedPheromone:
        output_handle.write(">"+sequence.id+"\n"+str(sequence.seq)+"\n")
    output_handle.close()
    print "\nDone writing the refined pheromone candidate list with", len(refinedPheromone), "(of the original", len(putativePheromone), ") sequences. Good Luck!"


#Function to BLAST against Saccharomyces cerevisiae genome in order to prioritize the results based on that. I am going to also think about how to combine the different prioritization approaches I have.
#Function to use the abundance of amino-acids against the other predicted mature pheromones that I have identified thus far from the literature. I am going to use these to score the 10amino-acids upstream of Cys for each of the initial hits.
def scorePheromone(inputFASTA,knownPheromone,speciesName):
    ##This is to record the hydrophobicity score using a dictionary/hash.
    hydrophobicityKyteDoolittle = {
        'A': 1.800,
        'R': -4.500,
        'N': -3.500,
        'D': -3.500,
        'C': 2.500,
        'Q': -3.500,
        'E': -3.500,
        'G': -0.400,
        'H': -3.200,
        'I': 4.500,
        'L': 3.800,
        'K': -3.900,
        'M': 1.900,
        'F': 2.900,
        'P': -1.600,
        'S': -0.800,
        'T': -0.700,
        'W': -0.900,
        'Y': -1.300,
        'V': 4.200
    }

    ##I am going to use the presented knownPheromone FASTA list as a metric to quantify the aa-distribution in known pheromones and construct a proportional scoring system. I rely on this since it seems quite clear that the approach of using the KyteDoolittle scale is not appropriate.
    readHandle = open(knownPheromone, "rU")
    knownPheromone = list(SeqIO.parse(readHandle, "fasta"))
    readHandle.close()

#    print knownPheromone[1].seq[-16:-4].upper()
#    exit()

    ##I am going to use the length of the Saccharomycetales pheromones to judge the aa-distribution in the knownPheromones. I am a little worried that I am going to bias my search but lets see how this goes. Importantly this script assumes I am definitely going to have this averagePheromoneLength upstream of all candidate pheromones.
    averagePheromoneLength = 12
    all_aas = defaultdict(int)
    aaDistributionKnownPheromone = ScoringSet()
    totalCount = len(knownPheromone)*averagePheromoneLength
    for record in knownPheromone:
        analysisStructure = ProteinAnalysis(str(record.seq[-(4+averagePheromoneLength):-4].upper()))
        for aa, count in analysisStructure.count_amino_acids().iteritems():
            all_aas[aa] += count
    
    for aa in all_aas:
        aaDistributionKnownPheromone[aa] = float(all_aas[aa])/totalCount

#    print all_aas, "\n", aaDistributionKnownPheromone, "\n", totalCount

    ##I am going to start by importing the FASTA list of all the putative pheromones from the output of CAAXLocator. 
    records = open(inputFASTA, "rU")
    putativePheromone = list(SeqIO.parse(records, "fasta"))
    records.close()

    ##This is the hash for the pheromone scores.
    pheromoneScores = []

    ##I am going to run a BLASTX against Scerevisiae locally in order to find the best-hit e-values and also use a hydrophobicity scoring scale for the 8aa upstream of the Cys. One option is to score the probabilities based on observed pheromones and the other is to use a standard hydrophobicity scoring scale.
    outputXML = str(inputFASTA)[:-6]+"-blastx.xml"
#    print outputXML
#    exit()
    blastxCline = NcbiblastxCommandline(query=str(inputFASTA), db="$HOME/Research/BLASTdb/Scer", outfmt=5, out=outputXML)
    stdOut, stdErr = blastxCline()

    ##I am going to parse the output XML from the BLASTX of the putative pheromones and using the the e-value of the top result.
#    print outputXML
    try:
        readHandle = open(outputXML)
        BLASTXResults = NCBIXML.parse(readHandle)
        BLASTXResultsList = list(BLASTXResults)
        readHandle.close()
    except:
        print "The output XML file is missing!"
        exit()

    ##I am going to run a BLASTN against each of the species' identified CDS  locally in order to find the best-hit e-values and try and eliminate any of the sequences that match identified CDS. I am following a similar scoring system compared to that I have already constructed for the BLASTX search.
    outputXML = str(inputFASTA)[:-6]+"-blastn.xml"
    blastDB = "$HOME/Research/BLASTdb/"+str(speciesName)
#    print outputXML
#    exit()
    blastnCline = NcbiblastnCommandline(query=str(inputFASTA), db=blastDB, outfmt=5, out=outputXML)
    stdOut, stdErr = blastnCline()

    ##I am going to parse the output XML from the BLASTN of the putative pheromones and using the the e-value of the top result.
    try:
        readHandle = open(outputXML)
        BLASTNResults = NCBIXML.parse(readHandle)
        BLASTNResultsList = list(BLASTNResults)
        readHandle.close()
    except:
        print "The output XML file is missing!"
        exit()

#    print len(BLASTXResultsList)
#    pprint (vars(BLASTXResultsList[0]))
    loopCount = 0
    for BLASTXresult in BLASTXResultsList:
#        print BLASTXresult.alignments[0].hsps[0].expect
        sortedBLASTXOutput = []
        sortedBLASTXOutput = sorted(BLASTXresult.descriptions, key=attrgetter('e'))
        try:
            eScore = sortedBLASTXOutput[0].e
        except:
            eScore = 1000000
#        queryID = str(BLASTXresult.query)
        queryID = str(putativePheromone[loopCount].id)
#        querySeq = str(BLASTXresult.alignments[0].hsps[0].query)
        querySeq = putativePheromone[loopCount].seq
#        print "The sequence id is:\t", str(queryID), "&\t", str(putativePheromone[loopCount].id), "\n", str(querySeq), "\n", str(putativePheromone[loopCount].seq)
#        exit()
        innerloopCounter = 0
        hydrophobicityScore = 0
        similarityScore = 0
#        print querySeq.translate()[-13:-5]
        while innerloopCounter <= averagePheromoneLength:
            hydrophobicityScore = hydrophobicityScore + hydrophobicityKyteDoolittle[querySeq.translate()[-(5+averagePheromoneLength)+int(innerloopCounter)]]
            similarityScore = similarityScore + aaDistributionKnownPheromone[querySeq.translate()[-(5+averagePheromoneLength)+int(innerloopCounter)]]
            innerloopCounter += 1
        pheromoneScores.append(pheromoneScoreClass(queryID,querySeq,eScore,0,hydrophobicityScore,similarityScore,0))
#        print "The scoring is as follows:\n", pheromoneScores[-1].id, "\ne-value:\t", pheromoneScores[-1].blastXevalue, "\tHydrophobicity Kyte-Doolittle score:\t", pheromoneScores[-1].hydrophobicity, "\tSimilarity to known Pheromones:\t", pheromoneScores[-1].similarityScore, "\n", pheromoneScores[-1].seq
        loopCount += 1
#    print len([record.hydrophobicity for record in pheromoneScores])

    ##This is to get the BLASTN results against the species CDS.
    idList = []
    idList = [record.id for record in pheromoneScores]
#    print len(idList)
    for BLASTNresult in BLASTNResultsList:
#        print BLASTNresult.alignments[0].hsps[0].expect
        sortedBLASTNOutput = []
        sortedBLASTNOutput = sorted(BLASTNresult.descriptions, key=attrgetter('e'))
        try:
            eScore = sortedBLASTNOutput[0].e
        except:
            eScore = 1000000
#        queryID = str(putativePheromone[loopCount].id)
        listIndex = idList.index(str(BLASTNresult.query))
        pheromoneScores[listIndex].blastNevalue = eScore

    ##I am going to calculate the compositeScore for the list of putativePheromones. I am normalizing the BLASTX e-value to be 1 for most dissimilar to Scer hits, and for the highest hydrophobicity score relative to hits from this genome. On [160414] I implemented the similarity scoring system based on the aaDistribution in known pheromones. I am going to normalize this by just dividing by the averagePheromoneLength since each aa-score is a probability. NOTE: Theoretically I should weight this by the aa-distribution of the proteome of each species, but I am not implementing that now. 
    for record in pheromoneScores:
        normalizedBLASTX = 1 - math.exp(-(record.blastXevalue))
        normalizedBLASTN = 1 - math.exp(-(record.blastNevalue))
        normalizedHydrophobicity = record.hydrophobicity/max([item.hydrophobicity for item in pheromoneScores])
        normalizedSimilarityScore = record.similarityScore/averagePheromoneLength
        record.compositeScore = (normalizedBLASTX+normalizedBLASTN+normalizedSimilarityScore)/3
#        print record.id, record.compositeScore

    ##This is to print a tab-delimited file of the sorted list of scoring.
    sortedPheromoneScores = []
    sortedPheromoneScores = sorted(pheromoneScores, key=attrgetter('compositeScore'), reverse=True)
    
    outputFile = str(inputFASTA)[:-6]+"_scores.dat"
    output_handle = open(outputFile, "w+")
    output_handle.write("Candidate ID"+"\t"+"Composite score"+"\t"+"BLASTX e-Value"+"\t"+"BLASTN e-Value"+"\t"+"Known pheromone similarity score"+"\t"+"Kyte-Doolittle hydrophobicity"+"\t"+"Protein Sequence"+"\t"+"DNA Sequence"+"\n") 
    for record in sortedPheromoneScores:
        output_handle.write(">"+record.id+"\t"+str(record.compositeScore)+"\t"+str(record.blastXevalue)+"\t"+str(record.blastNevalue)+"\t"+str(record.similarityScore)+"\t"+str(record.hydrophobicity)+"\t"+str(record.seq.translate())+"\t"+str(record.seq)+"\n")
    output_handle.close()

    ##This is to write the pheromoneScores to an XML output. Importantly I would like to sort these hits using an expected composite score.
    scoreTreeRoot = ET.Element("scoreTreeRoot")
    for record in sortedPheromoneScores:
        id = ET.SubElement(scoreTreeRoot, 'id')
        id.text = record.id
        seq = ET.SubElement(id, 'seq')
        seq.text = str(record.seq)
        blastXevalue = ET.SubElement(id, 'blastXevalue')
        blastXevalue.text = str(record.blastXevalue)
        blastNevalue = ET.SubElement(id, 'blastNevalue')
        blastNevalue.text = str(record.blastNevalue)
        hydrophobicity = ET.SubElement(id, 'hydrophobicity')
        hydrophobicity.text = str(record.hydrophobicity)
        similarityScore = ET.SubElement(id, 'similarityScore')
        similarityScore.text = str(record.similarityScore)
        compositeScore = ET.SubElement(id, 'compositeScore')
        compositeScore.text = str(record.compositeScore)

    scoreTree = ET.ElementTree(scoreTreeRoot)
    print len(scoreTreeRoot)
#    ET.dump(scoreTreeRoot)
    outputScoringXML = str(inputFASTA)[:-6]+"_scores.xml"
    scoreTree.write(outputScoringXML)
 
#Function to pull out the promoters of candidates
def promoterExtract(input_file,genome_file):
    records = []
    #This is where the new file is read.
    sequence = open(input_file, "rU")
    records = list(SeqIO.parse(sequence, "fasta"))
    sequence.close()

    genome = []
    #This is where the new file is read.
    sequence = open(genome_file, "rU")
    genome = list(SeqIO.parse(sequence, "fasta"))
    sequence.close()
    print "Working with ", len(records), " candidates to extract promoters!"
    pheromonePromoters = []
    for record in records:
        #print record.id
        ID = re.split('[|]{2}', record.id)
        #print(ID[-1][1:])
        IDstep2 = re.split('[\[\]\-]', ID[-1])
        #print ID[0], "\n", IDstep2[1], "\t", IDstep2[3], "\n", IDstep2[5], "\t", IDstep2[6], "\t", IDstep2[7]
        if record.id.rfind("|||") != -1:
            chrID = ID[0]+"|"
        else:
            chrID = ID[0]
        #print chrID
        strand = IDstep2[1]
        candidateNo = IDstep2[3]
        start = int(IDstep2[5])
        end = int(IDstep2[7])
#        sys.exit() 

        for chromosome in genome:
            #print chromosome.id[:-1]
            if (chrID==chromosome.id):
                if (strand=="forward"):
                    try:
                        promoterSeq=chromosome.seq[-1001+start:start-1].upper()
                        promoterStart = -1001+start
                    except:
                        promoterSeq=chromosome.seq[:start-1].upper()
                        promoterStart = 0
                    #print promoterSeq
                    if start < 1001 and end+1000 < len(chromosome.seq):
                        flankSeq=chromosome.seq[:end+1000].upper()
                    if start > 1001 and end+1000 > len(chromosome.seq):
                        flankSeq=chromosome.seq[-1001+start:].upper()
                    if start < 1001 and end+1000 > len(chromosome.seq):
                        flankSeq=chromosome.seq[:].upper()
                    else:
                        flankSeq=chromosome.seq[-1001+start:end+1000].upper()
                    
                    promoterEnd = start-1
                elif (strand=="reverse"):
                    try:
                        promoterSeq=chromosome.seq[start:start+1000].reverse_complement().upper()
                        promoterStart = start+1000
                    except:
                        promoterSeq=chromosome.seq[start:].reverse_complement().upper()
                        promoterEnd = len(chromosome.seq)-1
                    #print promoterSeq
                    if start+1000 > len(chromosome.seq) and end > 1001:
                        flankSeq=chromosome.seq[-1001+end:].reverse_complement().upper()
                    if start+1000 < len(chromosome.seq) and end < 1001:
                        flankSeq=chromosome.seq[:start+1000].reverse_complement().upper()
                    if start+1000 < len(chromosome.seq) and end < 1001:
                        flankSeq=chromosome.seq[:].reverse_complement().upper()
                    else:
                        flankSeq=chromosome.seq[-1001+end:start+1000].reverse_complement().upper()
                    
                    promoterEnd = start
                    
        #I am now going to finish the appending to the promoterClass
        #print promoterSeq
        #test = pheromoneCandidateClass(record.id,record.seq,str(start),str(end),promoterSeq,flankSeq,str(promoterStart),str(promoterEnd))
        pheromonePromoters.append(pheromoneCandidateClass(record.id,record.seq,str(start),str(end),promoterSeq,flankSeq,str(promoterStart),str(promoterEnd)))

##I am going to write out the promoters in a separate file.
    #I am going to make an Output folder to store the output files for promoters.
    if not(os.path.isdir("PromoterExtract")):
        os.mkdir("PromoterExtract")

    outputFile = "PromoterExtract/"+input_file.split("/")[-1]+".promoters"
    output_handle = open(outputFile, "w+")
    for record in pheromonePromoters:
        output_handle.write(">"+str(record.id)+"\n"+str(record.promoterSeq)+"\n")
    output_handle.close()

    outputFile = "PromoterExtract/"+input_file.split("/")[-1]+".flank"
    output_handle = open(outputFile, "w+")
    for record in pheromonePromoters:
        output_handle.write(">"+str(record.id)+"\n"+str(record.flankSeq)+"\n")
    output_handle.close()

    outputFile = "PromoterExtract/"+input_file.split("/")[-1]+".promoters.csv"
    output_handle = open(outputFile, "w+")
    output_handle.write("Candidate ID"+"\t"+"DNA Sequence"+"\t"+"Protein Sequence"+"\t"+"ORF Start"+"\t"+"ORF End"+"\t"+"Promoter Sequence"+"\t"+"Promoter Start"+"\t"+"Promoter End"+"\n") 
    for record in pheromonePromoters:
        output_handle.write(str(record.id)+"\t"+str(record.seq)+"\t"+str(record.seq.translate())+"\t"+str(record.start)+"\t"+str(record.end)+"\t"+str(record.promoterSeq)+"\t"+str(record.promoterStart)+"\t"+str(record.promoterEnd)+"\n")
    output_handle.close()

    ##This is to write the pheromoneScores to an XML output. Importantly I would like to sort these hits using an expected composite score.
    pheromonePromoterRoot = ET.Element("pheromonePromoterRoot")
    for record in pheromonePromoters:
        id = ET.SubElement(pheromonePromoterRoot, 'id')
        id.text = record.id
        seq = ET.SubElement(id, 'seq')
        seq.text = str(record.seq)
        proteinSeq = ET.SubElement(id, 'proteinSeq')
        proteinSeq.text = str(record.seq.translate())
        ORFStart = ET.SubElement(id, 'ORFStart')
        ORFStart.text = str(record.start)
        ORFEnd = ET.SubElement(id, 'ORFStart')
        ORFEnd.text = str(record.end)
        promoterSeq = ET.SubElement(id, 'promoterSeq')
        promoterSeq.text = str(record.promoterSeq)
        flankSeq = ET.SubElement(id, 'flankSeq')
        flankSeq.text = str(record.flankSeq)
        promoterStart = ET.SubElement(id, 'promoterStart')
        promoterStart.text = str(record.promoterStart)
        promoterEnd = ET.SubElement(id, 'promoterEnd')
        promoterEnd.text = str(record.promoterEnd)

    pheromonePromoterTree = ET.ElementTree(pheromonePromoterRoot)
    #print len(pheromonePromoterRoot)
#    ET.dump(pheromonePromoterRoot)
    outputPromoterXML = "PromoterExtract/"+input_file.split("/")[-1]+".promoters.xml"
    pheromonePromoterTree.write(outputPromoterXML)

    print "Finished printing all output!"

##This is a function to print out the statistics of all the genomes analyzed and the numbers of candidates identified after each panning round.
def candidateStats(genomeFolder, CAAXFolder, AsnCAAXFolder):
    #This is to check the statistics of the Asn...CAAX candidates from my algorithm.
    if (not(os.path.isdir(genomeFolder)) or not(os.path.isdir(CAAXFolder)) or not(os.path.isdir(AsnCAAXFolder))):
        print "One of the inputs provided is not a folder. Please use the folders containing the genomes (.fasta), CAAX candidates (.caax) and Asn-CAAX candidates (.asn) respectively."
    else:
        genomes = [f for f in listdir(genomeFolder) if (f.endswith(".fas") or f.endswith(".fasta") or f.endswith(".fa"))]
        #print(genomes)
        candidateCAAXList = [f[:f.rfind(".caax")] for f in listdir(CAAXFolder) if (f.endswith(".caax"))]
        #print(candidateCAAXList)
        candidateAsnCAAXList = [f[:f.rfind(".caax.asn")] for f in listdir(AsnCAAXFolder) if (f.endswith(".asn"))]
        #print(candidateAsnCAAXList)
        
        if ((set(genomes) != set(candidateCAAXList)) or (set(genomes) != set(candidateAsnCAAXList))):
            print "The files don't match up!\n"
            if (set(genomes) != set(candidateCAAXList)):
                print "There are extra CAAX-STOP candidate files:"
                extra = []
                extra = [f for f in candidateCAAXList not in genomes]
                print"\n".join(extra)
                print"\nThere are missing results for the genome files:"
                missing = []
                missing = [f for f in genomes not in candidateCAAXList]
                print "\n".join(extra)
            if (set(genomes) != set(candidateAsnCAAXList)):
                print "There are extra Asn...CAAX-STOP candidate files:"
                extra = []
                extra = [f for f in candidateAsnCAAXList not in genomes]
                print "\n".join(extra)
                print "\nThere are missing results for the genome files:"
                missing = []
                missing = [f for f in genomes not in candidateAsnCAAXList]
                print "\n".join(extra)
        else:
            genomeStats = []
            for g in genomes:
                fileHandle = open(genomeFolder+"/"+g, "r")
                scaffolds = list(SeqIO.parse(fileHandle, "fasta"))
                fileHandle.close()
                
                print "Genome file:\t\t\t\t\t", g
                print "\tScaffolds:\t\t\t\t", str(len(scaffolds))
                genomeSize = 0;
                for s in scaffolds:
                    genomeSize = genomeSize + len(s.seq);
                print "\tGenome Size:\t\t\t\t", str(genomeSize), "bp"
                
                fileHandle = open(CAAXFolder+"/"+g+".caax", "r")
                candidateCAAXs = list(SeqIO.parse(fileHandle, "fasta"))
                fileHandle.close()
                print "\tCAAX-Stop candidates:\t\t\t", str(len(candidateCAAXs))
                
                fileHandle = open(AsnCAAXFolder+"/"+g+".caax.asn", "r")
                candidateAsnCAAXs = list(SeqIO.parse(fileHandle, "fasta"))
                fileHandle.close()
                print "\tAsn...CAAX-Stop candidates:\t\t", str(len(candidateAsnCAAXs))
                
                uniqueCandidates = []
                skipCount = 0
                for candidate in reversed(candidateAsnCAAXs):
                    #print(candidate.id)
                    if (uniqueCandidates != []):
                        if (candidate.id.split("||")[0] == uniqueCandidates[-1].id.split("||")[0]) and (candidate.id.split("||")[1].split("[")[1] == uniqueCandidates[-1].id.split("||")[1].split("[")[1]):
                            if candidate.id.split("||")[1].split("[")[3].split("-")[2] == uniqueCandidates[-1].id.split("||")[1].split("[")[3].split("-")[2]:
                                skipCount += 1
                            else:
                                uniqueCandidates.append(candidate)
                        else:
                            uniqueCandidates.append(candidate)
                    else:
                        uniqueCandidates.append(candidate)
                        
                print "\tCandidates from overcounted loci:\t", skipCount
                uniqueCandidates.reverse()
                print "\tCandidates at unique loci:\t\t", str(len(uniqueCandidates)), "\n"
                
                genomeStats.append([g.split(".")[0], len(scaffolds), genomeSize, len(candidateCAAXs), len(candidateAsnCAAXs), len(uniqueCandidates)])
                
            genomeStatsFrame = pd.DataFrame.from_records(genomeStats, columns=["Species", "No. of scaffolds", "Genome Size (bp)", "No. of CAAX-Stop candidates", 
                                                                               "No. of Asn...CAAX-Stop candidates", "No. of unique locus Asn...CAAX-Stop candidates"])
            genomeStatsFrame.to_csv("./genomeStats.csv", index=False)
            #I think I'm almost done!
            print "Done inspecting ", str(genomeFolder), " for ", str(len(genomes)), " genomes and their Asn...CAAX-Locator results from other input folders. Good luck!"

##Writing the entire code as a function that I can copy straight into pheromoneFinder.
##This function is modified to account for multiple candidates from the same STOP where one of them (and not just the largest candidate) might represent a candidate that has multiple
## copies. I've set this up by estimating pairwise scores of all candidates and then setting the score to 0 when the pair considered are candidates from the same STOP.
##I am also going to ouput the histogram of pairwise distances after correcting for candidates from the same locus and also 
def candidateCopyNumber(AsnCAAXFile):
    ##Importing the sequences in the Asn...CAAX candidate output file.
    try:
        #print "Does the code start?"
        inputHandle = open(AsnCAAXFile, "U")
    except:
        print "The Input file, ", AsnCAAXFile, " does not seem to exist. Check the input file to make sure it is the output of the AsnLocator (*.caax.asn)."
    else:
        ##Importing 
        candidatesAsnCAAX = list(SeqIO.parse(inputHandle, "fasta"))
        inputHandle.close()
        
        if (AsnCAAXFile.rfind(".caax.asn") == -1 or len(candidatesAsnCAAX) == 0):
            print "The Input file does not seem to be the Asn...CAAX file (no .caax.asn extension) or there seems to be no candidates in the file."
        else:
            print "Analyzing ", len(candidatesAsnCAAX), " Asn...CAAX candidates from ", AsnCAAXFile.split("/")[-1].split(".")[0], " for candidates that are potential gene copies."
            
            pairScoresList = []
            reducePairScoresList = []
            alignScores = np.zeros((len(candidatesAsnCAAX), len(candidatesAsnCAAX)))
            for iterator in range(len(candidatesAsnCAAX)):
                for iterator2 in range(iterator):
                    ##This line uses the pairwise identity of candidates as measured with the region upstream of the found Asn/N to the end of the candidate, which is expected to be the mature pheromone + the CAAX motif.
                    alignScores[iterator, iterator2] = 100*pairwise2.align.globalxx(candidatesAsnCAAX[iterator].seq.translate()[max(candidatesAsnCAAX[iterator].seq.translate()[:-10].rfind("N")-2,0):], candidatesAsnCAAX[iterator2].seq.translate()[max(candidatesAsnCAAX[iterator2].seq.translate()[:-10].rfind("N")-2,0):], score_only=True)/max(len(candidatesAsnCAAX[iterator].seq.translate()[max(candidatesAsnCAAX[iterator].seq.translate()[:-10].rfind("N")-2,0):]), len(candidatesAsnCAAX[iterator2].seq.translate()[max(candidatesAsnCAAX[iterator2].seq.translate()[:-10].rfind("N")-2,0):]))
                    ##This line measures the pairwise identity of candidates across the entire candidate. This turns out is not great given I couldn't find S.cerevisiae pheromones.
                    #alignScores[iterator, iterator2] = 100*pairwise2.align.globalxx(candidatesAsnCAAX[iterator].seq.translate(), candidatesAsnCAAX[iterator2].seq.translate(), score_only=True)/max(len(candidatesAsnCAAX[iterator].seq.translate()), len(candidatesAsnCAAX[iterator2].seq.translate()))
                    #pairScoresList.append([candidatesAsnCAAX[iterator].id, candidatesAsnCAAX[iterator2].id, candidatesAsnCAAX[iterator, iterator2]])
                    
                    if ((candidatesAsnCAAX[iterator].id.split("||")[0] == candidatesAsnCAAX[iterator2].id.split("||")[0]) and (re.split(r'-', re.split(r'[\[\]]', candidatesAsnCAAX[iterator].id.split("||")[1])[5])[1] == re.split(r'-', re.split(r'[\[\]]', candidatesAsnCAAX[iterator2].id.split("||")[1])[5])[1])):
                        #print "Skipping ", candidatesAsnCAAX[iterator].id, " and ", candidatesAsnCAAX[iterator2].id, " they are from the same locus!"
                        reducePairScoresList.append([candidatesAsnCAAX[iterator].id, candidatesAsnCAAX[iterator2].id, 0])
                        alignScores[iterator, iterator2] = 0
                    else:
                        reducePairScoresList.append([candidatesAsnCAAX[iterator].id, candidatesAsnCAAX[iterator2].id, alignScores[iterator, iterator2]])
            #print(alignScores)
            
            pairScores = pd.DataFrame.from_records(reducePairScoresList, columns=["candidate1", "candidate2", "pairwiseScore"])
            #print(pairScores[pairScores['pairwiseScore'] > 85])
            
            #Checking for output directory, and making it if necessary.
            if not(os.path.isdir("./CopyNumber")):
                os.mkdir("./CopyNumber")
            
            ##Plotting the histogram of pairwise scores.
            listPairID = np.tril(alignScores).flatten()
            listPairID_new = listPairID[np.isfinite(listPairID)]
            binPlot = range(int(math.ceil(max(listPairID_new)))+1)
            fig = mpl.figure.Figure(figsize=(8,5))
            ax = sns.distplot(pairScores['pairwiseScore'], bins=binPlot, kde=False, hist_kws={"histtype": "step", "linewidth": 1, "alpha": 1, "color": "black"})
            ax.axvspan(85, 100, alpha=0.3, color='red')
            plt.ylabel("Counts")
            plt.yscale("log")
            plt.tight_layout()
            plt.savefig("./CopyNumber/"+AsnCAAXFile.split("/")[-1].split(".")[0]+"_pairID-hist.pdf", dpi=300, transparent=True)
            
            #Output the names of the candidates that have a copy among other candidates (as judged by >85% match in an alignment).
            copyNumList = list(set(list(itertools.chain(pairScores[pairScores['pairwiseScore'] > 85].candidate1.unique(), pairScores[pairScores['pairwiseScore'] > 85].candidate2.unique()))))
            outputFile = "./CopyNumber/"+AsnCAAXFile.split("/")[-1]+".multicopy"
            outputHandle = open(outputFile, "w")
            for element in copyNumList:
                for record in candidatesAsnCAAX:
                    if (record.id == element):
                        outputHandle.write(">"+str(element)+"\n"+str(record.seq)+"\n")
            outputHandle.close()
            
            print "Done processing ", len(candidatesAsnCAAX), " candidates from ", AsnCAAXFile.split("/")[-1].split(".")[0], ". Good luck!\n"
            #return pairScores


##Function to concatenate Asn...CAAX-Stop candidates that are judged to be within a phylogenetic group where pheromones are expected to be conserved.
def phyloGroupCandidates(listFile, AsnCAAXFolder):
    ##The following is to read the listFile for the species candidates that I need to pool as a Phylo Group.
    with open(listFile, "r") as inputHandle:
        speciesFileList = []
        #print inputHandle.readline()
        speciesFileList = inputHandle.read().split("\r\n")
        
    speciesFileList.remove("")

    print speciesFileList, "\n which contains ", len(speciesFileList), " species!\n"
    phyloGroupCandidates = []
    for speciesFile in speciesFileList:
        #This is to pool only the largest candidate per CAAX-Stop locus of Asn...CAAX-Stop candidates from speies within a PhyloGroup.
        inputHandle = open(AsnCAAXFolder+"/"+str(speciesFile+".unique"), "r")
        #This is when I pool all Asn...CAAX-Stop candidates from species within a PhyloGroup.
        #inputHandle = open(AsnCAAXFolder+"/"+str(speciesFile), "r")
        records = list(SeqIO.parse(inputHandle, "fasta"))
        inputHandle.close()
        
        ##Extending master list.
        phyloGroupCandidates.extend(records)
    
    print "Concatenated ", len(phyloGroupCandidates), " candidates in the ", listFile.split("/")[-1].split("."), " phylogenetic group.\n"
    
    ##Making Output directory.
    if not(os.path.isdir("./PhyloGroup_Asn-CAAX-Candidates")):
        os.mkdir("./PhyloGroup_Asn-CAAX-Candidates")
    
    outputFile = "./PhyloGroup_Asn-CAAX-Candidates/"+listFile.split("/")[-1].split(".")[0]+".fas.caax.asn"
    outputHandle = open(outputFile, "w")
    for candidate in phyloGroupCandidates:
        outputHandle.write(">"+str(candidate.id)+"\n"+str(candidate.seq)+"\n")
    outputHandle.close()

##Function to output the largest candidate at each locus from the identified Asn-CAAX candidates.
def uniqueCandidateLoci(AsnCAAXFile):
    inputHandle = open(AsnCAAXFile, "r")
    candidates = list(SeqIO.parse(inputHandle, "fasta"))
    inputHandle.close()
    
    uniqueCandidates = []
    skipCount = 0
    for candidate in reversed(candidates):
        #print(candidate.id)
        if (uniqueCandidates != []):
            if (candidate.id.split("||")[0] == uniqueCandidates[-1].id.split("||")[0]) and (candidate.id.split("||")[1].split("[")[1] == uniqueCandidates[-1].id.split("||")[1].split("[")[1]):
                if candidate.id.split("||")[1].split("[")[3].split("-")[2] == uniqueCandidates[-1].id.split("||")[1].split("[")[3].split("-")[2]:
                    skipCount += 1
                else:
                    uniqueCandidates.append(candidate)
            else:
                uniqueCandidates.append(candidate)
        else:
            uniqueCandidates.append(candidate)
    
    print("Skipped %i candidates" % skipCount)
    uniqueCandidates.reverse()
    print(len(uniqueCandidates))
    
    #Checking for output directory, and making it if necessary.
    if not(os.path.isdir("./UniqueLoci")):
        os.mkdir("./UniqueLoci")
    
    ##printing output.
    outputHandle = open("./UniqueLoci/"+AsnCAAXFile.split("/")[-1]+".unique", "w")
    for candidate in uniqueCandidates:
        outputHandle.write(">"+str(candidate.id)+"\n"+str(candidate.seq)+"\n")
    outputHandle.close()
    
##Function to run hmmbuild on the alignments of the mating genes that are needed to find homologs across the available proteomes using hmmsearch.This code-base was tested localling in my WSL Ubuntu 18.04. The code uses the subprocess module to call the command line binary of hmmer and outputs the results into the same "hmm_dir" that contains the alignments.
###The function will take as inputs:
    ###hmm_dir: Folder containing the HMM profiles of mating genes that are descibed in "hmm_dict_file". File extenstion of HMM profile files are ".hmm"
    ###hmm_dict_file: Dictionary file of all mating gene homologs with built HMMs. Full file path (not relational) path needs to be provided ideally.
def hmmbuildRun(hmm_dir, hmm_dict_file):
    ####Base variables defined for the entire code.
    data_dir = "./"
    
    ####Defining the dictionary file for mating genes.
    #hmm_dict_file = hmm_dir+"Scer-matingGenes_phmmer_200416.csv"
    hmm_dict = pd.read_csv(hmm_dict_file)
    
    ####Looping through the dictionary file to build all HMM profiles.
    for index, record in hmm_dict.iterrows():
        input_file = os.path.join(hmm_dir, str(record['FileName'])+".afa")
        #print(input_file)
        output_file = os.path.join(hmm_dir, str(record['GeneName'])+"_"+str(record['FileName'])+".hmm")
        log_file = os.path.join(hmm_dir, str(record['GeneName'])+"_"+str(record['FileName'])+".log")
        
        hmmbuild_cmd = subprocess.call(["hmmbuild", "--amino", "-n", str(record['GeneName']), "-o", log_file, output_file, input_file])
    
    print("Check "+hmm_dir+" for HMM profiles built by hmmbuild for the "+str(len(hmm_dict))+" gene alignments!\n")

##Function to run hmmsearch of the HMM-profiles saved in the "hmm" folder against the proteomes of the species that are stored in the "protein_dir". This code-base was tested localling in my WSL Ubuntu 18.04. The code uses the subprocess module to call the command line binary of hmmer and outputs the results into new folders "hmmsearch_raw" and "hmmsearch_best".
###The variables that need to be hard-coded are "protein_dir" and "genome_dir".
###The function will take as inputs:
    ###hmm_dir: Folder containing the HMM profiles of mating genes that are descibed in "hmm_dict_file". File extenstion of HMM profile files are ".hmm"
    ###hmm_dict_file: Dictionary file of all mating gene homologs with built HMMs. Full file path (not relational) path needs to be provided ideally.
    ###species_prot: File name for the processed proteins file in the "protein_dir" folder. Makes the function ideal for run in an SLURM array job on Harvard's Odyssey cluster.
    ###hmmer_protein_dir: Folder containing fasta (".max.pep") file of proteins for each species. Copy of the base dir is saved at hmmer_protein_dir (defnied below) in order to use esl-sfetch (with indexing of the protein fasta). protein_dir="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/332_genome_annotations/pep_Processed"
def hmmsearchRun(hmm_dir, hmm_dict_file, species_prot, protein_dir ="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/332_genome_annotations/pep_Processed"):
    ####Base variables defined for the entire code.
    data_dir = "./"
    hmmer_protein_dir = "/home/ssrikant/Research/20200415_pheromoneFinder/seqfiles/hmmer_proteins"
    
    ####Defining the dictionary file for mating genes.
    hmm_dict = pd.read_csv(hmm_dict_file)
    
    ####Running hmmsearch against the protein files of the various species files.
    for index, record in hmm_dict.iterrows():
        input_prot_database = os.path.join(protein_dir, species_prot)
        hmm_file = os.path.join(hmm_dir, str(record['GeneName'])+"_"+str(record['FileName'])+".hmm")
        
        if not(os.path.isdir(os.path.join(data_dir, "hmmsearch_raw"))):
            os.mkdir(os.path.join(data_dir, "hmmsearch_raw"))
        hmmsearch_output_file = os.path.join(data_dir, "hmmsearch_raw", str(record['GeneName'])+"_0_"+species_prot.split(".max.pep")[0]+"_hmmsearch.out")
        hmmsearch_output_tblout = os.path.join(data_dir, "hmmsearch_raw", str(record['GeneName'])+"_0_"+species_prot.split(".max.pep")[0]+"_hmmsearch.tbl")
        
        hmmsearch_cmd = subprocess.call(["hmmsearch", "-o", hmmsearch_output_file, "--tblout", hmmsearch_output_tblout, hmm_file, input_prot_database])
    
    ####Running esl-sfetch (indexing and search) on the protein fasta files in the "hmmer_protein_dir" which is a copy of the "protein_dir" folder.
    eslsfetch_index_cmd = subprocess.call(["esl-sfetch", "--index", os.path.join(hmmer_protein_dir, species_prot)])
    
    if not(os.path.isdir(os.path.join(data_dir, "hmmsearch_best"))):
        os.mkdir(os.path.join(data_dir, "hmmsearch_best"))
    output_fasta_file = os.path.join(data_dir, species_prot.split(".max.pep")[0]+"_hmmsearch.fasta")
    input_prot_database = os.path.join(hmmer_protein_dir, species_prot)
    print("Working on finding homologous mating genes of "+species_prot.split(".max.pep")[0]+".\n")
    
    homolog_protein_collection = []
    
    for index, record in hmm_dict.iterrows():
        hmmsearch_output_tblout = os.path.join(data_dir, "hmmsearch_raw", str(record['GeneName'])+"_0_"+species_prot.split(".max.pep")[0]+"_hmmsearch.tbl")
        temp_output_file = os.path.join(data_dir, "hmmsearch_best", "temp_0_"+species_prot.split(".max.pep")[0]+".fasta")
        #temp = "grep -v \"^#\" "+hmmsearch_output_tblout+" | awk '{print $1}' | esl-sfetch"+" -f "+input_prot_database+" - > "+temp_output_file
        #eslsfetch_cmd = subprocess.call([temp])
                
        command = '|'.join([
            "grep -v \"^#\" "+hmmsearch_output_tblout,
            "awk '{print $1}'",
            "esl-sfetch -f "+input_prot_database+" - > "+temp_output_file
        ])
        #print command
        ##This works exactly as needed but python strongly discourages using shell=True since it is not secure. I am going to go ahead with this temporarily.
        eslsfetch_cmd = subprocess.call([command], shell=True)
        #pipeline = Pipeline(command)
        #if not pipeline.run():
        #    print "ERROR: Pipeline failed"
        #else:
        #    print pipeline.output
        #sys.exit()
        
        ##The next step is to open the results file and append the top hit to a class for all mating homologous genes into a single file. 
        ##I can then delete the temp file of the output file.
        input_temp_seq_file = temp_output_file
        file_handle = open(input_temp_seq_file, "r")
        temp_candidates = list(SeqIO.parse(file_handle, "fasta"))
        file_handle.close()
        
        #print(temp_candidates[0].seq)
        #print(str(record['GeneName'])+"|/|"+temp_candidates[0].id+"\n"+temp_candidates[0].seq)
        if not temp_candidates:
            print("No homolog found for "+str(record['GeneName'])+" in "+species_prot.split(".max.pep")[0]+". Skipping!\n")
        else:
            homolog_protein_collection.append(FASTAStructureMinimal(str(record['GeneName'])+"|/|"+temp_candidates[0].id, temp_candidates[0].seq))
            print("Homolog found for "+str(record['GeneName'])+" in "+species_prot.split(".max.pep")[0]+". "+homolog_protein_collection[-1].id+" has been saved!")
        
        os.remove(temp_output_file)
    
    print("There are "+str(len(homolog_protein_collection))+" mating gene homologs identified in "+species_prot.split(".max.pep")[0]+". Good luck!")
    for homolog in homolog_protein_collection:
        print("\t"+homolog.id.split("|/|")[0])
        
    homolog_byspecies_output_file = os.path.join(data_dir, "hmmsearch_best", "homologs_0_"+species_prot.split(".max.pep")[0]+".fasta")
    output_file_handle = open(homolog_byspecies_output_file, "w")
    for homolog in homolog_protein_collection:
        output_file_handle.write(">"+str(homolog.id)+"\n"+str(homolog.seq)+"\n")
    output_file_handle.close()
    print("\nBest homologs saved to:"+ homolog_byspecies_output_file)

##Function to run makeblastdb on fungal genomes in order to use tBLASTN OR BLASTN with the python wrapper. This will be performed on a copy of the directory containing genome FASTA files since the original folder is often cycled through in my code. This code-base was tested localling in my WSL Ubuntu 18.04. The code uses the python wrapper of maskeblastdb from the packages provided by NCBI, and outputs the results into the same folder blastdb_genome_dir. This code is not meant to be run using an array since the processing and memory needed should not be high.
###The function will take as inputs:
    ###species_genome_list_file: A .txt file listing the contents of the blastdb_genome_dir for which BLAST databases are to constructed.
    ###blastdb_genome_dir: Folder containing the genome FASTA files that are listed in the .txt input file list. Output will be in the same folder. Ensure write-permissions are allowed. default="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/blastdb_332_genome_assemblies_Processed"
def makeblastdbRun(genome_file_list_file, blastdb_genome_dir="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/blastdb_332_genome_assemblies_Processed"):
    ####Defining base variables.
    data_dir = "./"
    with open(genome_file_list_file, 'r') as f:
        genome_file_list = f.read().splitlines()
    
    for species_genome in genome_file_list:
        blastdb_input_file = os.path.join(blastdb_genome_dir, species_genome)
        blastdb_input_file_type = "fasta"
        blastdb_dbtype = "nucl"
        blastdb_genome_title = species_genome.split(".f")[0]
        
        makeblastdb_cmd = subprocess.Popen(["makeblastdb", "-dbtype", blastdb_dbtype, "-input_type", blastdb_input_file_type, "-title", blastdb_genome_title, "-in", blastdb_input_file], 
                                      stdout=subprocess.PIPE)
        stdout = makeblastdb_cmd.communicate()
        print("Done making the blast database with the genome file of "+species_genome.split(".f")[0]+", the blastdb folder is in "+blastdb_genome_dir)
        
    print("\nFinished running the makeblastdb command on "+str(len(genome_file_list))+" fungal genomes. Good luck!\n")
    
##Function to run tBLASTN on the fungal genomes with the best homologs of mating genes found in each proteome by "hmmsearchRun" above. This function outputs the promoter (upstream 1000bp) region of the homologs of the mating genes. his code-base was tested localling in my WSL Ubuntu 18.04. This code uses the python wrapper of tBLASTN from the packages provided by NCBI, and makes a new output dir "tblastn" for the tabular output of tBLASTN and "promoters" for the FASTA file outputs of the promoters of the homologous mating genes. This function is written to 
###The function will take as inputs:
    ###species_genome: File name of the genome FASTA file that is going to be searched for the mating gene homologs.
    ###hmmsearch_dir: Folder containing the results of hmmsearch to find the best homologs of mating genes in fungal genomes by the hmmsearchRun function above.
    ###genome_dir: Folder containing the genome FASTA files that are listed in the .txt input file list. default="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/332_genome_assemblies_Processed"
    ###blastdb_genome_dir: Folder containing the genome FASTA files that are listed in the .txt input file list. Output will be in the same folder. Ensure write-permissions are allowed. default="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/blastdb_332_genome_assemblies_Processed"
def matingGeneHomologPromoters(species_genome, hmmsearch_dir, 
                               genome_dir = "/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/332_genome_assemblies_Processed",
                               blastdb_genome_dir="/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/blastdb_332_genome_assemblies_Processed"):
    ####Defining base variables.
    data_dir = "./"
    #blastdb_genome_dir = "/n/murraylab/Users/ssrikant/Cluster/fungalDB/DNAFungalDB/0_332yeast_genomes/blastdb_332_genome_assemblies_Processed"
    
    print("Finding the tBLASTN hits of hmmsearch hits of the homologous mating proteins against the corresponding genomes.\n")
    tblastn_query_file = os.path.join(hmmsearch_dir, "homologs_0_"+species_genome.split(".f")[0]+".fasta")    
    
    if not(os.path.isdir(os.path.join(data_dir, "tblastn"))):
        os.mkdir(os.path.join(data_dir, "tblastn"))
    tblastn_output_file = os.path.join(data_dir, "tblastn", species_genome.split(".")[0]+"_tblastn.tbl")
    
    tblastn_db = os.path.join(blastdb_genome_dir, species_genome)
    #print outputXML
    #exit()
    ##Output format 5 corresponds to an xml, while format 6 is tabular output (without comments)
    tblastn_cmd = NcbitblastnCommandline(query=tblastn_query_file, db=tblastn_db, outfmt=6, out=tblastn_output_file)
    stdOut, stdErr = tblastn_cmd()

    ##I am going to parse the output XML from the tBLASTN of the putative pheromones and using the the e-value of the top result.
    ###tblastn_results_raw = NCBIXML.parse(read_tblastn_handle)
    ###tblastn_results = list(tblastn_results_raw)
    ##The tabular output of tBLASTN has the results in the following columns [qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore]
    try:
        #read_tblastn_handle = open(tblastn_output_file)
        tblastn_column_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        tblastn_results = pd.read_csv(tblastn_output_file, sep="\t", names=tblastn_column_names)
        #read_tblastn_handle.close()
    except:
        print("The tabular output of tBLASTN is missing for "+species_genome.split(".f")[0]+"!")
        sys.exit()
        
    ####Importing the genome sequence and storing it to a dictionary.
    genome_dictionary = {}
    genome_fasta_handle = open(os.path.join(genome_dir, species_genome), "r")
    genome_fasta_read = list(SeqIO.parse(genome_fasta_handle, "fasta"))
    genome_fasta_handle.close()
    
    for chromosome in genome_fasta_read:
        genome_dictionary.update({chromosome.id: chromosome.seq})
    
    #print(tblastn_results['qseqid'].unique())
    homolog_promoters_fasta = []
    for query_seqid in tblastn_results['qseqid'].unique():
        temp_result_frame = pd.DataFrame()
        temp_result_frame = tblastn_results[tblastn_results['qseqid'] == query_seqid].sort_values(by=['evalue'], ascending=True).copy()
        #print(temp_result_frame)
        
        if (temp_result_frame['length'].iloc[0] == temp_result_frame['length'].max()):
            if (temp_result_frame['sstart'].iloc[0] < temp_result_frame['send'].iloc[0]):
                orf_start = temp_result_frame['sstart'].iloc[0] - 3*(temp_result_frame['qstart'].iloc[0] - 1)
                #print(orf_start)
                if (orf_start <= 1000):
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][0:orf_start-1]
                else:
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start-1001:orf_start-1]
            elif (temp_result_frame['sstart'].iloc[0] > temp_result_frame['send'].iloc[0]):
                orf_start = temp_result_frame['sstart'].iloc[0] + 3*(temp_result_frame['qstart'].iloc[0] - 1)
                if (orf_start >= (len(genome_dictionary[temp_result_frame['sseqid'].iloc[0]]) - 1000)):
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start:].reverse_complement()
                else:
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start:orf_start+1000].reverse_complement()
            elif (temp_result_frame['sstart'].iloc[0] > temp_result_frame['send'].iloc[0]):
                orf_start = temp_result_frame['sstart'].iloc[0] + 3*(temp_result_frame['qstart'].iloc[0] - 1)
                promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start:orf_start+1000].reverse_complement()
        else:
            print("The alignment with the lowest evalue is saved, but the longest match is:\n"
                  +str(temp_result_frame[temp_result_frame['length'] == temp_result_frame['length'].max()]['sseqid'].values)
                  +"\t"+str(temp_result_frame[temp_result_frame['length'] == temp_result_frame['length'].max()]['sstart'].values)
                  +"-"+str(temp_result_frame[temp_result_frame['length'] == temp_result_frame['length'].max()]['send'].values))
            if (temp_result_frame['sstart'].iloc[0] < temp_result_frame['send'].iloc[0]):
                orf_start = temp_result_frame['sstart'].iloc[0] - 3*(temp_result_frame['qstart'].iloc[0] - 1)
                if (orf_start <= 1000):
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][0:orf_start-1]
                else:
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start-1001:orf_start-1]
            elif (temp_result_frame['sstart'].iloc[0] > temp_result_frame['send'].iloc[0]):
                orf_start = temp_result_frame['sstart'].iloc[0] + 3*(temp_result_frame['qstart'].iloc[0] - 1)
                if (orf_start >= (len(genome_dictionary[temp_result_frame['sseqid'].iloc[0]]) - 1000)):
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start:].reverse_complement()
                else:
                    promoter_seq = genome_dictionary[temp_result_frame['sseqid'].iloc[0]][orf_start:orf_start+1000].reverse_complement()
        
        homolog_promoters_fasta.append(FASTAStructureMinimal(temp_result_frame['qseqid'].iloc[0], promoter_seq))
        print("Saving the 1000bp promoter region:\n"+homolog_promoters_fasta[-1].id+"\n"+homolog_promoters_fasta[-1].seq+"\n")
        
    if not(os.path.isdir(os.path.join(data_dir, "promoters"))):
        os.mkdir(os.path.join(data_dir, "promoters"))
    promoter_output_file = os.path.join(os.path.join(data_dir, "promoters", "promoters_0_"+species_genome.split(".f")[0]+".fasta"))
    output_file_handle = open(promoter_output_file, "w")
    for homolog in homolog_promoters_fasta:
        output_file_handle.write(">"+str(homolog.id)+"\n"+str(homolog.seq)+"\n")
    output_file_handle.close()
    
    print("Done printing promoters of "+str(len(homolog_promoters_fasta))+" mating homologous genes in "+species_genome.split(".f")[0]+".\n")
    
##Function to sort through exceptions
def formatExceptionInfo(maxTBlevel=5):
    cla, exc, trbk = sys.exc_info()
    excName = cla.__name__
    try:
        excArgs = exc.__dict__["args"]
    except KeyError:
        excArgs = "<no args>"
    excTb = traceback.format_tb(trbk, maxTBlevel)
    return (excName, excArgs, excTb)

##Main code starts here##
#The code has 1 functions:
#    'CAAXLocator': Identifies CAAX-Stop motif reference FASTA file.
#    'genomeIdProcess': Processes genomes by prefixing the genome file name to the chromosome IDs.
#    'genomeSplitter': Splits genome into individual FASTA files for the chromosomes.
#    'concatPutativePheromone': Concatenates results from the CAAX Locator per chromosome into a single FASTA file of putative pheromones.
#    'AsnLocator': This is to confirm the putative pheromones have an Asn in the region upstream of the CAAX box. 
#    'scorePheromone': This is to score the putative pheromone hits from CAAX Locator with BLASTX against Scerevisiae proteome & use a hydrophobicity scale for the 8 amino-acids upstream of Cys (mature region of pheromone).
#    'promoterExtract': This is to extract the pheromone sequence of the candidates after the CAAX-AsnLocator.
#    'candidateStats': This function analyzes the genomes, CAAX-Stop candidates and Asn...CAAX-Stop candidates to output a csv of stats that I can use for plotting.
#    'candidateCopyNumber': This function analyzes the Asn...CAAX-Stop candidates from an input file to output a histogram of pairwiseIDs of candidates (not from the same locus) and a csv list of the candidates with a copy in the genome.
#    'phyloGroupCandidates': This function takes as input lists of species that are expected to have conserved pheromone sequences and the folder with all the Asn...CAAX-Stop candidates. The output is a concatenated list of all candidates for the phylogenetic groups.
#    'uniqueCandidateLoci': This function takes as input the fasta file containing all the identified Asn...CAAX candidates from a genome and outputs the list of the largest candidates per CAAX locus. I need this function to test the pooled 'phyloGroupCandidates'. 
#    'hmmbuildRun': This function takes as input the directory containing the alignments of mating genes in Saccharomycetales[taxid:4892] (saved as .afa) with the details of mating genes stored in a dictionary file. Ouput HMM profiles are saved in the same directory.
#    'hmmsearchRun': This function takes as input the HMM profiles of mating genes and a dictionary (CSV file) describing the details of the genes. The output is a file with the protein sequence of the best hit homolog in each proteome of the species listed in prot_dir, saved to to the output directory "hmmserach_best" in the base directory.
#    'makeblastdbRun': This function takes as input the genome_file_list of the directory containing genome FASTA files. This directory is a copy of the original genome directory since I have other scripts that cycle through that directory, and this script will add files to the folder.
#    'matingGeneHomologPromoters': This function takes as input the FASTA files of the best homologs found by hmmsearch and uses tBLASTN to find these ORFs in the genome, and then extract the promoter (1000bp upstream) regions of these genes. This function is written to be run as a SLURM array job on Harvard Odyssey per genome. Output FASTA files containing promoter regions are output into a new folder 'promoters' in the base directory.

userParameters=sys.argv[1:]
print userParameters
try:
    if (userParameters[0]=="help"):
        print "The code has 12 function:"
        print "+\'caax\'\tIdentifies all the CAAX-Stop motifs in the FASTA file."
        print "+\'seqIdProcess\'\tProcesses genomes by prefixing the sequence file name, i.e. species name to the sequence record IDs."
        print "+\'genomeSplit\'\tSplits a genome into individual FASTA files for each chromosome, written into the Temp folder."
        print "+\'conCat\'\tConcatenates the individual FASTA files generated from the pheromone search run for  each chromosome/contig, written into the Temp folder. Be sure to preserve the same output name to prevent errors."
        print "+\'scaffoldShortlist\'\tShortlist scaffolds that are larger than the provided cutoff (provided as an integer). This helps prioritize scaffolds based on the size for limited computational resources."
        print "+\'scaffoldRunCheck\'\tChecks if the output files have been generated for all the scaffolds that are present in the input scaffold list. Generates a list of missing outputs that can be used as an input for the CAAXLocator SLURM script."
        print "+\'scorePheromone\'\tUses the results of the CAAXLocator to score the putative pheromones using BLASTXagainst the Scerevisiae proteome & also a hydrophobicity scale for the 8 amino-acids upstream of Cys."
        print "+\'promoterExtract\'\tUses the candidates from the CAAX-Asn-Locators to extract upstream 1000bp as a pheromone."  
        print "+\'candidateStats\'\tGenerates statistics of genomes and candidates of all genomes present in the input folder, given that output files .caax and .caax.asn are present in their corresponding folders." 
        print "+\'candidateCopyNumber\'\tIdentifies Asn...CAAX-Stop candidates (from input .caax.asn) file that have multiple copies in the genome as judged by pairwiseID of candidates (not counting candidates from the same CAAX-Stop)." 
        print "+\'phyloGroupCandidates\'\tIdentifies Asn...CAAX-Stop candidates from identified phylogenetic groups as identified by species within a divergence threshold from Shen 2018." 
        print "+\'uniqueCandidateLoci\'\tThis function takes as input the fasta file containing all the identified Asn...CAAX candidates from a genome and outputs the list of the largest candidates per CAAX locus. I need this function to test the pooled \'phyloGroupCandidates\'" 
        print "+\'hmmbuildRun\'\tThis function takes as input the directory containing the alignments of mating genes in Saccharomycetales[taxid:4892] (saved as .afa) with the details of mating genes stored in a dictionary file. Ouput HMM profiles are saved in the same directory."
        print "+\'hmmserachRun\'\tThis function takes as input the HMM profiles of mating genes and a dictionary (CSV file) describing the details of the genes. The output is a file with the protein sequence of the best hit homolog in each proteome of the species listed in prot_dir, saved to to the output directory \'hmmserach_best\' in the base directory."
        print "+\'maskblastdbRun\'\tThis function takes as input the genome_file_list of the directory containing genome FASTA files. This directory is a copy of the original genome directory since I have other scripts that cycle through that directory, and this script will add files to the folder."
        print "+\'matingGeneHomologPromoters\'\tThis function takes as input the FASTA files of the best homologs found by hmmsearch and uses tBLASTN to find these ORFs in the genome, and then extract the promoter (1000bp upstream) regions of these genes. This function is written to be run as a SLURM array job on Harvard Odyssey per genome. Output FASTA files containing promoter regions are output into a new folder \'promoters\' in the base directory."
 
    elif (userParameters[0]=="seqIdProcess"):
        if len(userParameters)==3:
            seqIdProcess(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for Sequence ID processing, for genomes or proteins!"
    elif (userParameters[0]=="genomeSplit"):
        if len(userParameters)==2:
            genomeSplitter(userParameters[1])
        else:
            print "There is a problem with the inputs for Genome Splitter!"
    elif (userParameters[0]=="conCat"):
        if len(userParameters)==3:
            concatPutativePheromone(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for result concatenator!" 
    elif (userParameters[0]=="scaffoldShortlist"):
        if len(userParameters)==3:
            scaffoldShortlist(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for shortlisting scaffolds!"  
    elif (userParameters[0]=="scaffoldRunCheck"):
        if len(userParameters)==3:
            scaffoldRunCheck(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for validating CAAX Locator results of scaffolds!" 
    elif (userParameters[0]=="caax"):
        if len(userParameters)==3:
            CAAXLocator(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for CAAX locator!" 
    elif (userParameters[0]=="Asn"):
        if len(userParameters)==2:
            AsnLocator(userParameters[1])
        else:
            print "There is a problem with the inputs for Asn locator!" 
    elif (userParameters[0]=="scorePheromone"):
        if len(userParameters)==4:
            scorePheromone(userParameters[1],userParameters[2],userParameters[3])
        else:
            print "There is a problem with the inputs for pheromone scoring!" 
    elif (userParameters[0]=="promoterExtract"):
        if len(userParameters)==3:
            promoterExtract(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for pheromone promoter extraction!"  
    elif (userParameters[0]=="candidateStats"):
        if len(userParameters)==4:
            candidateStats(userParameters[1],userParameters[2],userParameters[3])
        else:
            print "There is a problem with the inputs for evaluating candidate statistics!"  
    elif (userParameters[0]=="candidateCopyNumber"):
        if len(userParameters)==2:
            candidateCopyNumber(userParameters[1])
        else:
            print "There is a problem with the inputs to identify candidates with multiple copies!"  
    elif (userParameters[0]=="phyloGroupCandidates"):
        if len(userParameters)==3:
            phyloGroupCandidates(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs to pool Asn...CAAX-Stop candidates from species within the same phylogenetic group!"  
    elif (userParameters[0]=="uniqueCandidateLoci"):
        if len(userParameters)==2:
            uniqueCandidateLoci(userParameters[1])
        else:
            print "There is a problem with the inputs to identify the Asn...CAAX-Stop candidates with unique CAAX-Stop positions!"
    elif (userParameters[0]=="hmmbuildRun"):
        if len(userParameters)==3:
            hmmbuildRun(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs to run hmmbuild on the alignments of mating genes!"
    elif (userParameters[0]=="hmmsearchRun"):
        if len(userParameters)==5:
            hmmsearchRun(userParameters[1],userParameters[2],userParameters[3],userParameters[4])
        elif len(userParameters)==4:
            hmmsearchRun(userParameters[1],userParameters[2],userParameters[3])
        else:
            print "There is a problem with the inputs to run hmmsearch on fungal proteomes to find homologs of mating genes!"
    elif (userParameters[0]=="makeblastdbRun"):
        if len(userParameters)==3:
            makeblastdbRun(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs to run makeblastdb on the genome FASTA files of fungal genomes!"
    elif (userParameters[0]=="matingGeneHomologPromoters"):
        if len(userParameters)==5:
            matingGeneHomologPromoters(userParameters[1],userParameters[2],userParameters[3],userParameters[4])
        else:
            print "There is a problem with the inputs to run tBLASTN to find the best homologs of mating genes in fungal genomes and extract promoters and output into FASTA files per fungal genomes!"
 
    else:
        print "There is a problem with inputs!"
except:
    print "Aaah! There is a problem!", sys.exc_info()[0], "/n", formatExceptionInfo()

