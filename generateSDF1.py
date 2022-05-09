# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:15:08 2021

@author: SmirnygaTotoshka
"""

#! /usr/bin/env python
# -*- coding: utf-8 -*-
#generate sdf file containing peptide with input sequence
#reference file must be special,ask it to author, they should be with script
import argparse
import os
import re
import sys
import copy
import pandas as pd

def isCorrectSequence(seq):
    alphabet = list("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.strip()
    for i in range(0,len(seq)):
        if seq[i] not in alphabet:
            return False
    return True

def parseReference(amino_acids, referencePath):
    atoms = []
    bonds = []
    counts = []
    #get all atoms and bonds of future peptide from reference file
    #reference files have to named after one-letter name amino acid and be "mol" format
    for acid in amino_acids:
        ref_name = acid + ".mol"
        ref = open(os.path.join(referencePath,ref_name),"r",encoding="utf-8").readlines()
        f = False
        f1 = False
        tmp_atoms = []
        tmp_bonds = []
        for line in ref:
            if "M  V30 BEGIN ATOM" in line.strip():
                f = True
                continue
            if "M  V30 END ATOM" in line.strip():
                f = False
                atoms.append(tmp_atoms)
            if "M  V30 BEGIN BOND" in line.strip():
                f1 = True
                continue
            if "M  V30 END BOND" in line.strip():
                f1 = False
                bonds.append(tmp_bonds)
            if ("M  V30 COUNTS" in line.strip()):
                counts.append(line.strip())
            if (f):
                tmp_atoms.append(line.strip())
            if (f1):
                tmp_bonds.append(line.strip())
    return (atoms,bonds,counts)

#change parameter in line
#returning edited line
#line - input line
#index - paramenter place
#cur_val - new value of parameter
#isCoord - true for ccordinate, write with format 4 digits after point
def setParam(line,index,cur_val,isCoord = False):
    characteristics = re.split("\s+",line)
    out="M  V30 "
    for i in range(2,len(characteristics)):
        if i == index:
            if isCoord:
                out += "{:.4f}".format(cur_val)+ " "
            else:
                out += str(cur_val)+ " "
        else:
            out += characteristics[i] + " "
    return out

#get value of parameter from line
#index must be > 2
def getParam(line,index):
    if index < 2:
        raise IndexError("parameters starts from position 2")
    characteristics = re.split("\s+",line)
    return characteristics[index]

#need delete OH from COOH for build protein skeleton
#+need shift number of atoms after OH. This OH always have number 5 in reference files
def shiftAtoms(atoms):
    shift_ref_atom = copy.deepcopy(atoms)
    for i in range(0,len(shift_ref_atom)):
        flag = False
        j = 0
        while j < len(shift_ref_atom[i]):
            if (int(getParam(shift_ref_atom[i][j],2)) == 5 and not flag and i != len(shift_ref_atom)-1):
                flag = True
                del shift_ref_atom[i][j]
                j -= 1
            elif flag:
                n = int(getParam(shift_ref_atom[i][j],2)) - 1
                shift_ref_atom[i][j] = setParam(shift_ref_atom[i][j],2,n)
            j += 1
    return shift_ref_atom

#According to atom shifting, we shift bonds and bind first C in ref molecule and alpha N in next ref
def shiftBonds(bonds,shift_atoms):
    shifted_bond = copy.deepcopy(bonds)
    shift = 1
    for i in range(0,len(shifted_bond)):
        j = 0
        while j < len(shifted_bond[i]):
            if (int(getParam(shifted_bond[i][j],2)) == 4):
                if (i != len(shifted_bond)-1):
                    shift += len(shift_atoms[i])
                    shifted_bond[i][j] = setParam(shifted_bond[i][j],5,shift)
                j+=1
                continue
            first = int(getParam(shifted_bond[i][j],4))
            second = int(getParam(shifted_bond[i][j],5))
            if(i != len(shifted_bond) - 1):
                if (first >=5):
                    first -= 1
                if (second >=5):
                    second -= 1
            shifted_bond[i][j] = setParam(shifted_bond[i][j],4,first)
            shifted_bond[i][j] = setParam(shifted_bond[i][j],5,second)
            j += 1
    return shifted_bond


#increase x coordinate for nice view and renumber atoms.
def niceView(shift_atoms):
    ready_atom = copy.deepcopy(shift_atoms)
    shift_num = len(shift_atoms[0])+1
    delta_x = 5
    for i in range(1,len(ready_atom)):
        j = 0
        while j < len(ready_atom[i]):
            ready_atom[i][j] = setParam(ready_atom[i][j],2,shift_num)
            cur_x = float(getParam(ready_atom[i][j],4))
            shift_x = cur_x + delta_x*i
            ready_atom[i][j] = setParam(ready_atom[i][j],4,shift_x,True)
            shift_num += 1
            j += 1
    return ready_atom

#renumber atoms in bonds
def renumberBonds(shift_bonds,ready_atoms):
    ready_bond = copy.deepcopy(shift_bonds)
    shift_num = len(ready_atoms[0])+1
    n = shift_num
    for i in range(1,len(ready_bond)):
        j = 0
        while j < len(ready_bond[i]):
            first = int(getParam(shift_bonds[i][j],4)) + shift_num-1
            second = int(getParam(shift_bonds[i][j],5)) + shift_num-1
            ready_bond[i][j] = setParam(ready_bond[i][j],4,first)
            if not (int(getParam(shift_bonds[i][j],2)) == 4 and i != len(shift_bonds)-1):
                ready_bond[i][j] = setParam(ready_bond[i][j],5,second)
            ready_bond[i][j] = setParam(ready_bond[i][j],2,n)
            n += 1
            j += 1
        shift_num += len(ready_atoms[i])
    return ready_bond

#generate_counts
def generateCounts(ready_atoms,ready_bonds):
    atoms = bonds = 0
    for i in range(0,len(ready_atoms)):
        atoms += len(ready_atoms[i])
    for i in range(0,len(ready_bonds)):
        bonds += len(ready_bonds[i])
    edit_count = "M  V30 COUNTS "
    edit_count += str(atoms) + " "+ str(bonds) +" 0 0 0\n"
    return edit_count

#write amino acid sequence as mol structure to sdf
def writeStructure(out,seq, ref):
    if isCorrectSequence(seq):
        ref_atoms, ref_bonds, ref_counts = parseReference(seq,ref)
        shift_atoms = shiftAtoms(ref_atoms)
        shift_bonds = shiftBonds(ref_bonds, shift_atoms)
        ready_atoms = niceView(shift_atoms)
        ready_bonds = renumberBonds(shift_bonds,ready_atoms)
        edit_count = generateCounts(ready_atoms,ready_bonds)

        out.write(seq + "\n")
        out.write("Generated by the script of A. Smirnov, Dept. Bioinformatics, Pirogov RNRMU.\n\n")
        out.write("  0  0  0     0  0            999 V3000\n")
        out.write("M  V30 BEGIN CTAB\n")
        out.write(edit_count)
        out.write("M  V30 BEGIN ATOM\n")
        for mol in ready_atoms:
            for atom in mol:
                out.write(atom.strip()+"\n")
        out.write("M  V30 END ATOM\n")
        out.write("M  V30 BEGIN BOND\n")
        for mol in ready_bonds:
            for bond in mol:
                out.write(bond.strip()+"\n")
        out.write("M  V30 END BOND\n")
        out.write("M  V30 END CTAB\n")
        out.write("M  END\n")
    else:
        print("Check your sequence " + seq + ". It must be in first line and contain one-letter names of amino acids.")
        sys.exit(1)

def writeRecord(out,tbl,i,seq_col,ref):
    sequence = tbl.loc[i,seq_col]
    writeStructure(out, sequence, ref)
    for c in tbl.columns:
        if c != seq_col:
            out.write(">  <" + c + ">\n")
            out.write(str(tbl.loc[i,c]) + "\n\n")
    out.write("$$$$\n")

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", help="directory containing reference mol files with amino acids")
    parser.add_argument("input", help="path to csv file with records")
    parser.add_argument("output", help="path to output SDF")
    parser.add_argument("sequence",help = "Column name with sequences. They will be transformed into mol format")
    args = parser.parse_args()

    referencePath = args.reference
    database = pd.read_csv(args.input,sep = ";",header = 0)
    sequence = str(args.sequence).strip()
    all = len(database.index)
    j = 0
    printProgressBar(j, all, prefix = 'Progress:', suffix = 'Complete', length = 50)
    with open(args.output, "w", encoding="utf-8") as out:
        printProgressBar(j, all, prefix = 'Progress:', suffix = 'Complete', length = 50)
        for i in database.index:
            writeRecord(out,database,i,sequence,referencePath)
            j = j + 1
            printProgressBar(j, all, prefix = 'Progress:', suffix = 'Complete', length = 50)
