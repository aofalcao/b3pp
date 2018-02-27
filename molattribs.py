#!/usr/bin/python
"""
Copyright (C) 2012-2017, Andre Falcao, University of Lisbon - LaSIGE
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
Please cite the authors in any work or product based on this material: 
Martins IF, Teixeira AL, Pinheiro L, Falcao AO. 2012. A Bayesian approach to in silico blood-brain barrier penetration modeling. J. Chem. Inf. Model., 2012, 52 (6), pp 1686 - 1697 DOI: 10.1021/ci300124c


"""
import re
import pybel
import openbabel
import sys

   
def matrix_fp(mol):
    params = ""
    atom_dict = mult_dict()
    mw = calc_mw(mol)
    params = params + str(mw) + ' '
    na = natomsH(mol)
    params = params + str(mw/na) + ' '
    nring = countring(mol)
    params = params + str(nring) + ' '
    nringb = countringb(mol)
    params = params + str(nringb) + ' '
    mol.OBMol.DeleteHydrogens()
    atom_multiplicity = atommultiplicity(mol)
    natoms = numatoms(mol)
    for a in range(0, natoms):
        atom = mol.OBMol.GetAtom(a+1)
        atype = atom_type(atom)
        atom_m = atom_multiplicity[a][0]
        akey = str(atype)+'-'+str(atom_m)
        i = 0
        for k, v in atom_dict:
            if k == akey:
                atom_dict[i][1] += 1
            i = i + 1
    for k, v in atom_dict:
        params = params + str(v)+ ' '

    fps = mol.calcfp(fptype='FP2')
    bits=fps.bits
    for b in range(1,1024):
        if b in bits: params+="1 "
        else: params+="0 "
    return params


    

def mult_dict():
    mult_dict = [["20-0",0], ["53-1",0], ["15-4",0], ["16-4",0], ["35-1",0], ["35-0",0], ["11-0",0], ["6-4",0], ["6-2",0], ["6-3",0], ["6-1",0], ["9-1",0], ["5-3",0], ["7-3",0], ["7-2",0], ["7-1",0], ["7-4",0], ["16-2",0], ["16-3",0], ["16-1",0], ["17-1",0], ["17-0",0], ["8-0",0], ["8-1",0], ["8-2",0]]
    return mult_dict

        
def atom_type(atom):
    atype = atom.GetAtomicNum() 
    return atype

def calc_mw(mol):
    molmw = mol.molwt   
    return molmw   

def natomsH(mol):
    mol.addh()
    numatoms = len(mol.atoms)
    return numatoms

def numatoms(mol):
    '''Count the number of atoms in a molecule; param: molecule; returns the number of atoms '''
    natoms = mol.OBMol.NumAtoms()
    return natoms

def numbonds(mol):
    '''Count the number of bonds in a molecule; param: molecule; returns the number of bonds '''
    nbonds = mol.OBMol.NumBonds()
    return nbonds

def countring(mol):
    cont_ring = len(mol.OBMol.GetSSSR())
    return cont_ring

def countringb(mol):
    '''Count the number of ring bonds in a molecule; param: molecule; returns the number of ring bonds '''
    cont_ringb = 0
    nbonds = numbonds(mol)
    for i in range(0, nbonds):
        bond = mol.OBMol.GetBond(i)
        if bond.IsInRing():
            cont_ringb = cont_ringb + 1

    return cont_ringb


def atommultiplicity(mol):

    natoms = numatoms(mol)
    nbonds = numbonds(mol)

    if nbonds > 0:

    #list of atoms multiplicity (position 0 - multiplicity of the atom od the bond, position
    #1 - max bond order and position 2 - min bond order)

        atom_multiplicity = [[0]*3 for i in range(natoms)]

        for i in range(0, nbonds):
                bond = mol.OBMol.GetBond(i)
                atom1 = bond.GetBeginAtomIdx()
                atom2 = bond.GetEndAtomIdx()
                atom_multiplicity [atom1 - 1][0] += 1
                atom_multiplicity [atom2 - 1][0] += 1


        for i in range(0, natoms):
                atom = mol.OBMol.GetAtom(i+1)
                max = 0
                if atom.HasBondOfOrder(3):
                        max = 3
                elif atom.HasBondOfOrder(2):
                        max = 2
                elif atom.HasBondOfOrder(1):
                        max = 1
                atom_multiplicity [i][1] = max
                min = 0
                if atom.HasBondOfOrder(1):
                        min = 1
                elif atom.HasBondOfOrder(2):
                        min = 2
                elif atom.HasBondOfOrder(3):
                        min = 3
                atom_multiplicity [i][2] = min


        return atom_multiplicity


def get_attribs(mol, typ):
	try:
		mymol=pybel.readstring(typ, mol)
		s_descs=matrix_fp(mymol)
	
		descs=mymol.calcdesc()
		s_descs+="%s %s %s %s " % (descs['sbonds'], descs['dbonds'], descs['tbonds'], descs['abonds'])
		s_descs+="%s %s %s\n" % (descs['logP'], descs['MR'], descs['TPSA'])
		return s_descs
	except:
		return "ERRO"
    

#here's where the mol string enters along with the type!
mol = sys.argv[1]
if mol[0:6]=="InChI=":
	print get_attribs(mol, "inchi")
else:
	print get_attribs(mol, "smi")
