import sys
import numpy as np

# Contains functions to deal with pdb files. 
# PDB file format ATOM record
#   Columns     Data                            Justification    Data Type\n
#    1-4         “ATOM”                          character\n
#    7-11        Atom serial number              right            integer\n
#    13-16       Atom name                       left             character\n
#    17          Alternate location indicator                     character\n
#    18-20       Residue name                    right            character\n
#    22          Chain identifier                                 character\n
#    23-26       Residue sequence number         right            integer\n
#    27          Code for residue insertion                       character\n
#    31-38       X orthogonal Å coordinate       right            real (8.3)\n
#    39-46       Y orthogonal Å coordinate       right            real (8.3)\n
#    47-54       Z orthogonal Å coordinate       right            real (8.3)\n
#    55-60       Occupancy                       right            real (6.2)\n
#    61-66       Temperature factor              right            real (6.2)\n
#    73-76       Segment identifier              left             character\n
#    77-78       Element symbol                  right            character\n

## Function to parse line in pdb
#
# The parsing is based on the right alignment as detailed in the file descriptions. 
# The element might be empty depending on the file format.
#
# @return Returns the atom name, residue name, residue id, atom coordinates and element for each entry
def parse_line(line):
    atom_name = line[12:16].strip()
    res_name = line[17:20].strip()
    res_id = int(line[22:26])
    atom_xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    element = line[76:78].strip()
    return [atom_name, res_name, res_id, atom_xyz, element]

## Function to open pdb and parse it line by line
#
# Returns the atom data as dictionary,
# the number of atoms (natom) and a list of terminal atoms 
# Recognition of termini is based on TER statements in pdb file
#
# @param[in] inpfile - input file name
#
# @return Dictionary containg all pdb data, numebr of atoms and a list of terminal atoms
def get_pdb_data(inpfile):
    usemodels = False
    data = dict()
    termini = [1] #first atom is always a terminus
    natom = 0
    with open(inpfile,"r") as f:
        lines = f.readlines()
        for line in lines:
            if line[:5] == "MODEL":
                if usemodels:
                    return data,termini,natom
                else:
                    usemodels = True
            #parse atom information
            if line[:4] == "ATOM":
                natom += 1
                data[natom] = parse_line(line)
            #add entry for every terminal atom
            elif line[:3] == "TER":
                termini.append(natom)
            else:
                continue
    return data,termini,natom

## Obtain residue information
#
# Function to obtain the number of residues (resid, dual used as counter),
# the first and last atom of each residue as dictionary (res), a list of
# all residue names (resnames), and the coordinates.
#
# @param[in] natom - Number of atoms
# @param[in] data - pdb information from parsing
def get_residues(natom,data):
    res = {1: [1,0]}
    resid = 1
    for idx in range(natom):
        if resid < data[idx+1][2]:
            res[resid][1] = idx
            resid += 1
            res[resid] = [idx+1,0]
    res[resid][1] = natom
    resnames = list()
    for idx in range(resid):
        resnames.append(data[res[idx+1][0]][1])
    coords = list()
    for ridx in range(resid):
        for idx in range(res[ridx+1][0],res[ridx+1][1]+1):
            coords.append(data[idx][3])
    return resid, res, resnames, coords

def shift_resids(natom,data):
    resid = data[1][2]
    shift = resid - 1
    prev = 1
    for i in range(natom):
        curr = data[i+1][2] - shift
        if prev<curr:
            if prev+1==curr:
                prev = curr
            else:
                shift = shift + curr - prev - 1
                curr = data[i+1][2] - shift
                prev = curr
        data[i+1][2] = curr
    return data

## Function to get a list of elements from data dictionary
def get_elements(natom,data):
    elements = list()
    for idx in range(natom):
        elements.append(data[idx+1][4])
    return elements

## Function to get a list of atom names from data dictionary
def get_atomnames(natom,data):
    atomnames = list()
    for idx in range(natom):
        atomnames.append(data[idx+1][0])
    return atomnames

## Function to make sure we have correct termini
#
# The final entry should be the last atom
# In addition, we added the last atom of each molecule, but not the first
# Here, we fix this problem, and then transfer the list into a set to
# remove duplicates and then back into a sorted list to be useful
def fix_termini(natom,termini):
    if termini[-1] != natom:
        termini.append(natom)
    new_term = list()
    for t in termini:
        if t==1 or t==natom:
            continue
        else:
            new_term.append(t+1)
    return sorted(set(termini+new_term))

## Function to parse file and obtain all relevant data
def parse_pdb(inpfile):
    data, termini, natom = get_pdb_data(inpfile)
    data = shift_resids(natom,data)
    nres, res, resnames, coords= get_residues(natom,data)
    elements = get_elements(natom,data)
    atomnames = get_atomnames(natom,data)
    if elements[0] == "":
        for idx,name in enumerate(atomnames):
            elements[idx] = name[0]
    termini = fix_termini(natom, termini)
    return natom,nres,atomnames,elements,res,resnames,coords,termini


def find_peptide_bonds(nres,res,resnames,atomnames,termini):
    bonds = list()
#    print(termini)
    for i in range(nres-1):
        idx1 = i+1
        idx2 = i+2
        if (resnames[i]=="ACE") or (resnames[i+1]=="NME"):
            print("Skip caps, residue {} is {} and {} is {}".format(idx1,resnames[i],idx2,resnames[i+1]))
            continue
        if (res[idx1][1] in termini) and (res[idx2][0] in termini):
            print("Residues {}  and {} are both termini - no peptide bond".format(idx1,idx2))
            continue
        cidx = 0
        nidx = 0
        for j in range(res[idx1][0]-1,res[idx1][1]):
            if atomnames[j] == "C":
                cidx = j
                exit
        for j in range(res[idx2][0]-1,res[idx2][1]):
            if atomnames[j] == "N":
                nidx = j
                exit       
        if cidx!=0 and nidx!=0:
            bonds.append((cidx,nidx))
    return bonds

def atom_masses(elements):
    ## Defined atom masses
    atom_masses = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 16.000,"P": 30.974, "S": 32.06 }
    return [atom_masses[elem] for elem in elements]
    
def get_com(coords,mass):
    com = np.zeros(3)
    total_mass = 0.0
    for xyz,m in zip(coords,mass):
        com += m*np.array(xyz)
        total_mass += m
    com /= total_mass
    return com

def translate_to_origin(coords,mass):
    com = get_com(coords,mass)
    new_coords = list()
    for xyz in coords:
        new_coords.append(xyz - com)
    return new_coords

def calculate_principal_axes(coords,mass):
    inertia_tensor = np.zeros((3, 3))

    for r,m in zip(coords,mass):
        inertia_tensor[0, 0] += m * (r[1]**2 + r[2]**2)
        inertia_tensor[1, 1] += m * (r[0]**2 + r[2]**2)
        inertia_tensor[2, 2] += m * (r[0]**2 + r[1]**2)
        inertia_tensor[0, 1] -= m * r[0] * r[1]
        inertia_tensor[0, 2] -= m * r[0] * r[2]
        inertia_tensor[1, 2] -= m * r[1] * r[2]

     # Symmetrize the inertia tensor
    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    # Calculate eigenvectors (principal axes)
    eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)

    return eigenvectors  

def rotate_molecule(coords,masses):
    principal_axes = calculate_principal_axes(coords,masses)
    #Reorder them such that the z axis is the largest component
    principal_axes = principal_axes[:, [1, 2, 0]]
    new_coords = list()
    for r in coords:
        r = np.array(r)
        # Rotate the position using the principal axes
        new_r = np.dot(principal_axes.T, r)
        new_coords.append((new_r[0], new_r[1], new_r[2]))
    return new_coords

# the indexing is a bit over the place in this script - bonds are zero indexed
def get_bond_vecs(coords,bonds):
    vecs = list()
    for bond in bonds:
        at1 = bond[0]
        at2 = bond[1]
        vec = coords[at2,:] - coords[at1,:]
        vecs.append(vec)
    return vecs

def compute_angles(vec):
    z = np.array([0.0, 0.0, 1.0])
    magv = np.linalg.norm(vec)

    with np.errstate(divide="ignore", invalid="ignore"):
        cosp = np.dot(vec, z) / magv  

    cosp = np.clip(cosp, -1.0, 1.0)
    pitch = np.arccos(cosp)

    rot = np.arctan2(vec[1], vec[0])

    return pitch, rot

def get_suscept_tensor(p,r):
    sinp = np.sin(p)
    cosp = np.cos(p)
    sinr = np.sin(r)
    cosr = np.cos(r)
    sinp2 = sinp**2
    sinp3 = sinp**3
    cosp2 = cosp**2
    cosp3 = cosp**3
    sinr2 = sinr**2
    sinr3 = sinr**3
    cosr2 = cosr**2
    cosr3 = cosr**3

    X = np.zeros((3,3,3))

    X[0,0,0] = np.sum(-sinp3 * sinr3)
    X[0,0,1] = np.sum(cosr*sinp3*sinr2)
    X[0,0,2] = np.sum(cosp*sinp2*sinr2)
    X[0,1,0] = X[0,0,1]
    X[0,1,1] = np.sum(-cosr2*sinp3*sinr)
    X[0,1,2] = np.sum(-cosp*cosr*sinp2*sinr)
    X[0,2,0] = X[0,0,2]
    X[0,2,1] = X[0,1,2]
    X[0,2,2] = np.sum(-cosp2*sinp*sinr)

    X[1,0,0] = np.sum(cosr*sinp3*sinr2)
    X[1,0,1] = np.sum(-cosr2*sinp3*sinr)
    X[1,0,2] = X[0,1,2]
    X[1,1,0] = X[1,0,1]
    X[1,1,1] = np.sum(cosr3*sinp3)
    X[1,1,2] = np.sum(cosp*cosr2*sinp2)
    X[1,2,0] = X[0,1,2]
    X[1,2,1] = X[1,1,2]
    X[1,2,2] = np.sum(cosp2*cosr*sinp)

    X[2,0,0] = np.sum(cosp*sinp2*sinr2)
    X[2,0,1] = X[0,1,2]
    X[2,0,2] = np.sum(-cosp2*sinp*sinr)
    X[2,1,0] = X[0,1,2]
    X[2,1,1] = np.sum(cosp*cosr2*sinp2)
    X[2,1,2] = np.sum(cosp2*cosr*sinp)
    X[2,2,0] = X[2,0,2]
    X[2,2,1] = X[2,1,2]
    X[2,2,2] = np.sum(cosp3)

    X1 = X/len(p)
    X2 = X/X[2,0,0]

    return X1, X2

def get_data(coords,peptide_bonds):
    vecs = get_bond_vecs(coords,peptide_bonds)
    rot = list()
    pitch = list()
    for vec in vecs:
        p,r = compute_angles(vec)
        rot.append(r)
        pitch.append(p)
    rot = np.asarray(rot)
    pitch = np.asarray(pitch)
    meanp = np.mean(pitch)
    stdp = np.std(pitch)
    print("Mean pitch: ", meanp, " standard deviation: ", stdp)
    X1,X2 = get_suscept_tensor(pitch,rot)
    #print(X1)
    #print(X2)   
    rho1 = X2[2,2,2]
    rho2 = np.mean(np.cos(pitch)**3)/np.mean(np.cos(pitch)*np.sin(pitch)**2*np.sin(rot)**2)
    return rho1, rho2

if __name__ == "__main__":
    #parse data
    print("Parsing pdb file {}".format(sys.argv[1]))
    natom,nres,atomnames,elements,res,resnames,coords,termini = parse_pdb(sys.argv[1])
    peptide_bonds = find_peptide_bonds(nres,res,resnames,atomnames,termini)
    print("{} peptide bonds found".format(len(peptide_bonds)))
    #print(peptide_bonds)
    masses = atom_masses(elements)
    #translate to origin and then align to principal axes
    coords = translate_to_origin(coords,masses)
    coords = rotate_molecule(coords,masses)
    coords = np.array(coords)
    print("Max and min coordinates after alignment to z axis")
    print("X: ", np.max(coords[:,0]),np.min(coords[:,0]))
    print("Y: ", np.max(coords[:,1]),np.min(coords[:,1]))
    print("Z: ", np.max(coords[:,2]),np.min(coords[:,2]))
    rho1, rho2 = get_data(coords,peptide_bonds)
    cosp1 = np.sqrt(rho1/(rho1+2))
    cosp2 = np.sqrt(rho2/(rho2+2))
    print("Rho from X: {} | cosp from rho X: {} | Rho from approx.: {} | cosp from rho approx: {}".format(rho1, cosp1, rho2, cosp2))


    


