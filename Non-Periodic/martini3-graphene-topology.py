import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
import numpy as np
import argparse
import math
from MDAnalysis import transformations

# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str, default='martini_graphene',
                    help='Name of the output, default = martini_graphene')
parser.add_argument("-x", "--xlength", type=float,
                    help='dimension along x in nm', default = 12.0)
                    
parser.add_argument("-y", "--ylength", type=float,
                    help='dimension along y in nm', default = 12.0)
parser.add_argument("-z", "--zlength", type=float,
                    help='dimension along z in nm', default = 12.0)

args = parser.parse_args()

output = args.output
xlength = args.xlength
ylength = args.ylength
zlength = args.zlength

rows = int(round(ylength/0.2217))
raw_columns = int(round(xlength/0.256))
if (rows >=3) and (rows % 2 != 0):
    rows = rows
else:
    rows = rows +1

if (raw_columns % 3 == 0):
    columns = raw_columns - 1
else:
    n = 0
    while (raw_columns % 3 != 0):
        raw_columns = raw_columns +1
        n = n +1
    columns = raw_columns - 1


#----------------#
# Structure File #
#----------------#


positions = []
dist = 0  # setting y = 0 in the beginning


for i in range(rows):
    if i % 2 == 0:  # if the row is indexed even number

        if (i == rows-1) or (i == 0):  # special treatment for the first and the last row
            # AP-> tn = a+(n-1)d
            x = np.linspace(0, (columns-1) * 2.56, num=columns)
            # array for deleting the elements from the first and the last row which are virtual sites
            d_array = np.array(np.arange(3, columns+1, 3)-1)
            # deleting those virtual sites from the first and the last row
            x = np.delete(x, d_array)
        else:  # if the row is even. but are not the first and the last one
            # AP-> tn = a+(n-1)d
            x = np.linspace(0, (columns-1) * 2.56, num=columns)

        for k in range(len(x)):  # Now looping every x-coordinate
            # Adding the x,y,z coordinate for the position
            positions.append([x[k], dist, 0])

        dist += 2.217  # Going 2.17 Angstroms down/up along y

    else:  # if the row is indexed odd number
        x = np.linspace(-1.28, -1.28+columns*2.56,
                        num=columns+1)  # AP-> tn = a+(n-1)d

        for k in range(len(x)):
            positions.append([x[k], dist, 0])

        dist += 2.217  # Going 2.13 Angstroms down/up along y


# Creating an empty Universe with atoms, and trajectory is True for writing positions
w = mda.Universe.empty(n_atoms=len(positions), trajectory=True)
w.atoms.positions = positions
w.add_TopologyAttr('names')
w.atoms.names = [f'B{i}' for i in range(1, w.atoms.n_atoms+1)]
w.add_TopologyAttr('resnames')
w.residues.resnames = ['GRA']
d = np.unique(w.atoms.positions[:,1])
x_atom1 = w.atoms[w.atoms.positions[:,1] == d[1]][0]
x_atom2 = w.atoms[w.atoms.positions[:,1] == d[1]][-1]
y_atom1 = w.atoms[w.atoms.positions[:,1] == d[0]][0]
y_atom2 = w.atoms[w.atoms.positions[:,1] == d[-2]][0]
dim = [calc_bonds(x_atom1.position, x_atom2.position) + 0.24595 * 10,  calc_bonds(y_atom1.position, y_atom2.position) + 0.2130*2*10, zlength*10, 90, 90, 90]
w.trajectory.add_transformations(transformations.boxdimensions.set_dimensions(dim))
dim = w.trajectory.ts.triclinic_dimensions
box_center = np.sum(dim, axis = 0)/2
w.atoms.translate(box_center - w.atoms.select_atoms('resname GRA').center_of_geometry())
w.atoms.write(output+'.gro')

w.atoms.write(output+'.gro')


u = mda.Universe(output+'.gro')
# finding the unique y-coordinate so that we can loop through each row of atoms
c = np.unique(u.atoms.positions[:, 1])
u.atoms.masses = 36    # Since it's a TC5 bead, we use a mass of 36 for each


'''
Identifying the indices of the virtual sites, and setting their mass to zero
'''


for j in range(len(c)):  # Looping through each row, defined by the y-coordinate
    if (j != 0) and (j != len(c)-1):
        if j % 2 != 0:
            group = u.atoms[u.atoms.positions[:, 1] == c[j]]
            gr = np.arange(1, len(group), 3)

            for k in gr:
                u.atoms[group[k].index].mass = 0

        else:
            group = u.atoms[u.atoms.positions[:, 1] == c[j]]
            gr = np.arange(2, len(group), 3)
            for k in gr:
                u.atoms[group[k].index].mass = 0


def hexagon(universe):
    ''' returns the list with index for the virtual site, the corresponding indices for hexagon '''

    b = u.atoms[u.atoms.masses == 0].indices
    hexagon_indices = []
    for i in b:
        empty = []
        for j in u.atoms.indices:
            if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                empty.append(j)
        empty.append(i)
        hexagon_indices.append(empty)
    return hexagon_indices

vs = hexagon(u) # list of virtual site, and the corresponding hexagon indices
hexagons = [hex[:-1] for hex in vs] # List of hexagon indices
def virtual_sites(vs):
    virtual_indices = []
    for hex in vs:
        virtual_indices.append([hex[-1], hex[0], hex[1], hex[4], hex[5]])
    return virtual_indices

v_site = np.array(virtual_sites(vs)) + 1 # numpy array for exclusions, and virtual sites

    




def bonds(hexagons):
    bonds = []
    for hex in hexagons:
        bond1 = sorted([hex[0], hex[1]])
        bond2 = sorted([hex[0], hex[2]])
        bond3 = sorted([hex[2], hex[4]])
        bond4 = sorted([hex[4], hex[5]])
        bond5 = sorted([hex[5], hex[3]])
        bond6 = sorted([hex[1], hex[3]])

        if bond1 not in bonds:
            bonds.append(bond1)
        if bond2 not in bonds:
            bonds.append(bond2)
        if bond3 not in bonds:
            bonds.append(bond3)
        if bond4 not in bonds:
            bonds.append(bond4)
        if bond5 not in bonds:
            bonds.append(bond5)
        if bond6 not in bonds:
            bonds.append(bond6)

    return bonds


"""
def add_bonds(hexagons):
    bonds = []
    for hex in hexagons:
        bond1 = sorted([hex[0], hex[5]])
        bond2 = sorted([hex[1], hex[4]])
        
        if bond1 not in bonds:
            bonds.append(bond1)
        if bond2 not in bonds:
            bonds.append(bond2)
        
    return bonds
    
"""


def add_bonds(hexagons):
    bonds = []
    for hex in hexagons:
        bond1 = sorted([hex[0], hex[4]])
        bond2 = sorted([hex[1], hex[5]])
        bond3 = sorted([hex[1], hex[2]])
        bond4 = sorted([hex[3], hex[4]])
        bond5 = sorted([hex[0], hex[3]])
        bond6 = sorted([hex[2], hex[5]])

        if bond1 not in bonds:
            bonds.append(bond1)
        if bond2 not in bonds:
            bonds.append(bond2)
        if bond3 not in bonds:
            bonds.append(bond3)
        if bond4 not in bonds:
            bonds.append(bond4)
        if bond5 not in bonds:
            bonds.append(bond5)
        if bond6 not in bonds:
            bonds.append(bond6)

    return bonds


def long_bond(hexagons):
    bonds = []
    for hex in hexagons:
        bond1 = sorted([hex[2], hex[3]])
        bond2 = sorted([hex[0], hex[5]])
        bond3 = sorted([hex[1], hex[4]])

        if bond1 not in bonds:
            bonds.append(bond1)
        if bond2 not in bonds:
            bonds.append(bond2)
        if bond3 not in bonds:
            bonds.append(bond3)
    return bonds
    




def find_angles(hexagons):
    angles = []
    for hex in hexagons:
        angle1 = [hex[0], hex[2], hex[4]]
        angle2 = [hex[2], hex[4], hex[5]]
        angle3 = [hex[4], hex[5], hex[3]]
        angle4 = [hex[5], hex[3], hex[1]]
        angle5 = [hex[3], hex[1], hex[0]]
        angle6 = [hex[1], hex[0], hex[2]]
        angles.append(angle1)
        angles.append(angle2)
        angles.append(angle3)
        angles.append(angle4)
        angles.append(angle5)
        angles.append(angle6)
    return angles




### Uniformly distributing the masses
number_of_real_particles = u.atoms[u.atoms.masses != 0].n_atoms
total_mass = u.atoms.n_atoms * 36
mass_of_each_real_particle = total_mass/number_of_real_particles


real_particle_index = u.atoms[u.atoms.masses != 0].indices
for i in real_particle_index:
    u.atoms[i].mass = mass_of_each_real_particle


# Below two functions are for the impropers

def get_angles(hexagons):
    angles = []
    for hex in hexagons:
        angle1 = [hex[2], hex[4], hex[0], hex[5]]
        angle2 = [hex[4], hex[0], hex[5], hex[1]]
        angle3 = [hex[0], hex[5], hex[1], hex[3]]
        angles.append(angle1)
        angles.append(angle2)
        angles.append(angle3)
    return angles


def get_between_hexagons_angles(hexagons, beads_per_row):
    def get_common_side(hexa1, hexa2):
        common = list(set(hexa1).intersection(set(hexa2)))
        common.sort()
        return common

    angles = []
    treated_side = []
    gap_for_after = int(beads_per_row / 3)
    gap_for_below = int(beads_per_row / 3) + int(beads_per_row/3 - 1)
    gap_for_before = int(beads_per_row / 3) - 1

    for i in range(len(hexagons)):
        hexa = hexagons[i]
        try:
            hexa_after = hexagons[i+gap_for_after]
        except IndexError:
            hexa_after = None
        try:
            hexa_below = hexagons[i+gap_for_below]
        except IndexError:
            hexa_below = None

        if i > 0:
            try:
                hexa_before = hexagons[i+gap_for_before]
            except IndexError:
                hexa_before = None
        else:
            hexa_before = None

        if hexa_below == hexa_before:
            hexa_before = None

        if hexa_after:
            side_after = get_common_side(hexa, hexa_after)
            if side_after:
                idx1 = hexa.index(side_after[0])
                idx2 = hexa_after.index(side_after[1])
                side = (side_after[0], side_after[1])
                if side not in treated_side:
                    treated_side.append(side)
                    angle = [hexa[idx1 - 2], side_after[0],
                             side_after[1], hexa_after[idx2 + 2]]
                    angles.append(angle)

        if hexa_below:
            side_below = get_common_side(hexa, hexa_below)
            if side_below:
                idx1 = hexa.index(side_below[0])
                idx2 = hexa_below.index(side_below[1])
                side = (side_below[0], side_below[1])
                if side not in treated_side:
                    treated_side.append(side)
                    angle = [hexa[idx1 - 4], side_below[0],
                             side_below[1], hexa_below[idx2 + 4]]
                    angles.append(angle)

        if hexa_before:
            side_before = get_common_side(hexa, hexa_before)
            if side_before:
                idx1 = hexa.index(side_before[0])
                idx2 = hexa_before.index(side_before[1])
                side = (side_before[0], side_before[1])
                if side not in treated_side:
                    treated_side.append(side)
                    angle = [hexa[idx1 - 2], side_before[0],
                             side_before[1], hexa_before[idx2 + 2]]
                    angles.append(angle)

    return angles


impropers = np.vstack(
    (get_angles(hexagons), get_between_hexagons_angles(hexagons, int(raw_columns))))
impropers = impropers + 1


#---------------#
# Topology File #
#---------------#


# Open the file for writing

topology_file = open(output+".itp", 'w')

# Variables


# Header

topology_file.write(
    "; \n;  Graphene topology\n; for the Martini3 force field\n;\n; created by martini3-graphene-topology.py\n;\n")
topology_file.write("; Roshan Shrestha\n; CNRS\n;\n\n")

topology_file.write("[ moleculetype ]\n")
topology_file.write("; molname	 nrexcl\n")
topology_file.write("  GRA           1")


# Atoms

topology_file.write("\n[ atoms ]\n")
topology_file.write("; nr	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n")
for i in range(1, u.atoms.n_atoms+1):
    topology_file.write(
        f"  {i:<5}     TC5     0     GRA     B{i:<5}     {i:<5}     0     {(u.atoms[i-1].mass):.2f}\n")
"""
# constraints
topology_file.write("\n[ constraints ]\n")
for i in constraints(hexagons):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.4919\n")
"""

# Bonds

topology_file.write("\n[ bonds ]\n")
topology_file.write("; i	 j	  funct	 length	 kb\n")
for i in bonds(hexagons):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.24595     30000\n")



topology_file.write("; short_bonds_across\n")
for i in add_bonds(hexagons):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.4260     30000\n")

topology_file.write("; long_bonds_across\n")
for i in long_bond(hexagons):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.4919     30000\n")


# Angles

topology_file.write("\n[ angles ]\n")
topology_file.write("; i	 j	 k	 funct	 angle	 force_k\n")
for i in find_angles(hexagons):
    topology_file.write(
        f"{i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     1    120     300\n")


# Improper Dihedrals

topology_file.write("\n[ dihedrals ]\n")
topology_file.write("; i	 j	 k	 l     funct	 ref.angle     force_k\n")
for i in impropers:
    topology_file.write(
        f"  {i[0]:<3}     {i[1]:<3}     {i[2]:<3}     {i[3]:<3}    2     180     200\n")


# Virtual sites

topology_file.write("\n[ virtual_sitesn ]\n")
topology_file.write("; site	 funct	 constructing atom indices\n")
for i in v_site:
    topology_file.write(
        f"  {i[0]:<3}     1     {i[1]:<3}     {i[2]:<3}     {i[3]:<3}    {i[4]:<3}\n")

# Exclusions
topology_file.write("\n[ exclusions ]\n")
for i in vs:
    topology_file.write(
        f"{i[-1]+1:<3}     {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}    {i[3]+1:<3}     {i[4]+1:<3}     {i[5]+1:<3}\n")

topology_file.close()
