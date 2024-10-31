import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
import argparse
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
columns = int(round(xlength/0.256))
rows = int(round(ylength/0.2217))
if (columns >= 3) and (columns % 3 == 0):
    columns = columns
else:
    n = 0
    while (columns % 3 !=0):
        columns += 1
        n += 1
    
if (rows >= 4) and (rows % 2 == 0):
    rows = rows
else:
    rows += 1
 



positions = []
dist = 0
for i in range(rows):
    if i % 2 == 0:
        x = np.linspace(0, (0+(columns-1)*2.56), columns)
        for k in range(len(x)):
            positions.append([x[k], dist, 0])
        dist += 2.217

    else:
        x = np.linspace(-1.28, (-1.28+(columns-1)*2.56), columns)
        for k in range(len(x)):
            positions.append([x[k], dist, 0])
        dist += 2.217
w = mda.Universe.empty(n_atoms=len(positions), trajectory=True)
w.atoms.positions = positions
w.add_TopologyAttr('names')
w.atoms.names = [f'B{i}' for i in range(1, w.atoms.n_atoms + 1)]
w.add_TopologyAttr('resnames')
w.residues.resnames = ['GRA']
d = np.unique(w.atoms.positions[:,1])
x_atom1 = w.atoms[w.atoms.positions[:,1] == d[1]][0]
x_atom2 = w.atoms[w.atoms.positions[:,1] == d[1]][-1]
y_atom1 = w.atoms[w.atoms.positions[:,1] == d[0]][0]
y_atom2 = w.atoms[w.atoms.positions[:,1] == d[-2]][0]
dim = [calc_bonds(x_atom1.position, x_atom2.position) + 0.24595 * 10,  calc_bonds(y_atom1.position, y_atom2.position) + 0.2130*2*10, 100, 90, 90, 90]
w.trajectory.add_transformations(transformations.boxdimensions.set_dimensions(dim))
dim = w.trajectory.ts.triclinic_dimensions
box_center = np.sum(dim, axis = 0)/2
w.atoms.translate(box_center - w.atoms.select_atoms('resname GRA').center_of_geometry())
w.atoms.write(output+'.gro')


u = mda.Universe(output+'.gro')
c = np.unique(u.atoms.positions[:, 1])

u.atoms.masses = 36



# Setting mass of the virtual-site 0
for j in range(len(c)):
    if j == 0:
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        idx = group.atoms.indices
        gr = np.arange(idx[2], idx[-1]+1, 3)
        for k in gr:
            u.atoms[k].mass = 0
    elif (j != len(c) - 1) and (j % 2 != 0) and (j != 0):
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        gr = np.arange(1, len(group), 3)

        for k in gr:
            u.atoms[group[k].index].mass = 0
    elif (j != len(c) - 1) and (j % 2 == 0) and (j != 0):
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        gr = np.arange(2, len(group), 3)

        for k in gr:
            u.atoms[group[k].index].mass = 0

    elif j == len(c) - 1:
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        idx = group.atoms.indices
        gr = np.arange(idx[1], idx[-1]+1, 3)
        for k in gr:
            u.atoms[k].mass = 0


def regroup_lower(lst):
    final = []

    threes = [lst[i:i+3] for i in range(0, len(lst), 3)]
    for l, r in zip(threes, threes[1:]):
        final.append(l)
        final.append([l[-1], r[0]])
    if len(threes[-1]) > 1:
        final.append(threes[-1])

    data = []
    for i in final:
        if len(i) == 3:
            data.append(i)
        else:
            data.append(i)
    return data

def regroup_upper(lst):
    sublists = []
    data = []

    for idx in range(0, len(lst), 3):
        pair = lst[idx:idx+2]
        triplet = lst[idx+1:idx+4]

        if len(pair) == 2:
            sublists.append(pair)

        if len(triplet) == 3:
            sublists.append(triplet)
    for i in sublists:
        if len(i) == 3:
            data.append(i)
        else:
            data.append(i)
    return data

def regroup_side(lst):
    sub_lst = lst[1:]
    data = []
    idx = 0
    while idx < len(sub_lst):
        if len(sub_lst[idx: idx + 3]) == 3:
            data.append(sub_lst[idx: idx + 3])
        idx +=2
    return data


def hexagon(universe):
    idx_top = u.atoms[u.atoms.positions[:, 1] == c[0]]
    idx_bottom = u.atoms[u.atoms.positions[:, 1] == c[-1]]

    sel = u.atoms[(u.atoms.masses == 0)]
    sel1 = sel.atoms[sel.atoms.positions[:, 0]
                     == np.max(sel.atoms.positions[:, 0])]
    idx = u.atoms - idx_top - idx_bottom - sel1

    b = idx.atoms[idx.atoms.masses == 0].indices

    hexagon_indices = []
    for i in b:
        empty = []
        for j in u.atoms.indices:
            if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                empty.append(j)
        empty.append(i)
        hexagon_indices.append(empty)

    '''
    Hexagons along the horizontal, or the x-axis
    '''
    upper = u.atoms[u.atoms.positions[:, 1] == c[0]]
    lower = u.atoms[u.atoms.positions[:, 1] == c[-1]]
    upper_down = u.atoms[u.atoms.positions[:, 1] == c[1]]
    lower_up = u.atoms[u.atoms.positions[:, 1] == c[-2]]

    groups_upper = regroup_upper(upper.atoms.indices)
    groups_lower = regroup_lower(lower.atoms.indices)
    groups_upper_down = regroup_lower(upper_down.atoms.indices)
    groups_lower_up = regroup_upper(lower_up.atoms.indices)
    for i in range(len(groups_upper)):
        if i % 2 == 0:
            hexagon_indices.append([groups_lower_up[i][0], groups_lower_up[i][1], groups_lower[i]
                                   [0], groups_lower[i][2], groups_upper[i][0], groups_upper[i][1], groups_lower[i][1]])

        else:
            hexagon_indices.append([groups_lower[i][0], groups_lower[i][1], groups_upper[i]
                                   [0], groups_upper[i][2], groups_upper_down[i][0], groups_upper_down[i][1], groups_upper[i][1]])
            
    '''
    Hexagons along the y-axis, or the vertical
    '''
    right = [u.atoms[u.atoms.positions[:,1] == c[i]][-2].index if i % 2 == 0 else u.atoms[u.atoms.positions[:,1] == c[i]][-1].index for i in range(len(c))]
    left = [u.atoms[u.atoms.positions[:,1] == i][0].index for i in c]
    vss = [u.atoms[u.atoms.positions[:,1] == c[i]][-1].index for i in range(len(c)) if i % 2 == 0 ]
    vss_use = vss[1:]
    for i in range(len(regroup_side(right))):
        hexagon_indices.append([regroup_side(right)[i][0], regroup_side(left)[i][0], regroup_side(right)[i][1], regroup_side(left)[i][1], regroup_side(right)[i][2], regroup_side(left)[i][2], vss_use[i]])
        
    '''
    Last hexagon to be added for the remaining virtual site
    '''
    hexagon_indices.append([u.atoms[u.atoms.positions[:,1] == c[-1]][-1].index, u.atoms[u.atoms.positions[:,1] == c[-1]][0].index\
               ,u.atoms[u.atoms.positions[:,1] == c[0]][-2].index, u.atoms[u.atoms.positions[:,1] == c[0]][0].index, u.atoms[u.atoms.positions[:,1] == c[1]][-1].index, u.atoms[u.atoms.positions[:,1] == c[1]][0].index,\
               u.atoms[u.atoms.positions[:,1] == c[0]][-1].index])
            
    return hexagon_indices
    


vs = hexagon(u) # list of virtual sites, and the hexagon indices
hexagons = [hex[:-1] for hex in vs] # List of hexagon indices
def virtual_sites(vs):
    virtual_indices = []
    for hex in vs:
        virtual_indices.append([hex[-1], hex[0], hex[1], hex[4], hex[5]])
    return virtual_indices

v_site = np.array(virtual_sites(vs)) + 1 #numpy array for exclusions

"""
Uniformly distributing the masses
"""
number_of_real_particles = u.atoms[u.atoms.masses != 0].n_atoms
total_mass = u.atoms.n_atoms * 36
mass_of_each_real_particle = total_mass/number_of_real_particles

print(mass_of_each_real_particle)
real_particle_index = u.atoms[u.atoms.masses !=0].indices
for i in real_particle_index:
    u.atoms[i].mass = mass_of_each_real_particle 







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


def get_neighbors_angles(hexagons):
    
    def get_common_side(hexa1, hexa2):
        common = list(set(hexa1).intersection(set(hexa2)))
        if common:
            first_idx = hexa1.index(common[0])
            second_idx = hexa1.index(common[1])
            if first_idx > second_idx:
                common.reverse()
        return common
    
    angles = []
    for i in range(len(hexagons)):
        for j in range(i + 1, len(hexagons)):
            hexa1 = hexagons[i]
            hexa2 = hexagons[j]
            common_side = get_common_side(hexa1, hexa2)

            if(common_side): 
                hexa1_idx = [hexa1.index(atom) for atom in common_side]
                if hexa1_idx == [0,1]: #Top neighbor   
                    hexa2_idx = hexa2.index(common_side[0])
                    angle = [hexa2[hexa2_idx - 4], common_side[0], common_side[1], hexa1[hexa1_idx[1] + 4]]

                elif hexa1_idx == [4,5]: #Bottom neighbor
                    hexa2_idx = hexa2.index(common_side[1])
                    angle = [hexa1[hexa1_idx[0] - 4], common_side[0], common_side[1], hexa2[hexa2_idx + 4]]

                elif hexa1_idx == [0,2] or hexa1_idx == [1,3]: #Top right/left neighbor
                    hexa2_idx = hexa2.index(common_side[1])
                    angle = [hexa1[hexa1_idx[0] + 4], common_side[0], common_side[1], hexa2[hexa2_idx - 4]]

                elif hexa1_idx == [3,5] or hexa1_idx == [2,4]: #Bottom left/right neighbor
                    hexa2_idx = hexa2.index(common_side[0])
                    angle = [hexa1[hexa1_idx[1] - 4], common_side[1], common_side[0], hexa2[hexa2_idx + 4]]

                angles.append(angle)
                
    return angles

impropers = np.vstack(
    (get_angles(hexagons), get_neighbors_angles(hexagons)))
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
        f"  {i:<5}     TC5     0     GRA     B{i:<5}     {i:<5}     0     {u.atoms[i-1].mass}\n")



# Bonds

topology_file.write("\n[ bonds ]\n")
topology_file.write("; i         j        funct  length  kb\n")
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
        f"  {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     1    120     300\n")

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
    topology_file.write(f"{i[-1]+1:<3}     {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     {i[3]+1:<3}     {i[4]+1:<3}     {i[5]+1:<3}\n")

topology_file.close()



