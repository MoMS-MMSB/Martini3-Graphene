import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds
import argparse
import os
import sys
from pathlib import Path
from MDAnalysis import transformations

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from graphene_shared import (
    Profiler,
    build_row_groups,
    collect_hexagons_from_virtual_sites,
    make_temp_path,
    resolve_output_paths,
    unique_bonds_from_template,
    warn_gro_atom_index_wrap,
)


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
parser.add_argument(
    "--force",
    action="store_true",
    help="Overwrite existing output files if they already exist.",
)
parser.add_argument(
    "--output-dir",
    type=str,
    default=None,
    help="Directory where output files are written. Defaults to the path in --output.",
)
parser.add_argument(
    "--quiet",
    action="store_true",
    help="Suppress non-essential stderr output.",
)
parser.add_argument(
    "--profile",
    action="store_true",
    help="Print stage timing information to stderr.",
)

args = parser.parse_args()

output = args.output
xlength = args.xlength
ylength = args.ylength
zlength = args.zlength

if xlength <= 0 or ylength <= 0 or zlength <= 0:
    parser.error("xlength, ylength, and zlength must be positive numbers.")

profiler = Profiler(args.profile)

gro_path, itp_path = resolve_output_paths(output=output, output_dir=args.output_dir)
tmp_gro_path = make_temp_path(gro_path)
tmp_itp_path = make_temp_path(itp_path)

if not args.force:
    existing_outputs = [str(p) for p in (gro_path, itp_path) if p.exists()]
    if existing_outputs:
        parser.error(
            "Refusing to overwrite existing output file(s): "
            + ", ".join(existing_outputs)
            + ". Use --force to overwrite."
        )

columns = max(3, int(round(xlength / 0.256)))
while columns % 3 != 0:
    columns += 1

rows = max(4, int(round(ylength / 0.2217)))
if rows % 2 != 0:
    rows += 1
 
BOND_TEMPLATE = ((0, 1), (0, 2), (2, 4), (4, 5), (5, 3), (1, 3))
SHORT_BOND_TEMPLATE = ((0, 4), (1, 5), (1, 2), (3, 4), (0, 3), (2, 5))
LONG_BOND_TEMPLATE = ((2, 3), (0, 5), (1, 4))
EDGE_TEMPLATE = ((0, 1), (0, 2), (2, 4), (4, 5), (5, 3), (1, 3))




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
profiler.mark("built positions")
w = mda.Universe.empty(n_atoms=len(positions), trajectory=True)
w.atoms.positions = positions
w.add_TopologyAttr('names')
w.atoms.names = ['B'] * w.atoms.n_atoms
w.add_TopologyAttr('resnames')
w.residues.resnames = ['GRA']
w.add_TopologyAttr('resids')
w.residues.resids = [1]
warn_gro_atom_index_wrap(w.atoms.n_atoms, args.quiet)
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
w.atoms.write(str(tmp_gro_path))
profiler.mark("wrote temporary .gro")


u = mda.Universe(str(tmp_gro_path))
c = np.unique(u.atoms.positions[:, 1])
row_groups = build_row_groups(u, c)

u.atoms.masses = 36



# Setting mass of the virtual-site 0
for j in range(len(c)):
    if j == 0:
        group = row_groups[j]
        idx = group.atoms.indices
        gr = np.arange(idx[2], idx[-1]+1, 3)
        for k in gr:
            u.atoms[k].mass = 0
    elif (j != len(c) - 1) and (j % 2 != 0) and (j != 0):
        group = row_groups[j]
        gr = np.arange(1, len(group), 3)

        for k in gr:
            u.atoms[group[k].index].mass = 0
    elif (j != len(c) - 1) and (j % 2 == 0) and (j != 0):
        group = row_groups[j]
        gr = np.arange(2, len(group), 3)

        for k in gr:
            u.atoms[group[k].index].mass = 0

    elif j == len(c) - 1:
        group = row_groups[j]
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


def hexagon():
    idx_top = row_groups[0]
    idx_bottom = row_groups[-1]

    sel = u.atoms[(u.atoms.masses == 0)]
    sel1 = sel.atoms[sel.atoms.positions[:, 0]
                     == np.max(sel.atoms.positions[:, 0])]
    idx = u.atoms - idx_top - idx_bottom - sel1

    b = idx.atoms[idx.atoms.masses == 0].indices

    hexagon_indices = collect_hexagons_from_virtual_sites(u, b)

    '''
    Hexagons along the horizontal, or the x-axis
    '''
    upper = row_groups[0]
    lower = row_groups[-1]
    upper_down = row_groups[1]
    lower_up = row_groups[-2]

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
    right = [group[-2].index if i % 2 == 0 else group[-1].index for i, group in enumerate(row_groups)]
    left = [group[0].index for group in row_groups]
    vss = [row_groups[i][-1].index for i in range(len(c)) if i % 2 == 0]
    vss_use = vss[1:]
    right_groups = regroup_side(right)
    left_groups = regroup_side(left)
    for i in range(len(right_groups)):
        hexagon_indices.append(
            [
                right_groups[i][0],
                left_groups[i][0],
                right_groups[i][1],
                left_groups[i][1],
                right_groups[i][2],
                left_groups[i][2],
                vss_use[i],
            ]
        )
        
    '''
    Last hexagon to be added for the remaining virtual site
    '''
    hexagon_indices.append(
        [
            row_groups[-1][-1].index,
            row_groups[-1][0].index,
            row_groups[0][-2].index,
            row_groups[0][0].index,
            row_groups[1][-1].index,
            row_groups[1][0].index,
            row_groups[0][-1].index,
        ]
    )
            
    return hexagon_indices
    


vs = hexagon() # list of virtual sites, and the hexagon indices
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

real_particle_index = u.atoms[u.atoms.masses !=0].indices
for i in real_particle_index:
    u.atoms[i].mass = mass_of_each_real_particle 







def bonds(hexagons):
    return unique_bonds_from_template(hexagons, BOND_TEMPLATE)

def add_bonds(hexagons):
    return unique_bonds_from_template(hexagons, SHORT_BOND_TEMPLATE)


def long_bond(hexagons):
    return unique_bonds_from_template(hexagons, LONG_BOND_TEMPLATE)

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


def get_neighbors_angles(hexagons):
    side_to_hexagons = {}
    for hexagon_index, hexa in enumerate(hexagons):
        for i, j in EDGE_TEMPLATE:
            a, b = hexa[i], hexa[j]
            key = (a, b) if a < b else (b, a)
            side_to_hexagons.setdefault(key, []).append(hexagon_index)

    neighbor_pairs = []
    for side_key, owners in side_to_hexagons.items():
        if len(owners) == 2:
            i, j = owners
            if i > j:
                i, j = j, i
            neighbor_pairs.append((i, j, side_key))
    neighbor_pairs.sort()

    angles = []
    for i, j, side_key in neighbor_pairs:
        hexa1 = hexagons[i]
        hexa2 = hexagons[j]
        atom_a, atom_b = side_key

        if hexa1.index(atom_a) < hexa1.index(atom_b):
            common_side = [atom_a, atom_b]
        else:
            common_side = [atom_b, atom_a]

        hexa1_idx = [hexa1.index(common_side[0]), hexa1.index(common_side[1])]
        angle = None
        if hexa1_idx == [0, 1]:  # Top neighbor
            hexa2_idx = hexa2.index(common_side[0])
            angle = [hexa2[hexa2_idx - 4], common_side[0], common_side[1], hexa1[hexa1_idx[1] + 4]]
        elif hexa1_idx == [4, 5]:  # Bottom neighbor
            hexa2_idx = hexa2.index(common_side[1])
            angle = [hexa1[hexa1_idx[0] - 4], common_side[0], common_side[1], hexa2[hexa2_idx + 4]]
        elif hexa1_idx == [0, 2] or hexa1_idx == [1, 3]:  # Top right/left neighbor
            hexa2_idx = hexa2.index(common_side[1])
            angle = [hexa1[hexa1_idx[0] + 4], common_side[0], common_side[1], hexa2[hexa2_idx - 4]]
        elif hexa1_idx == [3, 5] or hexa1_idx == [2, 4]:  # Bottom left/right neighbor
            hexa2_idx = hexa2.index(common_side[0])
            angle = [hexa1[hexa1_idx[1] - 4], common_side[1], common_side[0], hexa2[hexa2_idx + 4]]

        if angle is not None:
            angles.append(angle)

    return angles

impropers = np.vstack(
    (get_angles(hexagons), get_neighbors_angles(hexagons)))
impropers = impropers + 1
profiler.mark("built topology terms")




#---------------#
# Topology File #
#---------------#


with open(tmp_itp_path, "w") as topology_file:
    topology_file.write(
        "; \n; Graphene topology\n; for the Martini3 force field\n;\n"
    )
    topology_file.write("; Roshan Shrestha\n; CNRS\n;\n\n")

    topology_file.write("[ moleculetype ]\n")
    topology_file.write("; molname	 nrexcl\n")
    topology_file.write("  GRA           1")

    topology_file.write("\n[ atoms ]\n")
    topology_file.write("; nr	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n")
    for i in range(1, u.atoms.n_atoms + 1):
        topology_file.write(
            f"  {i:<5}     TC5     0     GRA     B{i:<5}     {i:<5}     0     {u.atoms[i-1].mass}\n"
        )

    topology_file.write("\n[ bonds ]\n")
    topology_file.write("; i         j        funct  length  kb\n")
    for i in bonds(hexagons):
        topology_file.write(
            f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.24595     30000\n"
        )

    topology_file.write("; short_bonds_across\n")
    for i in add_bonds(hexagons):
        topology_file.write(
            f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.4260     30000\n"
        )

    topology_file.write("; long_bonds_across\n")
    for i in long_bond(hexagons):
        topology_file.write(
            f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.4919     30000\n"
        )

    topology_file.write("\n[ angles ]\n")
    topology_file.write("; i	 j	 k	 funct	 angle	 force_k\n")
    for i in find_angles(hexagons):
        topology_file.write(
            f"  {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     1    120     300\n"
        )

    topology_file.write("\n[ dihedrals ]\n")
    topology_file.write("; i	 j	 k	 l     funct	 ref.angle     force_k\n")
    for i in impropers:
        topology_file.write(
            f"  {i[0]:<3}     {i[1]:<3}     {i[2]:<3}     {i[3]:<3}    2     180     200\n"
        )

    topology_file.write("\n[ virtual_sitesn ]\n")
    topology_file.write("; site	 funct	 constructing atom indices\n")
    for i in v_site:
        topology_file.write(
            f"  {i[0]:<3}     1     {i[1]:<3}     {i[2]:<3}     {i[3]:<3}    {i[4]:<3}\n"
        )

    topology_file.write("\n[ exclusions ]\n")
    for i in vs:
        topology_file.write(
            f"{i[-1]+1:<3}     {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     {i[3]+1:<3}     {i[4]+1:<3}     {i[5]+1:<3}\n"
        )

profiler.mark("wrote temporary .itp")
os.replace(tmp_gro_path, gro_path)
os.replace(tmp_itp_path, itp_path)
profiler.mark("finalized outputs")
