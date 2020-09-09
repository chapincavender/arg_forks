#!/usr/bin/env python3

import argparse
import loos
import numpy

# Command-line arguments
parser = argparse.ArgumentParser(description = 'Find nucleobase contacts.')
parser.add_argument(
    '-a', '--angle_cutoff',
    help = 'Cutoff for angle between the normal vectors of two planes in deg.',
    type = float,
    default = 30.0
)
parser.add_argument(
    '-b', '--backbone_suites',
    help = 'Suite definition file. If none is given, suites are not assigned.',
    type = str,
    default = None
)
parser.add_argument(
    '-d', '--distance_cutoff',
    help = 'Cutoff for H bond heavy atom distance in angstrom.',
    type = float,
    default = 3.5
)
parser.add_argument(
    '-p', '--protein',
    help = 'Selection string for protein atoms.',
    type = str,
    default = 'resname == "ARG"'
)
parser.add_argument(
    '-r', '--rna',
    help = 'Selection string for RNA atoms.',
    type = str,
    default = 'resname =~ "^[ACGU][35N]?$"'
)
parser.add_argument(
    '-s', '--stack_cutoff',
    help = 'Cutoff for distance between planar groups in angstrom.',
    type = float,
    default = 6.0
)
parser.add_argument(
    '-w', '--rotamer_width',
    help = 'Width of rotamer interval around common atom values in deg.',
    type = float,
    default = 35.0
)
parser.add_argument(
    'input_pdb',
    help = 'File containing PDB structure.',
    type = str
)
args = parser.parse_args()

# Check whether an arginine residue contains NE, NH1, and NH2 atoms
def is_valid_arg_residue(residue):
    return (len(loos.selectAtoms(residue, 'name == "NE"')) == 1
            and len(loos.selectAtoms(residue, 'name == "NH1"')) == 1
            and len(loos.selectAtoms(residue, 'name == "NH2"')) == 1)

# Check whether an RNA residue contains C2, C4, and C6 atoms
def is_valid_rna_residue(residue):
    return (len(loos.selectAtoms(residue, 'name == "C2"')) == 1
            and len(loos.selectAtoms(residue, 'name == "C4"')) == 1
            and len(loos.selectAtoms(residue, 'name == "C6"')) == 1)

# Calculate the square of the sine of the angle between two normalized vectors.
# sin() and sqrt() are both expensive, and sin() will be small for both parallel
# and anti-parallel vectors.
# Since the vectors are normalized, |A| = |B| = 1 and sin(theta) = |A cross B|
def sine_angle_sq(A, B):
    return numpy.sum(numpy.square(numpy.cross(A, B)))

# Wrap an angle into the range (-180, 180]
def wrap_180(angle):
    return angle - numpy.ceil(angle / 360.0 - 0.5) * 360.0

# Define Arg rotamer library (Lovell et al. (2000) PSFG, 40: 389-408)
# This is a dictionary with the key given by the character for the first three
# sidechain torsions and the value given by a tuple:
# (List of common atom angles for first three sidechain torsions,
#  List of possible common angle angles for the fourth sidechain torsion,
#  List of names associated with the fourth sidechain torsion)
arg_rotamer_library = {
    "ptp": ([  62,  180,   65], [  85, -175      ], [  "85",  "180"        ]),
    "ptt": ([  62,  180,  180], [  85,  180,  -85], [  "85",  "180",  "-85"]),
    "ptm": ([  62,  180,  -65], [ 175,  -85      ], [ "180",  "-85"        ]),
    "tpp": ([-177,   65,   65], [  85, -175      ], [  "85",  "180"        ]),
    "tpt": ([-177,   65,  180], [  85,  180      ], [  "85",  "180"        ]),
    "ttp": ([-177,  180,   65], [  85, -175, -105], [  "85",  "180", "-105"]),
    "ttt": ([-177,  180,  180], [  85,  180,  -85], [  "85",  "180",  "-85"]),
    "ttm": ([-177,  180,  -65], [ 105,  175,  -85], [ "105",  "180",  "-85"]),
    "mtp": ([ -67,  180,   65], [  85, -175, -105], [  "85",  "180", "-105"]),
    "mtt": ([ -67,  180,  180], [  85,  180,  -85], [  "85",  "180",  "-85"]),
    "mtm": ([ -67,  180,  -65], [ 105,  175,  -85], [ "105",  "180",  "-85"]),
    "mmt": ([ -62,  -68,  180], [  85,  180,  -85], [  "85",  "180",  "-85"]),
    "mmm": ([ -62,  -68,  -65], [ 175,  -85      ], [ "180",  "-85"        ])
}

# Assign Arg rotamer based on four sidechain torsions
def get_arg_rotamer(chi_1, chi_2, chi_3, chi_4, width = args.rotamer_width,
    library = arg_rotamer_library):

    # First part of rotamer is assigning first three torsions to gauche plus,
    # gauche minus, or trans
    rotamer = '!!!'
    for i in arg_rotamer_library.keys():

        if (numpy.abs(wrap_180(chi_1 - library[i][0][0])) < width
            and numpy.abs(wrap_180(chi_2 - library[i][0][1])) < width
            and numpy.abs(wrap_180(chi_3 - library[i][0][2])) < width):

            rotamer = i

    if rotamer not in library:
        return rotamer + "!"

    for i in range(len(library[rotamer][1])):
        if numpy.abs(wrap_180(chi_4 - library[rotamer][1][i])) < width:
            return rotamer + library[rotamer][2][i]

    # If no match was found for chi_4, this is an outlier
    return rotamer + "!!!!"

rad2deg = 180.0 / numpy.pi

# List of arginine hydrogen bond donors and associated hydrogens
arg_donors = ["NE", "NH1", "NH2"]
N_donor = len(arg_donors)

# Dictionary of RNA hydrogen bond acceptors by residue name
rna_acceptors = {
    "A": [ "N1",  "N3",  "N7", "O2'", "O3'", "O4'", "O5'", "OP1", "OP2"],
    "C": [ "N3",  "O2", "O2'", "O3'", "O4'", "O5'", "OP1", "OP2"],
    "G": [ "N3",  "N7",  "O6", "O2'", "O3'", "O4'", "O5'", "OP1", "OP2"],
    "U": [ "O2",  "O4", "O2'", "O3'", "O4'", "O5'", "OP1", "OP2"]
}
N_acceptor = {s: len(rna_acceptors[s]) for s in rna_acceptors}

# Compare to square of distance and angle cutoffs because sqrt() is expensive
# Compare to sin of angle cutoff because planar groups can have normal
# vectors that are close to parallel or close to anti-parallel
dist_cutoff_sq = numpy.square(args.distance_cutoff)
stack_cutoff_sq = numpy.square(args.stack_cutoff)
angle_cutoff_sq = numpy.square(numpy.sin(args.angle_cutoff / rad2deg))

# Create system model
system = loos.createSystem(args.input_pdb)

# Get arginine residues
arg_atoms = loos.selectAtoms(system, args.protein)
all_arg_residues = arg_atoms.splitByResidue()
arg_residues = [res for res in all_arg_residues if is_valid_arg_residue(res)]
N_arg = len(arg_residues)

# If there are no arginines, then exit
if N_arg == 0:
    quit()

# Get RNA residues
rna_atoms = loos.selectAtoms(system, args.rna)
all_rna_residues = rna_atoms.splitByResidue()
rna_residues = [res for res in all_rna_residues if is_valid_rna_residue(res)]
N_rna = len(rna_residues)

# Get resids and resnames for Arg and RNA
arg_chainids = [res[0].chainId() for res in arg_residues]
arg_resids = [res[0].resid() for res in arg_residues]
arg_resnames = [res[0].resname() for res in arg_residues]
rna_chainids = [res[0].chainId() for res in rna_residues]
rna_resids = [res[0].resid() for res in rna_residues]
rna_resnames = [res[0].resname() for res in rna_residues]

# Loop over Arg residues
arg_coords = numpy.zeros((N_arg, N_donor, 3))
arg_phi = []
arg_psi = []
arg_rotamer = []

for i in range(N_arg):

    res = arg_residues[i]

    # Construct arrays of coordinates for H bond donors: NE, NH1, & NH2
    for j in range(N_donor):

        coords = loos.selectAtoms(
            res, 'name == "%s"' % arg_donors[j])[0].coords()
        arg_coords[i][j] = numpy.array([coords.x(), coords.y(), coords.z()])

    # Get backbone torsions
    prev_co = loos.selectAtoms(
        system, 'chainid == "%s" && resid == %d && name == "%s"' % (
            res[0].chainId(), res[0].resid() - 1, "C"))
    n = loos.selectAtoms(res, 'name == "%s"' % "N")[0].coords()
    ca = loos.selectAtoms(res, 'name == "%s"' % "CA")[0].coords()
    co = loos.selectAtoms(res, 'name == "%s"' % "C")[0].coords()
    next_n = loos.selectAtoms(
        system, 'chainid == "%s" && resid == %d && name == "%s"' % (
            res[0].chainId(), res[0].resid() + 1, "N"))

    if len(prev_co) == 1:
        arg_phi.append('%11.6f' % loos.torsion(prev_co[0].coords(), n, ca, co))
    else:
        arg_phi.append("N")

    if len(next_n) == 1:
        arg_psi.append('%11.6f' % loos.torsion(n, ca, co, next_n[0].coords()))
    else:
        arg_psi.append("N")

    # Get sidechain rotamer
    cb = loos.selectAtoms(res, 'name == "%s"' % "CB")[0].coords()
    cg = loos.selectAtoms(res, 'name == "%s"' % "CG")[0].coords()
    cd = loos.selectAtoms(res, 'name == "%s"' % "CD")[0].coords()
    ne = loos.selectAtoms(res, 'name == "%s"' % "NE")[0].coords()
    cz = loos.selectAtoms(res, 'name == "%s"' % "CZ")[0].coords()

    chi_1 = loos.torsion(n, ca, cb, cg)
    chi_2 = loos.torsion(ca, cb, cg, cd)
    chi_3 = loos.torsion(cb, cg, cd, ne)
    chi_4 = loos.torsion(cg, cd, ne, cz)

    arg_rotamer.append(get_arg_rotamer(chi_1, chi_2, chi_3, chi_4))

# Construct local coordinate system for each Arg residue

# Origin of local coordinate system is centroid of NE, NH1, and NH2
arg_origin = numpy.sum(arg_coords, axis = 1) / 3.0

# x axis is opposite direction from origin to NE
arg_x = arg_origin - arg_coords[:, 0]

# z axis is cross product of the x axis with a coplanar vector to NH2
arg_z = numpy.cross(arg_x, arg_coords[:, 2] - arg_origin)

# y axis is cross product of z axis with x axis
arg_y = numpy.cross(arg_z, arg_x)

# Normalize unit vectors for protein local coordinate system
arg_x /= numpy.linalg.norm(arg_x, axis = 1)[:, None]
arg_y /= numpy.linalg.norm(arg_y, axis = 1)[:, None]
arg_z /= numpy.linalg.norm(arg_z, axis = 1)[:, None]

# Get boolean array of whether RNA residues are purines
purines = numpy.array([resname in ['A', 'G'] for resname in rna_resnames])

# Loop over RNA residues and construct arrays of coordinates for C2, C4, C6, and
# acceptor atoms and a boolean array of whether the residues are purines
c2 = numpy.zeros((N_rna, 3))
c4 = numpy.zeros((N_rna, 3))
c6 = numpy.zeros((N_rna, 3))
rna_coords = numpy.zeros((N_rna, numpy.max([i for i in N_acceptor.values()]), 3))

for i in range(N_rna):

    res = rna_residues[i]
    c2_coords = loos.selectAtoms(res, 'name == "C2"')[0].coords()
    c2[i] = numpy.array([c2_coords.x(), c2_coords.y(), c2_coords.z()])
    c4_coords = loos.selectAtoms(res, 'name == "C4"')[0].coords()
    c4[i] = numpy.array([c4_coords.x(), c4_coords.y(), c4_coords.z()])
    c6_coords = loos.selectAtoms(res, 'name == "C6"')[0].coords()
    c6[i] = numpy.array([c6_coords.x(), c6_coords.y(), c6_coords.z()])

    for j in range(N_acceptor[rna_resnames[i]]):

        atom = loos.selectAtoms(
            res, 'name == "%s"' % rna_acceptors[rna_resnames[i]][j])

        if len(atom) == 1:
            coords = atom[0].coords()
            rna_coords[i][j] = numpy.array([coords.x(), coords.y(), coords.z()])

# Construct local coordinate system for each RNA residue

# Origin of local coordinate system is centroid of C2, C4, and C6 atoms
rna_origin = (c2 + c4 + c6) / 3.0

# x axis is based on C2
rna_x = c2 - rna_origin

# z axis is based on C6 for purines and on C4 for pyrimidines
# z axis is cross product of the x axis with a coplanar vector to C4 or C6
rna_z = numpy.cross(
    rna_x, numpy.where(purines[:, None], c6 - rna_origin, c4 - rna_origin))

# y axis is cross product of z axis with x axis
rna_y = numpy.cross(rna_z, rna_x)

# Normalize unit vectors for RNA local coordinate system
rna_x /= numpy.linalg.norm(rna_x, axis = 1)[:, None]
rna_y /= numpy.linalg.norm(rna_y, axis = 1)[:, None]
rna_z /= numpy.linalg.norm(rna_z, axis = 1)[:, None]

# Calculate square of heavy atom distance for hydrogen bonds
h_bond_dist_sq = numpy.sum(
    numpy.square(
        arg_coords[:, :, None, None, :] - rna_coords[None, None, :, :, :]),
    axis = 4)

# Calculate square of distance between the origins of the local coordinate
# systems of the Arg-RNA planar groups
stack_dist_sq = numpy.sum(
    numpy.square(arg_origin[:, None, :] - rna_origin[None, :, :]), axis = 2)

# Calculate square of the sine of the angle between the Arg-RNA planar groups
angle_sq = numpy.array([[sine_angle_sq(a, b) for b in rna_z] for a in arg_z])

# Print header
if args.backbone_suites:
    print("# pdb_file donor_chain donor_resid donor_resname donor_atom "
        "acceptor_chain acceptor_resid acceptor_resid acceptor_atom "
        "h_bond_dist plane_angle arg_phi arg_psi arg_rotamer suite_name "
        "suiteness (stack_resname stack_resid stack_dist stack_angle)_N")
else:
 print("# pdb_file donor_chain donor_resid donor_resname donor_atom "
        "acceptor_chain acceptor_resid acceptor_resid acceptor_atom "
        "h_bond_dist plane_angle arg_phi arg_psi arg_rotamer "
        "(stack_resname stack_resid stack_dist stack_angle)_N")

# Loop over Arg residues
for i in range(N_arg):

    # Loop over H bond donors
    for ii in range(N_donor):

        # Loop over RNA residues
        for j in range(N_rna):

            # Loop over H bond acceptors
            for jj in range(N_acceptor[rna_resnames[j]]):

                if h_bond_dist_sq[i][ii][j][jj] < dist_cutoff_sq:

                    # Find RNA cation-pi stacks with this Arg
                    stack_indices = numpy.where(
                        numpy.logical_and(stack_dist_sq[i] < stack_cutoff_sq,
                                          angle_sq[i] < angle_cutoff_sq))[0]

                    # Print H bond atoms, distance, and plane angle
                    out = '%8s\t%1s\t%4d\t%3s\t%3s\t%1s\t%4d\t%3s\t%3s' % (
                        args.input_pdb, arg_chainids[i], arg_resids[i],
                        arg_resnames[i], arg_donors[ii], rna_chainids[j],
                        rna_resids[j], rna_resnames[j],
                        rna_acceptors[rna_resnames[j]][jj])
                    out += '\t%8.6f\t%9.6f' % (
                        numpy.sqrt(h_bond_dist_sq[i][ii][j][jj]),
                        numpy.arcsin(numpy.sqrt(angle_sq[i][j])) * rad2deg)

                    # Print Arg backbone torsions and sidechain rotamer
                    out += '\t%11s\t%11s\t%7s' % (
                        arg_phi[i], arg_psi[i], arg_rotamer[i])

                    # Print RNA backbone suites
                    if args.backbone_suites:

                        suite_name = "NN"
                        suiteness = 0.0
                        acceptor_name = rna_acceptors[rna_resnames[j]][jj]

                        if j != (N_rna - 1) and acceptor_name == "O3'":

                            suite = loos.RnaSuite(
                                rna_residues[j] + rna_residues[j + 1],
                                args.backbone_suites
                            )

                        elif j != 0 and acceptor_name in ["O5'", "OP1", "OP2"]:

                            suite = loos.RnaSuite(
                                rna_residues[j - 1] + rna_residues[j],
                                args.backbone_suites
                            )

                        suite.calculateBackboneDihedrals()
                        suite.assignSuitenameSuites()
                        suite_names = suite.getSuiteNames()

                        if len(suite_names) == 1:
                            suite_name = suite_names[0]
                            suiteness = suite.getSuitenessScores()[0]

                        out += '\t%2s\t%4.2f' % (
                            suite_name, suiteness)

                    # Print Arg-RNA stacks
                    for k in stack_indices:
                        if k == j:
                            continue
                        out += '\t%3s\t%4d\t%8.6f\t%9.6f' % (
                            rna_resnames[k], rna_resids[k],
                            numpy.sqrt(stack_dist_sq[i][k]),
                            numpy.arcsin(numpy.sqrt(angle_sq[i][k])) * rad2deg)

                    print(out)

