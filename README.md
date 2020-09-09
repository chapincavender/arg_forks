# arg_forks
Analyze arginine-RNA interactions

# Dependencies
You must install LOOS (download at [https://github.com/GrossfieldLab/loos]) to 
use this script. See `INSTALL.md` for installation instructions and dependencies.

# Usage
Basic usage is

    get_arg_rna_contacts.py XXXX.pdb

To assign RNA backbone suites, you must specify the location of the suite 
definition file from LOOS (i.e. suitename_definitions.dat) with the `-b` flag.

    get_arg_rna_contacts.py -b $CONDA_PREFIX/share/loos/suitename_defintions.dat \
        XXXX.pdb

The `-d`, `-s`, `-a`, and `-w` flags change the cutoffs for the hydrogen bond 
distance, the stacking distance, the stacking angle, and the rotamer width. For 
example, the following changes the rotamer width to 30 deg.

    get_arg_rna_contacts.py -b $CONDA_PREFIX/share/loos/suitename_defintions.dat \
        -w 30 XXXX.pdb

By default, atoms in the proteins are any atom with residue name "ARG". Atoms in 
the RNA are any atom with residue name "A", "C", "G", or "U" possibly followed 
by exactly one of "3", "5", or "N". These defaults can be changed with the `-p` 
and `-r` flags, respectively.

# Output
Output is tab-delimited with one row per arginine-RNA hydrogen bond. The output 
columns are the PDB file followed by the arginine donor chain ID, residue ID, 
residue name, and atom name and then the RNA acceptor chain ID, residue ID, 
residue name, and atom name. Next are the heavy atom hydrogen bond distance in 
angstroms, the angle between the normal vectors of the guanidinium group and the 
nucleobase in deg, the protein backbone torsions (phi and psi) in deg, and the 
arginine rotamer. If RNA backbone suites were requested, the assigne suited and 
suiteness score are given next. Finally, for each RNA nucleobase that stacks on 
the guanidinium group, the output will contain the RNA residue name, RNA residue 
ID, stacking distance in angstroms, and angle between the normal vectors in deg.

For arginine rotamers and RNA backbone suites, a string containing an 
exclamation mark "!" indicates that the residue was assigned as an outlier. For 
RNA backbone suites, "NN" indicates that there was no valid dinucleotide to 
which to assign a backbone conformation because of a chain termination.
