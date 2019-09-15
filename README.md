# Weighted_contact_map
Code to obtain all atom weighted contacts.
This folder contains 3 files.
README : all the instructions for using the code.
Makefile : This has the recipe to generate the executable.
AA_CA_weights.c : This is the C program which calcualtes the contacts and weights C alpha contacts according to all atom contacts.
The input is : a pdb file with mutliple frames / chains . Each frame/chain must be separated by a TER. The pdb file is expected to have an END as its last line.
The .pdb should not contain Hydrogen atoms.
The output is : AA_frame/chain_no.cont, this has the all atom contacts ; CA_frame/chain_no.cont this has the C alpha contacts and the atomic weights normalised to C alpha contacts.
