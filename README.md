# ABAQUSnonlocal
Subroutines to allow nonlocal simulations in ABAQUS

Contains the UEXTERNALDB and URDFIL subroutines for Abaqus to read in quantities from all Gauss points. 
Also contains nlmodule.f which can be called from the UMAT to access these quantities and an example of a python script.
Does not include the corresponding UMAT and other subroutines.
