# CP2K data files for HQmdstemkit

Put these files from the CP2K repository into examples/cp2k/data/ :

  BASIS_MOLOPT_UZH   (basis sets, includes DZVP-MOLOPT-SR-GTH-q11/q12 for Cu/Zn)
  POTENTIAL_UZH      (GTH pseudopotentials: GTH-LDA-q11/q12 and GTH-PBE-q11/q12)
  dftd3.dat          (DFT-D3 parameters, needed only for PBE-D3)

Download URLs:
  https://github.com/cp2k/cp2k/raw/master/data/BASIS_MOLOPT_UZH
  https://github.com/cp2k/cp2k/raw/master/data/POTENTIAL_UZH
  https://github.com/cp2k/cp2k/raw/master/data/dftd3.dat

On SAI you can also use the module-installed data directory, e.g.:
  module load cp2k/2025.1-cuda12.4-sm70-auto
  echo $CP2K_DATA_DIR
then pass it to the workflow with:
  HQmdstemkit.sh 3 312 SRC --func PBE-D3 --data $CP2K_DATA_DIR

Functional choices:
  LDA    -> GTH-LDA potentials, LDA (PADE)
  PBE    -> GTH-PBE potentials, PBE
  PBE-D3 -> GTH-PBE potentials, PBE + DFT-D3 (dftd3.dat)