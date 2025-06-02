# Open a wsl session and run the following:
# this is where the files are saved
# \\wsl.localhost\<distro name>
shortcut:
# \\wsl$\<Distribution>:
ex:
\\wsl$\Ubuntu


# this will copy the inla_rocker.def file to the wsl session so it can be run
cp /mnt/c/Users/nb221/'OneDrive - Imperial College London'/Documents/Panther_NEL/PANTHER_personal/'Sample Size'/HPC/inla_rocker.def ~/

# call apptainer and build the inla.sif (inmage) from the definition file inla_rocker.def
sudo apptainer build inla.sif inla_rocker.def

# then test by running the apptainer and seeing if the packages are working in r

apptainer shell inla.sif

R

library(INLA)

# If it work then save in the hpc with the following command:

scp ~/inla.sif nb221@login.hpc.ic.ac.uk:/rds/general/user/nb221/home/ICTU/INLA/