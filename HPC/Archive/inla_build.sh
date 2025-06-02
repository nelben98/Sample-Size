#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=00:01:00
#PBS -N Build_image
#PBS -o /rds/general/user/nb221/home/ICTU/binaries
#PBS -e /rds/general/user/nb221/home/ICTU/binaries

# Load Apptainer module (name may vary depending on your HPC)
module purge
module load apptainer
cd $PBS_O_WORKDIR

# Path to your INLA definition file
DEF_FILE=inla.def
SIF_FILE=inla.sif

# Check if fakeroot is supported
echo "Checking fakeroot support..."
apptainer build --help | grep fakeroot

# Build the container using fakeroot
echo "Starting build at $(date)"
apptainer build --fakeroot $SIF_FILE $DEF_FILE
echo "Build completed at $(date)"