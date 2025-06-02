

.libPaths( c( "/rds/general/user/nb221/home/anaconda3/lib/R/library" , .libPaths() ) )
"/usr/local/lib/R/site-library" 
"/usr/local/lib/R/library"  


Sys.setenv(R_HOME="/rds/general/user/nb221/home/anaconda3/lib/R")
Sys.setenv(R_INCLUDE_DIR="/rds/general/user/nb221/home/anaconda3/lib/R/include")
Sys.setenv(R_LIBS="/rds/general/user/nb221/home/anaconda3/lib/R/site-library:/rds/general/user/nb221/home/anaconda3/lib/R/library")
Sys.setenv(R_LIBS_SITE="/rds/general/user/nb221/home/anaconda3/R/site-library")




#!/bin/bash

# Define personal R library location on the host
LIBDIR="rds/general/user/nb221/home/anaconda3/lib/R/library"
mkdir -p "$LIBDIR"

# Path to your Apptainer container
CONTAINER="Panther.sif"

# Define matching path inside container
CONTAINER_LIB="/home/nb221/R/library"

# Check if the container file exists
if [ ! -f "$CONTAINER" ]; then
  echo "Error: Container '$CONTAINER' not found."
  exit 1
fi

# Run R inside the container using the personal library path
apptainer exec \
  --bind "$LIBDIR:$CONTAINER_LIB" \
  --env R_LIBS_USER="$CONTAINER_LIB" \
  "$CONTAINER" R


# this provides a temporary solution
  apptainer exec --bind "/rds/general/user/nb221/home/anaconda3/lib/R/library:/rds/general/user/nb221/home/anaconda3/lib/R/library" Panther.sif bash -c 'export R_LIBS_USER=/rds/general/user/nb221/home/anaconda3/lib/R/library && R'


  apptainer exec  --bind "/rds/general/user/nb221/home/anaconda3/lib/R/config/.Rprofile:/home/nb221/.Rprofile"   --bind "/rds/general/user/nb221/home/anaconda3/lib/R/library:/home/nb221/R/library"   Panther.sif R

  apptainer exec --bind "$LIBDIR:/home/nb221/R/library"  --bind "/rds/general/user/nb221/home/anaconda3/lib/R/config/.Rprofile:/home/nb221/.Rprofile"  --env R_LIBS_USER="/home/nb221/R/library" Panther.sif R