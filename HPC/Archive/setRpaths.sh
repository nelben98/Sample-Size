#!/bin/bash

# Define your host-side personal R lib path
HOST_LIB="$HOME/ICTU/R/library"
mkdir -p "$HOST_LIB"

# Define the container-side mount point (check your container user first)
CONTAINER_USER_HOME="/home/nb221"
CONTAINER_LIB="$CONTAINER_USER_HOME/R/library"

apptainer exec \
  --bind "$HOST_LIB":"$CONTAINER_LIB" \
  --env R_LIBS_USER="$CONTAINER_LIB" \
  Panther.sif bash -c "R"