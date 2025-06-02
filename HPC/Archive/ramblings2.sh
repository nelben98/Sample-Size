
mkdir -p "$HOME/ICTU/INLA/R/config"
echo '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))' > "$HOME/ICTU/INLA/R/config/.Rprofile"

apptainer exec \
  --bind "$HOME/ICTU/INLA/R/config/.Rprofile:/home/nb221/.Rprofile" \
  --bind "/rds/general/user/nb221/home/anaconda3/lib/R/library:/rds/general/user/nb221/home/anaconda3/lib/R/library" \
  Panther.sif \
  bash -c 'export R_LIBS_USER=/rds/general/user/nb221/home/anaconda3/lib/R/library && R'