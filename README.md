# Void Abundance
funcionRadio.c   for voids abundance

funcionMasa_2019.c   for halo abundance


# INPUTS
voids catalog:
ASCII format x, y, z, r in Mpc/h


halos catalog:
ASCII format x, y, z, m, r: xyzr in Mpc/h; m is the halo mass in 10^11 m_\sun units

# OUTPUTS
dnv_... voids abundance

dnh...  haloes abundance

# Parameters
parameters.h

includes cosmological parameters, simulation box size, voids and haloes catalog file name, number of bins and limits for splitting the halos or voids
