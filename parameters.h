// cosmological parameters
const double Om=0.267;
const double Ol=0.733;
const double H0=71.9;

// simulation box size
const double Lbox = 256.;		//	longitud de la caja en [Mpc/h]
const int    np   = 512;		//	numero de particulas en cada dimension


char voids_file[]="esferas_512x256_03.txt";			//	voids catalog
const int   acumulada=0;		//	1:computes cumulated abundance, 0:computes usual abundance


const double radio_min_todos=1.155162;  //	minimum radius considered
const double radio_max_todos=18.28353;  //	maximum radius considered
const int   bin=15;                    //	radial log bins for splitting the voids sample



char halos_file[]="halos_512x256_03.txt";			//	voids catalog

const double masa_min_todos=9.;		//	minimum halo mass considered in 10^11 M_\sun
const double masa_max_todos=9000.;	//	maximum halo mass considered in 10^11 M_\sun
//const int   bin=30;			//	number of log radial bins for splitting the halo sample

