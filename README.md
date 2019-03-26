# Caribbean
NEMO V4.0 regional configuration of the Caribbean

The following code was used in this configuration:

svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0

The initial conditions and boundary data can be downloaded from JASMIN:

http://gws-access.ceda.ac.uk/public/recicle/Caribbean/

NB This recipe has be written with the ARCHER HPC INTEL environment in mind.

```
# Change to some working directory of choice
export WORK_DIR='path_to_working_directory'
if [ ! -d "$WORK_DIR" ]; then
  mkdir $WORK_DIR
fi
cd $WORK_DIR

# Checkout the NEMO code from the SVN Paris repository 
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 nemo
cd nemo/cfgs

# Checkout configuration directory structure
git init .
git clone git@github.com:NOC-MSM/Caribbean.git
cd Caribbean
git checkout release-4.0
cd ../

# Add it to the configuration list
echo "Caribbean OCE" >> work_cfgs.txt
```

At this point you can checkout and compile XIOS or use a version you already have. If you're starting from scratch:

```
# Choose an appropriate directory for your XIOS installation
export XIOS_DIR='path_to_checkout_xios'
if [ ! -d "$XIOS_DIR" ]; then
  mkdir $XIOS_DIR
fi
cd $XIOS_DIR
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 xios
cd xios
cp $WORK_DIR/nemo/NEMOGCM/CONFIG/Caribbean/arch_xios/* ./arch
./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par --job 4

# Let's update the path to xios
export XIOS_DIR=$XIOS_DIR/xios
```

You can fold the ```make_xios``` command into a serial job. NB ```$NETCDF_DIR``` and ```$HDF5_DIR``` must be part of your environment. This should be the case if you've used ```modules``` to setup the netcdf and hdf5 e.g. 

```
module swap PrgEnv-cray PrgEnv-intel
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
```

Next, compile the NEMO code itself. First we copy the arch files into the appropriate directory.

```
cd $WORK_DIR/nemo/NEMOGCM/CONFIG/Caribbean
cp ARCH/* ../../ARCH
```

NB while ```$XIOS_DIR``` is in the current environment if you ever compile in a new session ```$XIOS_DIR``` will have to be redefined as ```../ARCH/arch-XC_ARCHER_Intel.fcm``` use this environment variable.

```
cd ../
./makenemo -n Caribbean -m XC_ARCHER_Intel -j 4
```

That should be enough to produce a valid executable. Now to copy the forcing data from JASMIN. 

```
cd Caribbean/EXP00
wget -r -np -nH --cut-dirs=3 -erobots=off --reject="index.html*" http://gws-access.ceda.ac.uk/public/recicle/Caribbean/
```

And finally link the XIOS binary to the configuration directory and create a restarts directory.

```
ln -s $XIOS_DIR/bin/xios_server.exe xios_server.exe
mkdir restarts
```

Edit and run the ```run_script.pbs``` script in ```../EXP00``` accordingly (namely enter a valid project code) and submit to the queue: ```qsub run_script.pbs```
