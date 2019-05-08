************
# Caribbean
************

This model configuration has been developed through the [Commonwealth Marine Economies Programme](http://projects.noc.ac.uk/cme-programme/)

"Enabling safe and sustainable marine economies across Commonwealth Small Island Developing States". 

The Commonwealth Marine Economies (CME) Programme was announced by the British Prime Minister in 2015 to help Commonwealth Small Island Developing States (SIDS) make the most of their natural maritime advantages, to enable sustainable economic growth and alleviate poverty.

********************************************
## NEMO regional configuration of the Caribbean
********************************************

The following code was used in this configuration:

svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk -r 8395

NB This recipe has be written with the [ARCHER](https://www.archer.ac.uk) HPC INTEL environment in mind.

```
# Change to some working directory of choice
export WORK_DIR='path_to_working_directory'
if [ ! -d "$WORK_DIR" ]; then
  mkdir $WORK_DIR
fi
cd $WORK_DIR

# Checkout the NEMO code from the SVN Paris repository 
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk -r 8395 nemo
cd nemo/NEMOGCM/CONFIG

# Checkout configuration directory structure
git init .
git clone git@github.com:NOC-MSM/Caribbean.git

# Add it to the configuration list
echo "Caribbean OPA_SRC" >> cfg.txt
```

You can fold the ```make_xios``` command into a serial job. NB ```$NETCDF_DIR``` and ```$HDF5_DIR``` must be part of your environment. This should be the case if you've used ```modules``` to setup the netcdf and hdf5 e.g. 

```
module swap PrgEnv-cray PrgEnv-intel
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
```

At this point you can checkout and compile XIOS or use a version you already have. If you're starting from scratch:

```
# Choose an appropriate directory for your XIOS installation
export XIOS_DIR='path_to_checkout_xios'
if [ ! -d "$XIOS_DIR" ]; then
  mkdir $XIOS_DIR
fi
cd $XIOS_DIR
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk@1242 xios
cd xios
cp $WORK_DIR/nemo/NEMOGCM/CONFIG/Caribbean/arch_xios/* ./arch
./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par --job 4

# Let's update the path to xios
export XIOS_DIR=$XIOS_DIR/xios
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

That should be enough to produce a valid executable. Now to extract the forcing data (copy the attached tar ball to the EXP00 directory): 

```
cd Caribbean/EXP00
tar xvfz caribbean_netcdf_files.tar.gz
```

Copy files generated using the ERA5 script (`SCRIPTS/ERA5_gen.py`) into `./SBC`:

```
cp path_to_era5_files SBC
```

And finally link the XIOS binary to the configuration directory and create a restarts directory.

```
ln -s $XIOS_DIR/bin/xios_server.exe xios_server.exe
mkdir restarts
```

Edit and run the ```run_script.pbs``` script in ```../EXP00``` accordingly (namely enter a valid project code) and submit to the queue: ```qsub run_script.pbs```
