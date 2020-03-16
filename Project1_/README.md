# Variational Monte Carlo
The structure of the program is based on the example structure given by the teaching staff in FYS4411. See: https://github.com/Schoyen/variational-monte-carlo-fys4411 . 

In this project we study particles trapped in different potentials. We start by looking at a non-interacting system, where the there is no forces between the particles in the system. After that we will introduce a new wavefunction where the interactions between the particles are considered. 

## Compilling and running the project
There are now several options you can use for compiling the project. If you use QT Creator, you can import this project into the IDE and point it to the `.pro`-file. If not, you can use CMake to create a Makefile for you which you can then run. You can install CMake through one of the Linux package managers, e.g., `apt install cmake`, `pacman -S cmake`, etc. For Mac you can install using `brew install cmake`. Other ways of installing are shown here: [https://cmake.org/install/](https://cmake.org/install/).

### Compiling the project using CMake
In a Linux/Mac terminal this can be done by the following commands
```bash
# Create build-directory
mkdir build

# Move into the build-directory
cd build

# Run CMake to create a Makefile
cmake ../

# Make the Makefile using two threads
make -j2

# Move the executable to the top-directory
mv vmc ..
```
Or, simply run the script `compile_project` via
```bash
./compile_project
```
and the same set of commands are done for you. Now the project can be run by executing
```bash
./vmc
```
in the top-directory.

#### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

#### Windows
Compilation of the project using Windows is still an open question to me, but please include a pull-request if you've got an example. CMake should be OS-independent, but `make` does not work on Windows.

