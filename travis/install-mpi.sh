#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Julia

set -e
set -x

os=`uname`
TRAVIS_ROOT="$1"
MPI_IMPL="$2"

# this is where updated Autotools will be for Linux
export PATH=$TRAVIS_ROOT/bin:$PATH

case "$os" in
    Darwin)
        echo "Mac"
        brew update
        case "$MPI_IMPL" in
            mpich)
                brew install gcc || brew upgrade gcc || true
                brew link --overwrite gcc || true
                brew install mpich || true
                ;;
            openmpi)
                brew info open-mpi || true
                brew install open-mpi || true
                ;;
            *)
                echo "Unknown MPI implementation: $MPI_IMPL"
                exit 10
                ;;
        esac
    ;;

    Linux)
        echo "Linux"
        case "$MPI_IMPL" in
            mpich)
                if [ ! -d "$TRAVIS_ROOT/mpich" ]; then
                    VERSION=3.3b3
                    wget --no-check-certificate http://www.mpich.org/static/downloads/$VERSION/mpich-$VERSION.tar.gz
                    tar -xzf mpich-$VERSION.tar.gz
                    cd mpich-3*
                    mkdir build && cd build
                    ../configure CFLAGS="-w" --prefix=$TRAVIS_ROOT/mpich --disable-fortran --disable-static
                    make -j2
                    make install
                else
                    echo "MPICH already installed"
                fi
                ;;
            openmpi)
                if [ ! -d "$TRAVIS_ROOT/open-mpi" ]; then
                    VERSION=3.1.2
                    wget --no-check-certificate https://www.open-mpi.org/software/ompi/v3.1/downloads/openmpi-$VERSION.tar.gz
                    tar -xzf openmpi-$VERSION.tar.gz
                    cd openmpi-$VERSION
                    mkdir build && cd build
                    ../configure CFLAGS="-w" --prefix=$TRAVIS_ROOT/open-mpi \
                                --without-verbs --without-fca --without-mxm --without-ucx \
                                --without-portals4 --without-psm --without-psm2 \
                                --without-libfabric --without-usnic \
                                --without-udreg --without-ugni --without-xpmem \
                                --without-alps --without-munge \
                                --without-sge --without-loadleveler --without-tm \
                                --without-lsf --without-slurm \
                                --without-pvfs2 --without-plfs \
                                --without-cuda --disable-oshmem \
                                --disable-mpi-fortran --disable-oshmem-fortran \
                                --disable-libompitrace \
                                --disable-mpi-io  --disable-io-romio \
                                --disable-static #--enable-mpi-thread-multiple
                    make -j2
                    make install
                else
                    echo "Open-MPI already installed"
                fi
                ;;
            *)
                echo "Unknown MPI implementation: $MPI_IMPL"
                exit 20
                ;;
        esac
        ;;
esac
