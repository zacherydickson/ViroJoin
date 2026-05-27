projectDir=$(dirname "$(dirname "$(readlink -f $0)")")
cd "$projectDir" || exit 1;
installDir="$projectDir/built_libs"
mkdir -p "$installDir"
htsLibVer=1.21
igraphVer=1.0.1
libleidenalgVer=0.12.0

#Build HTSLib
tar -xjf "htslib-$htsLibVer.tar.bz2"
cd "htslib-$htsLibVer"
./configure --prefix=$installDir
make
make install
cd ..
rm -rf "htslib-$htsLibVer"

##Build igraph
tar -xzf tar -xzf "igraph-$igraphVer.tar.gz"
cd "igraph-$igraphVer"
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX="$installDir" -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build .
cmake --install .
cd ../..
rm -rf "igraph-$igraphVer"

#Build libleidenalg
tar -xzf "libleidenalg-$libleidenalgVer.tar.gz" || exit 1;
cd "libleidenalg-$libleidenalgVer" || exit 1;
mkdir build && cd build || exit 1;
cmake .. -DCMAKE_INSTALL_PREFIX="$installDir" -DCMAKE_PREFIX_PATH="$installDir"
cmake --build .
cmake --install .
cd ../../
rm -rf "libleidenalg-$libleidenalgVer"
