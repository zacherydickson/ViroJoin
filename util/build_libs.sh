projectDir=$(dirname "$(dirname "$(readlink -f $0)")")
cd "$projectDir" || exit 1;
installDir="$projectDir/built_libs"
mkdir -p "$installDir"
htsLibVer=1.21

#Build HTSLib
tar -xjf "htslib-$htsLibVer.tar.bz2"
cd "htslib-$htsLibVer"
./configure --prefix=$installDir
make
make install
cd ..
rm -rf "htslib-$htsLibVer"
