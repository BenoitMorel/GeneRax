mkdir -p scalasca_build
cd scalasca_build
cmake ..

filesToReplace=( "src/JointSearch/CMakeFiles/JointSearch.dir/build.make" "src/JointSearch/CMakeFiles/JointSearch.dir/link.txt" "src/core/CMakeFiles/jointsearch-core.dir/build.make" "src/core/CMakeFiles/jointsearch-core.dir/link.txt" "src/generax/CMakeFiles/generax.dir/link.txt" "src/generax/CMakeFiles/generax.dir/build.make")

gccPrefix="scorep "
gccSuffix="-g "

for fileToReplace in "${filesToReplace[@]}"
do
  sed -i "s#/opt/rh/devtoolset-6/root/usr/bin/c++#$gccPrefix /opt/rh/devtoolset-6/root/usr/bin/c++ $gccSuffix #g" $fileToReplace
done



