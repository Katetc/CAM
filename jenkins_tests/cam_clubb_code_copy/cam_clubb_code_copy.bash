#!/bin/bash

# $1 = clubb branch $2 = cam directory

if [ "$#" -eq 2 ]; then
	clubb_branch=$1
	cam_dir=$2
	echo "Running locally with the CLUBB branch: $clubb_branch and the CAM directory: $cam_dir."
elif [ $1 = "-h" ]; then
	echo "Usage: ./camClubbCodeCopy.bash [CLUBB branch] [CAM directory]"
	echo "For help type -h."
	exit 0
else 
	echo "Invalid arguments, please rerun with the option -h to see how the script is used."
	exit 1
fi

git clone https://github.com/larson-group/clubb.git clubb_copy_code
cd clubb_copy_code/
git checkout $clubb_branch
cd ..

cp clubb_copy_code/src/CLUBB_core/* $cam_dir/src/physics/clubb
cp clubb_copy_code/src/SILHS/* $cam_dir/src/physics/silhs

rm -rf clubb_copy_code/

echo "CLUBB files successfully updated!"
exit 0
