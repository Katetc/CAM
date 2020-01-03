#!/bin/bash

### Function definitions
myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This script creates pdf documents from AMWG diagnostics for all subdirectories of"
echo "AMWG_parent_directory, using CAM_amwgplots2pdf.bash."
echo ""
echo "Dependencies: pdflatex"
echo "" 
echo "Usage: ./CAM_amwgplots2pdf_all.bash [options] AMWG_parent_directory"
echo ""
echo "Options:"
echo "-f -- force overwrite"
echo "-h -- show this page"
echo "------------------------------------------------------------------------------------"
}


### Parse options
force=""

while getopts hf opt
do
  case "$opt" in
    f) force='-f'; shift;;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option."; myhelp; exit;;
  esac
done


parent_dir=$1
parent_dir=${parent_dir%/}
log_file=${parent_dir}/amwgplots2pdf.log

script_dir='.'
pdfscript=${script_dir}/CAM_amwgplots2pdf.bash

if [ "${parent_dir}" = "" ]; then
   echo "Error: You must provide the path to the directory, containing the AMWG diagnostics"
   echo "subdirectories!"
   exit
fi

if [ -f ${log_file} ]; then
   rm ${log_file}
fi

touch ${log_file}

dirlist=`ls -d ${parent_dir}/*/`

for dir in ${dirlist[@]}; do
    #NOTE: the -f option will overwrite existing pdf's without asking
    $pdfscript $force $dir >> ${log_file} 
done
