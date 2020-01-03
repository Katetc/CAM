#!/bin/bash

### Function definitions
myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This script creates a pdf document that contains a subset of plots of CESM's AMWG diagnostics."
echo "See https://www2.cesm.ucar.edu/working-groups/amwg/amwg-diagnostics-package."
echo "The plots can be chose via the variable imglist in this script."
echo ""
echo "Dependencies: pdflatex"
echo "" 
echo "Usage: ./CAM_amwgplots2pdf.bash AMWG_directory"
echo ""
echo "Options:"
echo "-f -- force overwrite"
echo "-h -- show this page"
echo "------------------------------------------------------------------------------------"
}


### Parse options
force=false

while getopts hf opt
do
  case "$opt" in
    f) force=true; shift;;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option."; myhelp; exit;;
  esac
done

### Init
amwg_dir=$1 # Input for this script; AMWG diagnostics directory
amwg_dir=${amwg_dir%/} # Remove last slash from path

### Define list of images to be included in the pdf
declare -a imglist
declare -a txtlist
imglist=( "set5_6/set5_ANN_SWCF_CERES-EBAF_obsc.png" \
          "set5_6/set5_ANN_LWCF_CERES-EBAF_obsc.png" \
          "set5_6/set5_ANN_PRECT_GPCP_obsc.png" \
          "set5_6/set5_ANN_PREH2O_NVAP_obsc.png" \
          "set5_6/set5_ANN_TGCLDLWP_NVAP_obsc.png" \
          "set5_6/set6_ANN_STRESS_ERS_obsc.png" \
          "set5_6/set5_ANN_PSL_MERRA_obsc.png" \
          "set5_6/set5_ANN_U_200_ERAI_obsc.png" \
          "set5_6/set5_ANN_SHFLX_LARYEA_obsc.png" \
          "set5_6/set5_ANN_QFLX_LARYEA_obsc.png" \
          "set5_6/set5_ANN_TREFHT_WILLMOTT_obsc.png" \
          "set14/set14_ANN_SPACE_TIME_obsc.png")

txtlist=( "set1/table_GLBL_ANN_obs.asc" )

echo ""
echo "================== CAM_amwgplots2pdf.bash ================="

### Check if amwg_dir is genuine
if [ "${amwg_dir}" == "" ]; then
   echo "Error: You must provide the path to the AMWG diagnostics!"
   exit
elif [[ -f "${amwg_dir}" ]]; then
   # The given directory might be a tar file. Try to extract it.
   echo "Extracting archive"
   tar -xf "${amwg_dir}" -C "`dirname "${amwg_dir}"`" || { echo "The path provided was not a tar file." ; exit 1; }
   # Set the AMWG directory: simply strip the ".tar" extension.
   amwg_dir="${amwg_dir%.*}"
elif [ ! -d ${amwg_dir} ]; then
   echo "Error: Could not find the AMWG diagnostics directory: ${amwg_dir}!"
   exit
fi

outfile_dir=${amwg_dir}/UWM_diag # Output directory (where the pdf will be located)
outfile_name="${amwg_dir##*/}_UWM_diag" # Name of the pdf document (without extension)

### Let's create the pdf
output_tex="${outfile_dir}/tex/${outfile_name}.tex"

# Check if output directory already exists, if so remove it
if [ -d ${outfile_dir} ]; then
   echo "WARNING: The output directory ${outfile_dir} already exists."

   if [ "$force" = "false" ]; then 
       # wait for user input
#      while true; do
#          read -p "Do you want to overwrite ${outfile_dir}? (y/n): " yn
#          case $yn in
#               [Yy]* ) rm -r ${outfile_dir}; break;;
#               [Nn]* ) echo "Exiting script.."; exit;;
#               * ) echo "Please answer yes (y) or no (n).";;
#          esac
#      done
       echo "Exiting CAM_amwgplots2pdf.bash because directory already exists . . ."
       exit 1
   else # force = true; overwrite directories without asking
      rm -r ${outfile_dir}
   fi

else
   echo "Creating output directory: ${outfile_dir}.."
fi

# Create directory structure
mkdir ${outfile_dir}
mkdir ${outfile_dir}/tex
touch ${output_tex}

# Create tex file 
echo "Creating tex file.."

echo "\documentclass{beamer}" > ${output_tex}
echo "\usepackage{graphicx}" >> ${output_tex}
echo "\usepackage{verbatim}" >> ${output_tex}
echo "\begin{document}" >> ${output_tex}

# Add images from imglist
for file in ${imglist[@]}; do

    # Check if files exist
    if [ -f ${amwg_dir}/$file ]; then # Include in tex document
       echo "   \begin{figure}" >> ${output_tex}
       echo "      \centering" >> ${output_tex}
       echo "      \includegraphics[height=0.95\textheight]{"${amwg_dir}/${file}"}" >> ${output_tex}
       echo "   \end{figure}" >> ${output_tex}
       echo "   \clearpage" >> ${output_tex}
       echo "" >> ${output_tex} 
    else # Print warning message and ignore file
       echo "WARNING: Can not find $file! The image will not be included in the PDF."
    fi

done

# Add tables at the end of the document
for file in ${txtlist[@]}; do

    # Check if files exist
    if [ -f ${amwg_dir}/$file ]; then # Include in tex document
       echo "   {\tiny" >> ${output_tex}
       echo "      \verbatiminput{"${amwg_dir}/${file}"}" >> ${output_tex}
       echo "      \clearpage" >> ${output_tex}
       echo "   }" >> ${output_tex}
       echo "" >> ${output_tex} 
    else # Print warning message and ignore file
       echo "WARNING: Can not find $file! The file will not be included in the PDF."
    fi

done

echo "\end{document}" >> ${output_tex}

### Create pdf from tex code
echo "Creating pdf document.."
pdflatex -output-directory ${outfile_dir} ${output_tex} > ${outfile_dir}/pdflatex.log

echo "==========================================================="
echo "If the script was successful, you will find your pdf here:"
echo ${outfile_dir}
echo "If you can not find the pdf in the directory above, you might want to have a look at:"
echo ${outfile_dir}/pdflatex.log
echo "==========================================================="
echo ""
