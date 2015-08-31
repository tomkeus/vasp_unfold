#! /bin/bash

# This scripts packs the Python source files 
# into an executable ZIP archive. This is done
# separatly for the vasp_unfold and 

# Names for the executables
VASP_UNFOLD="vasp_unfold"
PLOT="fatplot"

# Command used to invoke Python
PYTHON_COMMAND="python"

# Command used to invoke env
ENV_COMMAND="/usr/bin/env"


SRC_FILES="__main__.py parse.py unfolding.py utils.py write.py errors.py"
PLOT_SRC_FILES="__main__.py"

# Change into source directory
cd src

# Zip the source files
cd unfolding; zip ../../${VASP_UNFOLD}.zip ${SRC_FILES}; cd ..;
cd plot;      zip ../../${PLOT}.zip ${PLOT_SRC_FILES};   cd ..;

cd ..

# Prepend shebang to the zip archive
echo "#! ${ENV_COMMAND} ${PYTHON_COMMAND}" | cat - ${VASP_UNFOLD}.zip > ${VASP_UNFOLD}
echo "#! ${ENV_COMMAND} ${PYTHON_COMMAND}" | cat - ${PLOT}.zip > ${PLOT}

# Remove zip archive
rm ${VASP_UNFOLD}.zip
rm ${PLOT}.zip

# Make it executable
chmod +x ${VASP_UNFOLD}
chmod +x ${PLOT}
