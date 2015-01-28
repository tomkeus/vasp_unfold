#! /bin/bash

# This scripts packs the Python source files 
# into an executable ZIP archive

# Name for the executable
VASP_UNFOLD="vasp_unfold"

# Command used to invoke Python
PYTHON_COMMAND="python"

# Command used to invoke env
ENV_COMMAND="/usr/bin/env"


SRC_FILES="__main__.py parse.py unfolding.py utils.py write.py"

# Change into source directory
cd src

# Zip source files
zip ../${VASP_UNFOLD}.zip ${SRC_FILES}

cd ..

# Prepend shebang to the zip archive
echo "#! ${ENV_COMMAND} ${PYTHON_COMMAND}" | cat - ${VASP_UNFOLD}.zip > ${VASP_UNFOLD}

# Remove zip archive
rm ${VASP_UNFOLD}.zip

# Make it executable
chmod +x ${VASP_UNFOLD}
