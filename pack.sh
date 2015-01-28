#! /bin/bash

VASP_UNFOLD="vasp_unfold"
PYTHON_COMMAND="python"
ENV_COMMAND="/usr/bin/env"

SRC_FILES="__main__.py parse.py unfolding.py utils.py write.py"

zip ${VASP_UNFOLD}.zip ${SRC_FILES}

echo "#! ${ENV_COMMAND} ${PYTHON_COMMAND}" | cat - ${VASP_UNFOLD}.zip > ${VASP_UNFOLD}

rm ${VASP_UNFOLD}.zip

chmod +x ${VASP_UNFOLD}
