#!/bin/bash

usage() {
  echo "Usage: ./setup.sh [-v=VERSION] [-h]"
  echo "  -v=VERSION     setup VERSION version of scdb. If omitted,"
  echo "                 the current version will be setup."
  echo "  -h,  --help    print this help."
  echo
}

if (( $# == 0 )); then
  export PYTHONPATH="/data/single_cell_database/src/current:$PYTHONPATH"
elif (( $# == 1 )); then
  arg=$1
  if [[ ${arg:0:2} == "-v" ]]; then
    arg=${arg:2:${#arg}}
    if [[ ${arg:0:1} == "=" ]]; then
      arg=${arg:1:${#arg}}
    fi
    if (( ${#arg} == 0 )); then
      echo "ERROR: No version was passed!"
    else
      if [[ $arg != "current" ]]; then
        arg="v$arg"
      fi
      if [[ -d "/data/single_cell_database/src/$arg" ]]; then
        export PYTHONPATH="/data/single_cell_database/src/$arg:$PYTHONPATH"
      else
        echo "ERROR: $arg is not a valid version!"
      fi
    fi
  elif [[ ${arg:0:2} == "-h" ]]; then
    usage
  else
    echo "ERROR: Invalid argument!"
    echo
    usage
  fi
else
  echo "ERROR: Only 0 or 1 argument allowed!"
  echo
  usage
fi
