#!/bin/bash
mtput() {
    if command -v tput > /dev/null; then
        tput $@
    fi
}

if uname -s | grep -Fqe CYGWIN ; then
    # TODO: Make CLI replacement
    echo "battlecode.sh won't work on windows! Use battlecode.bat :)"
    exit 1
fi
if uname -s | grep -Fqe MINGW ; then
    echo "battlecode.sh won't work on windows! Use battlecode.bat :)"
    exit 1
fi

# Use tput to show different colors in the terminal
mtput setaf 5
<<<<<<< HEAD
echo "$ pip3 install --user cffi tqdm werkzeug ujson psutil"
=======
>>>>>>> 6b0b8df56ef8ebdba88911a44d7374befbda1e30
mtput sgr0
pip3 install -q --user cffi tqdm werkzeug ujson psutil

RESULT=$?
if [ $RESULT -ne 0 ]; then
    echo "Warning: pip install failed!"
    echo "I'll keep going, but maybe try to fix whatever error you just got."
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH="$DIR/battlecode/python:$PYTHONPATH"
<<<<<<< HEAD
python3 $DIR/battlecode-manager/simple_cli.py "$@"
=======
export NODOCKER=1
python3 $DIR/battlecode-manager/simple_cli.py "$@"
>>>>>>> 6b0b8df56ef8ebdba88911a44d7374befbda1e30
