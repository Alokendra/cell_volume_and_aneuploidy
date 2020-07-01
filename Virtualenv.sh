#!/bin/bash

# Script to setup a virtual environment for viewing jupyter notebook
# virtual environment must be installed, python interpreter should
# be > 2.7.15

set -e

VIRTUALENVNAME=pyenv
NOTEBOOKFILE="Cell_Volume_Aneuploidy_Osmotic_Stress.md"

mkdir -p docs figures
virtualenv ./$VIRTUALENVNAME
echo "Virtual environment $VIRTUALENVNAME installed, switching"
source ./$VIRTUALENVNAME/bin/activate
echo "New python executable:"
which python
read -p "Install necessary packages?(y/n): " Res
if [ "$Res" != "y" ] ; then
    echo "Exiting...."
    exit 1
fi
echo "Installing packages, this may take a while...."
pip install -r requirements.txt > Piplog
Status=$?
if [ $Status -eq 0 ] ; then
    echo "Installation successful"
    echo "Syncing jupyter notebook..."
    jupytext --sync "$NOTEBOOKFILE"
    NOTEBOOKIPYNB="${NOTEBOOKFILE%.md}".ipynb
    read -p "Open notebook file $NOTEBOOKIPYNB ?(y/n): " Op
    if [ "$Op" = "y" ] ; then
	jupyter notebook "$NOTEBOOKIPYNB"
    else
	Command="source pyenv/bin/activate ; jupyter notebook $NOTEBOOKIPYNB"
	echo "Open the notebook with command: \"$Command\""
	echo "End of script.."
    fi
else
    echo "Installation unsuccessful, see the log for details"
fi


