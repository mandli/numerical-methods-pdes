#!/usr/bin/env
"""Basic test script that runs all notebooks."""

import os
import subprocess
import tempfile
import sys
import nbformat
import doctest
import glob

if sys.version_info >= (3,0):
    kernel = 'python3'
else:
    kernel = 'python2'

# Note we leave out the python intro as there are purposeful exceptions and
# the finite volume lecture as there are packages that are not installed.
notebooks = glob.glob("*.ipynb")
notebooks.remove('01_python.ipynb')
notebooks.remove('13_finite_volume.ipynb')

def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=60",
                "--ExecutePreprocessor.kernel_name="+kernel,
                "--output", fout.name, path]
        subprocess.check_call(args)

        fout.seek(0)
        nb = nbformat.reads(fout.read().decode('utf-8'), nbformat.current_nbformat)

    errors = [output for cell in nb.cells if "outputs" in cell
              for output in cell["outputs"]
              if output.output_type == "error"]

    return nb, errors

def run_tests():

    for filename in notebooks:
        if (filename.split('.')[-1] == 'ipynb'):
            nb, errors = _notebook_run(filename)
            if errors != []:
                raise(Exception)


if __name__ == '__main__':
    run_tests()