You need python2 to run the code (python3 will not work). It should be installed by default in UNIX systems.

The needed libraries can be installed with:
pip install -r requirements.txt

Or the standard Canopy (https://www.enthought.com/products/canopy/) should have them installed.

The code is organized in two main parts:
- runSimulations.py
- studyResults.py

They require running costBenefit.py once before to create the files containing the growth rates for all protein concentrations and stress levels (already created).

== runSimulations.py
Just uncomment in main() the data you want to create.

== studyResults.py
Just uncomment in main() the figure you want to create.
