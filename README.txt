spinthon
=====

spinthon is a lightweight python module for elementary spin dynamics calculations.

At this point, arbitrary Hamiltonian's may be set up in spinthon, but support for evolution and/or the density matrix formalism is not included.


Installation
-------------
spinthon is available at pip. It may be installed
with

::

	pip install spinthon



In order to obtain spin operators in the Zeeman product basis run the following code
..code:: python

	import numpy as np
	import spinthon

	from spinthon.spin.operators import SpinOperator, VectorSpinOperator
	from spinthon.spin.system import SpinSystem

	S1 = SpinSystem(["E", "14N"])

	# S spin operators
	Sx = SpinOperator(S1, 0, "Ix")
	SxM = Sx.getMatrixRepresentation()

	# I spin operators
	Ix = SpinOperator(S1, 1, "Ix")
	IxM = Ix.getMatrixRepresentation()

	Svec = VectorSpinOperator(S1, 0)
	Ivec = VectorSpinOperator(S1, 1)

SpinOperators and VectorSpinOperators are custom types in Spinthon, for which multiplication and addition are defined as appropriate.

..code:: python
	H = Svec@np.diag([20e6, 20e6, 100e6])@Ivec

In this way an arbitrary Hamiltonian may be set up. Once set up, the eigenvectors and eigenvalues may be obtained
..code:: python
	eigenVals, eigenVecs = H.getEigenValuesAndVectors()


The eigenvectors may be used to set up a basis in which the Hamiltonian is diagonal:
..code:: python
	S1.rebase(eigenVecs)


Following this transformation the matrix representation of H will be diagonal.



Further Notes
-------------
To install this package locally and install for development, use this command

>>> pip install -e ./

More info available here: https://packaging.python.org/tutorials/installing-packages/
