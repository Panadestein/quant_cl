# TITLE: A quantum interpreter written in Lisp

This repository contains `clq`, my implementation in Common Lisp of a quantum computer
interpreter. It is based on the excellent tutorial by Robert Smith, which you can find [here](https://www.stylewarning.com/posts/quantum-interpreter/).

The interpreter's correctness has been tested by comparing with the results of the
IBM's Qiskit Aer simulator, and also the results in a real quantum computer
(IBM's Eagle QPU).

I have made a literate programming implementation, so the code is embedded in the docs.
I guess if you are here I don't have to say it, but it is best to explore the
interpreter with Emacs, Slime and Org-mode.

## Overview

The interpreter's idea is simple: evolve the wave function $\Psi_0^n$ in the full
Hilbert space, that means $2^n \times 2^n$ if we use $n$ qubits. Then construct
the full gates by lifting:

$L_U = I \otimes \cdots U \cdots \otimes I$

and evolve it:

$$\Psi_m^n = \Psi_0^n \prod_i^m L_{U_i}$$

The measurement is done by sampling the CDF of the squared amplitudes of the states.
But please see the org-mode file for the gory details.
