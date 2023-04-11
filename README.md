Min Cost Flow library
=====================

This is a proof of concept library for solving the Min Cost Flow (MCF) problem
for linear and convex costs functions.
The code and API are highly specialized for the problem domain that involves the
*Uncertainty Networks*, a concept invented by Rene Pickhardt
and implemented in
[pickhardtpayments](github.com/renepickhardt/pickhardtpayments).

The algorithms used here are taken from Ahuja, et. al book "Network Flows:
Theory, Algorithms and Applications", 1993.

Here are some benchmark test I run on my solvers, compared against Google's OR-Tools library.
These results are produced with random graph generated using `networkx` (see `benchmark/gen.py`
for details) using a variable number *N* of nodes and a *M = 7.5 N*.
40 different graphs are generated for each value of *N*.
Costs and capacities for arcs are selected as random numbers from 0 to 200.
![](/assets/benchmark.png)

The C++ API will be out soon!

The C API is described [here](c-api.md).

The Python API is described [here](python-api.md).

Author: Eduardo Quintana Miranda

This work is funded by the [Summer of Bitcoin](www.summerofbitcoin.org) project 2023 edition.
