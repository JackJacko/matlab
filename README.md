# matlab
Repository for MATLAB code

Sections A1-A4 and B1-B2
------------------------

These files contain the code for my Mie-base size and RI determination procedure (A1-A4) and the code for Mie-based IOP forward modelling.

Note that the various sections need to be put together and are further dependent on a few functions to run properly, which I'm not including here. This is done on purpose as a soft-block to third-party appropriation (i.e. I'm ok with people asking for the code, and I guess I'm also ok with people that for whatever reason do not want to ask for it and instead decide to put in the time to construct a makeshift functioning version of it - kudos to you I suppose, but why?.)

circledrop.m
------------

A piece of code I made to solve a coding problem I found on the internet once. It just computes the final configuration of an arbitrary number of circles with random x positioning and an arbitrary radius as if they fell one by one from infinite height down to y = 0.

occlusion3d.m
-------------

This contains a potentially useful code I made after talking with my supervisor about fractal particles. As it is it plots power law distributions of spheres within a cubic 3D volume and their 2D projection, highlighting the occlusion that happens as particles hide behind each other.

cnw_B3_S23.m
------------

A clunky MATLAB version of of Conway’s Game of Life.

expungeText.m, cCipher.m and vigenere.m
---------------------------------------
 
Functions dedicated to Caesar and Vigenère cipher encoding and decoding.