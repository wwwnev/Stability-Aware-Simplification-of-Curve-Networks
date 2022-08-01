#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
// Defining constants here (instead of changing parameters in multiple different files)

const double radius = 2; 
const double young_mod = 1;

const double gc_stretching_stiffness = young_mod * M_PI * radius * radius; // N. 
const double gc_bending_stiffness = (gc_stretching_stiffness / 4) * radius * radius; // N * mm^2. 
//const double gc_stretching_stiffness = 500;
//const double gc_bending_stiffness = 1;


const bool gc_fix_endpoints = true;
const bool gc_do_bend = true;
const bool gc_do_stretch = true;

const double gc_endpoint_penalty_mu = 500; // N/mm^2.

const double gc_resample_length = 2.0;