// Updated: Dec 3, 2013

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _TOOL_H_
#define _TOOL_H_

#define PI 3.141592653
#define z_red 0.5
#define z_src 2
#define annulus 0.001 / (1 + z_red)
#define delta_x 0.001 / (1 + z_red)   // Resolution of the grids in angular diameter distance unit
#define delta_y 0.001 / (1 + z_red)
#define delta_z 0.001 / (1 + z_red)
#define crit_den 1
#define rad_x 4 * 150           // 150 * delta_x Mpc/h^2 is the radius of the halo, 4 is four times the Rvir
#define rad_y 4 * 150
#define rad_z 4 * 150
#define rad (int)(4 * 150)           // 1.5 is to include those pixels outside R_vir (For multipole expansion)
#define data_link "/Users/rickyccy/Documents/Research_weak_lensing/Data/"
#define data2_link "/Users/rickyccy/Documents/Research_weak_lensing/Data/final_check_Feb_25_2014/"

int max(float arr[], int size, float subtract);
int min(float arr[], int size, float subtract);
float max_index(float arr[], int size);
float min_index(float arr[], int size);
float trap_int(float a, float b, int interval_num, float f[]);
int count_line(int spot, FILE *countfile, int line);
void skip_line(int spot, FILE *readfile, int line, int skipped_line);

#endif
