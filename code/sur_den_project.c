// Parallel programming using OpenMP.

// This program projects the halos in the xy, yz and xz plane. All distance in Mpc/h
// Two particles with furthest distance from the halo center are found, then a circle is defined to have diameter being the distance between the two furthest particles. All particles at the edges of the box are excluded.
// Input: mass bin (argv[1]) half the number of pixel on length (argv[2]) (11.7, 150 pixels)
// Updated: Feb 25, 2014

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"
#include <time.h>


typedef struct {
    long int halo_id;
    float x0, y0, z0, Rvir, Reff, Mvir, partmass;
    int partnum;
} halo;

typedef struct {
    float x, y, z;
} part;


int main(int argc, char *argv[])
{
    FILE *halo_all, *halo_individual, *halo_smd_xy, *halo_smd_xz, *halo_smd_yz, *Rcut_file;
    int i, j, k, m, spot_range, individual_range, halo_count = -1, halo_number = 0, num_count = 0, halo_index = 0, pixel_Rvir_factor = 0, line_number = 0, skip_line = 0, particle_num = 0, pixel_R_size = 0, pixel_Rvir = 0;
    float Rcut = 0, **xy_grid, **xz_grid, **yz_grid, *x_coord, *y_coord, *z_coord;
    int min(float arr[], int size, float coord);
    char str[150], rootdir[150], R_str[150], rootdir2[150];
    
    
    // Input the halo catalog.
    strcpy (rootdir, data_link);
    strcpy (str, rootdir);
    strcat (str, "halo_catalog_0.2_");
    strcat (str, argv[1]);
    strcat (str, ".dat");
    
    halo_all = fopen(str, "r");
    
    // Making sure the file is present.
    if (!halo_all) {
        printf("No halo files found!\n");
        exit(1);
    }
    
    // Make sure the exact amount of arguments are passed.
    if (argc != 3) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Number of halos in that mass bin.
    halo_number = count_line(spot_range, halo_all, halo_count);
    
    halo *halos_comov = malloc(halo_number * sizeof *halos_comov);  // Comoving distance
    halo *halos_ang = malloc(halo_number * sizeof *halos_ang);      // Angular diameter distance
    

    // Get the names of the files, coordinates of the halo center, mass of each particle and virial radius of the halo.
    while(fscanf(halo_all, "%ld %f %f %f %d %f %f", &halos_comov[halo_index].halo_id, &halos_comov[halo_index].x0, &halos_comov[halo_index].y0, &halos_comov[halo_index].z0, &halos_comov[halo_index].partnum, &halos_comov[halo_index].Mvir, &halos_comov[halo_index].Rvir) == 7) {
        
        // Convert the coordinates and the radius to angular diameter distance.
        halos_ang[halo_index].halo_id = halos_comov[halo_index].halo_id;
        halos_ang[halo_index].x0 = halos_comov[halo_index].x0 / (1. + z_red);
        halos_ang[halo_index].y0 = halos_comov[halo_index].y0 / (1. + z_red);
        halos_ang[halo_index].z0 = halos_comov[halo_index].z0 / (1. + z_red);
        halos_ang[halo_index].partmass = halos_comov[halo_index].Mvir / halos_comov[halo_index].partnum;
        halos_ang[halo_index].Reff = halos_comov[halo_index].Rvir / (1. + z_red);
        
        halo_index++;
    }
    
    pixel_Rvir = (int)atof(argv[2]);
    
    halo_index = 4;
    
    // Create an array to store up the Rcutoff for the halos.
    float *Rcutoff = malloc((halo_index) * sizeof(float));
    
    double tstart, tstop, ttime;
    
    tstart = (double)clock()/CLOCKS_PER_SEC;
    
    #pragma omp parallel for default(shared) private(i,j,k, m,halo_individual,individual_range,line_number,skip_line,particle_num,Rcut,pixel_Rvir_factor, pixel_R_size, halo_smd_xy, halo_smd_yz,xy_grid,xz_grid,yz_grid,x_coord,y_coord,z_coord) schedule(dynamic)
    // For all of the halos.
    for (i = 0; i < halo_index; i++) { // for (i = 0; i < halo_index; i++)
        // Read the particles info in individual files.
        
        char str[150], skip[400], label[12];
        line_number = 0, skip_line = 0, particle_num = 0;
        strcpy (rootdir2, data2_link);
        strcpy (str, rootdir2);
        strcat (str, "Isolated_Halo_Particle_Data/");
        sprintf(label, "%ld", halos_comov[i].halo_id);
        strcat (str, label);
        strcat (str, ".csv");
        
        halo_individual = fopen(str, "r");
        
        if (!halo_individual) {
            printf("No halo files found!\n");
            exit(1);
        }
        
        // Count the number of particles.
        do {
            individual_range = fgetc(halo_individual);
            if(individual_range == '\n')
                line_number++;
        } while (individual_range != EOF);
        
        // last line doesn't end with a new line!
        // but there has to be a line at least before the last line
        if(individual_range != '\n' && line_number != 0)
            line_number++;
        
        rewind(halo_individual);
        
        num_count = line_number - 12;
        
        // Skipping the first 10 lines.  Dynamic arrays to store particle information.
        part *parts_comov = malloc(num_count * sizeof *parts_comov);
        part *parts_ang = malloc(num_count * sizeof *parts_ang);
        float *part_ang_x = malloc(num_count * sizeof(float));
        float *part_ang_y = malloc(num_count * sizeof(float));
        float *part_ang_z = malloc(num_count * sizeof(float));
        
        do {
            individual_range = fgetc(halo_individual);
            if(individual_range == '\n')
                skip_line++;
            if (skip_line == 10)
                break;
        } while (individual_range != EOF);
        
        while(fscanf(halo_individual, "%f,%f,%f", &parts_comov[particle_num].x, &parts_comov[particle_num].y, &parts_comov[particle_num].z) == 3) {
            parts_ang[particle_num].x = parts_comov[particle_num].x / (1. + z_red);
            parts_ang[particle_num].y = parts_comov[particle_num].y / (1. + z_red);
            parts_ang[particle_num].z = parts_comov[particle_num].z / (1. + z_red);
            part_ang_x[particle_num] = parts_ang[particle_num].x;
            part_ang_y[particle_num] = parts_ang[particle_num].y;
            part_ang_z[particle_num] = parts_ang[particle_num].z;
            particle_num++;
            if (particle_num == num_count)
                break;
        }

        // Cut-off radius for the halo (in angular diameter distance), and the factor of pixel_Rvir in which the halo is covered.
        Rcut = 0;
        pixel_Rvir_factor = 0;
        
        Rcut = ((max_index(part_ang_x, particle_num) - min_index(part_ang_x, particle_num)) + (max_index(part_ang_y, particle_num) - min_index(part_ang_y, particle_num)) + (max_index(part_ang_z, particle_num) - min_index(part_ang_z, particle_num))) / (2 * 3);
        pixel_Rvir_factor = (int)ceil(Rcut / (delta_x * pixel_Rvir));
        pixel_R_size = pixel_Rvir_factor * pixel_Rvir;
        
        Rcutoff[i] = Rcut;
        

        // Freeing the dynamic arrays.
        free(part_ang_x);
        free(part_ang_y);
        free(part_ang_z);
        part_ang_x = NULL;
        part_ang_y = NULL;
        part_ang_z = NULL;
        
        
        
        // Open a file to write the surface mass density.
        char smd_str_xy[150], smd_str_yz[150], smd_str_xz[150];
        
        // File for xy plane.
        strcpy (smd_str_xy, rootdir);
        strcat (smd_str_xy, "sur_den/");
        strcat (smd_str_xy, argv[1]);
        strcat (smd_str_xy, "/xy/");
        strcat (smd_str_xy, label);
        strcat (smd_str_xy, ".dat");
        halo_smd_xy = fopen(smd_str_xy, "w");
        
        // File for xz plane.
        strcpy (smd_str_xz,rootdir);
        strcat (smd_str_xz, "sur_den/");
        strcat (smd_str_xz, argv[1]);
        strcat (smd_str_xz, "/xz/");
        strcat (smd_str_xz, label);
        strcat (smd_str_xz, ".dat");
        halo_smd_xz = fopen(smd_str_xz, "w");
        
        // File for yz plane.
        strcpy (smd_str_yz, rootdir);
        strcat (smd_str_yz, "sur_den/");
        strcat (smd_str_yz, argv[1]);
        strcat (smd_str_yz, "/yz/");
        strcat (smd_str_yz, label);
        strcat (smd_str_yz, ".dat");
        halo_smd_yz = fopen(smd_str_yz, "w");
        

        // Grids
        xy_grid = (float **)malloc(2 * pixel_R_size * sizeof(float*));
        for (k = 0; k < 2 * pixel_R_size; k++) {
            xy_grid[k] = (float*)malloc(2 * pixel_R_size * sizeof(float));
        }
        xz_grid = (float **)malloc(2 * pixel_R_size * sizeof(float*));
        for (k = 0; k < 2 * pixel_R_size; k++) {
            xz_grid[k] = (float*)malloc(2 * pixel_R_size * sizeof(float));
        }
        yz_grid = (float **)malloc(2 * pixel_R_size * sizeof(float*));
        for (k = 0; k < 2 * pixel_R_size; k++) {
            yz_grid[k] = (float*)malloc(2 * pixel_R_size * sizeof(float));
        }
        // Coordinates of the grids.
        x_coord = (float *)malloc(2 * pixel_R_size * sizeof(float));
        y_coord = (float *)malloc(2 * pixel_R_size * sizeof(float));
        z_coord = (float *)malloc(2 * pixel_R_size * sizeof(float));
        
        
        
        // Put in the coordinates of the grids.
        for (j = 0; j < 2 * pixel_R_size; j++) {
            x_coord[j] = (halos_ang[i].x0 + (delta_x / 2) + (j - pixel_R_size) * delta_x);
        }
        for (j = 0; j < 2 * pixel_R_size; j++) {
            y_coord[j] = (halos_ang[i].y0 + (delta_y / 2) + (j - pixel_R_size) * delta_y);
        }
        for (j = 0; j < 2 * pixel_R_size; j++) {
            z_coord[j] = (halos_ang[i].z0 + (delta_z / 2) + (j - pixel_R_size) * delta_z);
        }
        
        
        
        // The grids
        for (m = 0; m < 2 * pixel_R_size; m++) {
            for (k = 0; k < 2 * pixel_R_size; k++) {
                xy_grid[m][k] = 0;
            }
        }
        
        for (m = 0; m < 2 * pixel_R_size; m++) {
            for (k = 0; k < 2 * pixel_R_size; k++) {
                xz_grid[m][k] = 0;
            }
        }
        
        for (m = 0; m < 2 * pixel_R_size; m++) {
            for (k = 0; k < 2 * pixel_R_size; k++) {
                yz_grid[m][k] = 0;
            }
        }
       
        // Project the surface mass density.
        for (j = 0; j < particle_num; j++) {
            if (pow((pow((parts_ang[j].x - halos_ang[i].x0), 2) + pow((parts_ang[j].y - halos_ang[i].y0), 2)), 0.5) <= Rcut) {
                xy_grid[min(y_coord, 2 * pixel_R_size, parts_ang[j].y)][min(x_coord, 2 * pixel_R_size, parts_ang[j].x)] += halos_ang[i].partmass / (delta_x * delta_y);
            }
        }
        
        
        for (j = 0; j < particle_num; j++) {
            // Project in xz plane.
            if (pow((pow((parts_ang[j].x - halos_ang[i].x0), 2) + pow((parts_ang[j].z - halos_ang[i].z0), 2)), 0.5) <= Rcut) {
                xz_grid[min(z_coord, 2 * pixel_R_size, parts_ang[j].z)][min(x_coord, 2 * pixel_R_size, parts_ang[j].x)] += halos_ang[i].partmass / (delta_x * delta_z);
            }
        }
        
        
        for (j = 0; j < particle_num; j++) {
            // Project in yz plane.
            if (pow((pow((parts_ang[j].y - halos_ang[i].y0), 2) + pow((parts_ang[j].z - halos_ang[i].z0), 2)), 0.5) <= Rcut) {
                yz_grid[min(z_coord, 2 * pixel_R_size, parts_ang[j].z)][min(y_coord, 2 * pixel_R_size, parts_ang[j].y)] += halos_ang[i].partmass / (delta_z * delta_y);
            }
        }
        
        for (m = 0; m < 2 * pixel_R_size; m++) {
            for (k = 0; k < 2 * pixel_R_size; k++) {
                fprintf(halo_smd_xy, "%f\t", xy_grid[m][k]);
            }
            fprintf(halo_smd_xy, "\n");
        }
        for (m = 0; m < 2 * pixel_R_size; m++) {
            for (k = 0; k < 2 * pixel_R_size; k++) {
                fprintf(halo_smd_xz, "%f\t", xz_grid[m][k]);
            }
            fprintf(halo_smd_xz, "\n");
        }
        for (m = 0; m < 2 * pixel_R_size; m++) {
            for (k = 0; k < 2 * pixel_R_size; k++) {
                fprintf(halo_smd_yz, "%f\t", yz_grid[m][k]);
            }
            fprintf(halo_smd_yz, "\n");
        }
        

        fclose(halo_smd_xy);
        fclose(halo_smd_xz);
        fclose(halo_smd_yz);
        fclose(halo_individual);
        
        free(parts_comov);
        free(parts_ang);
        free(xy_grid);
        free(xz_grid);
        free(yz_grid);
        free(x_coord);
        free(y_coord);
        free(z_coord);
        parts_comov = NULL;
        parts_ang = NULL;
        xy_grid = NULL;
        xz_grid = NULL;
        yz_grid = NULL;
        x_coord = NULL;
        y_coord = NULL;
        z_coord = NULL;
        
    }
    
    tstop = (double)clock()/CLOCKS_PER_SEC;
    
    ttime = tstop - tstart;
    
    printf("%f\n", ttime);
    
    strcpy (R_str, rootdir);
    strcat (R_str, argv[1]);
    strcat (R_str, "_cutoff.dat");
    
    Rcut_file = fopen(R_str, "w");
    
    if (!Rcut_file) {
        printf("No cut off file can be opened!\n");
        exit(1);
    }
    
    for (i = 0; i < halo_index; i++) {
        fprintf(Rcut_file, "%ld\t%f\n", halos_ang[i].halo_id, Rcutoff[i]);
    }

    
    
    free(halos_ang);
    free(halos_comov);
    free(Rcutoff);
    halos_ang = NULL;
    halos_comov = NULL;
    Rcutoff = NULL;
    
    
    fclose(halo_all);
    fclose(Rcut_file);
     
    
    return 0;
}
