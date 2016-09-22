// This program projects the halos in the xy, yz and xz plane, calculates the multipole terms at all radius bins. All distance in Mpc/h
// Input: mass bin (argv[1]) half the number of pixel on length (argv[2]) (11.7, 150 pixels)
// Updated: Mar 6, 2014

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"


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
    FILE *halo_all, *halo_individual, *halo_smd_xy, *halo_smd_xz, *halo_smd_yz, *halo_random, *Rcut_file;
    int i, j, k, m, spot_range, individual_range, cut_range, halo_number = -1, halo_index = 0, mini_xy, mini_yz, mini_xz, pixel_Rvir = 0, pixel_Rvir_factor = 0, pixel_R_size = 0, mult_R_bin = 0, cut_count = -1, cut_number = 0, cut_index = 0, line_number = 0, skip_line = 0, particle_num = 0, num_count = 0;
    float Rcut = 0, e = 0, rad_particle_xy, rad_particle_xz, rad_particle_yz, x_diff, y_diff, z_diff, xy_theta, yz_theta, xz_theta, *xy_sig_2_re, *xz_sig_2_re, *yz_sig_2_re, *xy_coord, *xz_coord, *yz_coord;
    double *xy_sig_0, *xz_sig_0, *yz_sig_0;
    int min(float arr[], int size, float coord);
    char str[150], rootdir[150], R_str[150], rootdir2[150];
    
    // Input the halo catalog.
    strcpy (rootdir, data2_link);
    strcpy (str, rootdir);
    strcat (str, "halo_catalog_0.2_");
    strcat (str, argv[1]);
    strcat (str, ".dat");
    
    halo_all = fopen(str, "r");
    
    if (!halo_all) {
        printf("No halo catalog found!\n");
        exit(1);
    }
    
    // Make sure the exact amount of arguments are passed.
    if (argc != 3) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Count the number of files.
    do {
        spot_range = fgetc(halo_all);
        if(spot_range == '\n')
            halo_number++;
    } while (spot_range != EOF);
    
    // last line doesn't end with a new line!
    // but there has to be a line at least before the last line
    if(spot_range != '\n' && halo_number != 0)
        halo_number++;
    
    rewind(halo_all);
    
    
    halo *halos_comov = malloc(halo_number * sizeof *halos_comov);  // Comoving distance
    halo *halos_ang = malloc(halo_number * sizeof *halos_ang);      // Angular diameter distance
    
    // Get the names of the files.
    while(fscanf(halo_all, "%ld %f %f %f %d %f %f", &halos_comov[halo_index].halo_id, &halos_comov[halo_index].x0, &halos_comov[halo_index].y0, &halos_comov[halo_index].z0, &halos_comov[halo_index].partnum, &halos_comov[halo_index].Mvir, &halos_comov[halo_index].Rvir) == 7) {
        
        // Convert the coordinates and the radius to angular diameter distance.
        halos_ang[halo_index].halo_id = halos_comov[halo_index].halo_id;
        halos_ang[halo_index].x0 = halos_comov[halo_index].x0 / (1. + z_red);
        halos_ang[halo_index].y0 = halos_comov[halo_index].y0 / (1. + z_red);
        halos_ang[halo_index].z0 = halos_comov[halo_index].z0 / (1. + z_red);
        halos_ang[halo_index].partmass = 1;
        halos_ang[halo_index].Reff = halos_comov[halo_index].Rvir / (1. + z_red);
        
        halo_index++;
    }
    
    pixel_Rvir = (int)atof(argv[2]);

    
    // Read in the R cutoff.
    strcpy (R_str, rootdir);
    strcat (R_str, argv[1]);
    strcat (R_str, "_cutoff.dat");
    
    Rcut_file = fopen(R_str, "r");
    
    if (!Rcut_file) {
        printf("No cut off file found!\n");
        exit(1);
    }

    // Count the number of files in the cut off file, and make sure it agrees with the number in the halo catalog.
    cut_number = count_line(cut_range, Rcut_file, cut_count);
    
    // Create an array to store up the file name and the cut off radius.
    long int *cut_halo = malloc(cut_number * sizeof(long int));
    float *Rcutoff = malloc(cut_number * sizeof(float));
    
    while(fscanf(Rcut_file, "%ld %f", &cut_halo[cut_index], &Rcutoff[cut_index]) == 2) {
        cut_index++;
    }

    fclose(Rcut_file);
    
    // For all of the halos.
    for (i = 0; i < halo_index; i++) { // for (i = 0; i < halo_index; i++)
        
        // Read the particles info in individual files.
        char str[150], skip[400], label[12];
        line_number = 0, skip_line = 0, particle_num = 0;
        strcpy (rootdir2, data_link);
        strcpy (str, rootdir2);
        strcat (str, "Isolated_Halo_Particle_Data/");
        sprintf(label, "%ld", halos_ang[i].halo_id);
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
        
        num_count = line_number - 12
        
        // Skipping the first 10 lines.  Defining the dynamic arrays.
        part *parts_comov = malloc(num_count * sizeof *parts_comov);
        part *parts_ang = malloc(num_count * sizeof *parts_ang);
        float *part_ang_x2 = malloc(num_count * sizeof(float));
        float *part_ang_y2 = malloc(num_count * sizeof(float));
        float *part_ang_z2 = malloc(num_count * sizeof(float));
        
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
            part_ang_x2[particle_num] = parts_ang[particle_num].x;
            part_ang_y2[particle_num] = parts_ang[particle_num].y;
            part_ang_z2[particle_num] = parts_ang[particle_num].z;
            particle_num++;

            if (particle_num == num_count)   // To avoid the last line #OK
                break;
            
        }        
        
        // Cut-off radius for the halo (in angular diameter distance), and the factor of pixel_Rvir in which the halo is covered.
        Rcut = Rcutoff[i];
        pixel_Rvir_factor = 0;

        // Number of pixels covering radius of the halo (in terms of R_vir)
        pixel_Rvir_factor = (int)ceil(Rcut / (delta_x * pixel_Rvir));
        // Number of pixels covering radius of the halo
        pixel_R_size = pixel_Rvir_factor * pixel_Rvir;
        
        // The number of radius bin for multipole expansion. (Only the shears in pixel_Rvir is to be trusted.)
        mult_R_bin = 8 * pixel_Rvir;
        
        free(part_ang_x2);
        free(part_ang_y2);
        free(part_ang_z2);
        part_ang_x2 = NULL;
        part_ang_y2 = NULL;
        part_ang_z2 = NULL;
        
        // Open a file to write the multipole expansion.
        char smd_str_xy[150], smd_str_yz[150], smd_str_xz[150];
        
        // File for xy plane.
        strcpy (smd_str_xy, rootdir);
        strcat (smd_str_xy, "multi_exp/");
        strcat (smd_str_xy, argv[1]);
        strcat (smd_str_xy, "/xy/");
        strcat (smd_str_xy, label);
        strcat (smd_str_xy, ".dat");
        halo_smd_xy = fopen(smd_str_xy, "w");
        
        // File for xz plane.
        strcpy (smd_str_xz, rootdir);
        strcat (smd_str_xz, "multi_exp/");
        strcat (smd_str_xz, argv[1]);
        strcat (smd_str_xz, "/xz/");
        strcat (smd_str_xz, label);
        strcat (smd_str_xz, ".dat");
        halo_smd_xz = fopen(smd_str_xz, "w");
        
        // File for yz plane.
        strcpy (smd_str_yz, rootdir);
        strcat (smd_str_yz, "multi_exp/");
        strcat (smd_str_yz, argv[1]);
        strcat (smd_str_yz, "/yz/");
        strcat (smd_str_yz, label);
        strcat (smd_str_yz, ".dat");
        halo_smd_yz = fopen(smd_str_yz, "w");
        
        // Creating arrays for the multipole terms
        xy_sig_0 = (double *)malloc(mult_R_bin * sizeof(double));
        yz_sig_0 = (double *)malloc(mult_R_bin * sizeof(double));
        xz_sig_0 = (double *)malloc(mult_R_bin * sizeof(double));
        xy_sig_2_re = (float *)malloc(mult_R_bin * sizeof(float));
        yz_sig_2_re = (float *)malloc(mult_R_bin * sizeof(float));
        xz_sig_2_re = (float *)malloc(mult_R_bin * sizeof(float));
        
        // Coordinates of the grids.
        xy_coord = (float *)malloc(mult_R_bin * sizeof(float));
        xz_coord = (float *)malloc(mult_R_bin * sizeof(float));
        yz_coord = (float *)malloc(mult_R_bin * sizeof(float));
        
        // Put in the radius of the annulus.
        for (j = 0; j < mult_R_bin; j++) {
            xy_coord[j] = (j + 0.5) * annulus;
            xz_coord[j] = (j + 0.5) * annulus;
            yz_coord[j] = (j + 0.5) * annulus;
        }

        // Initializing sigma0 and sigma2
        for (m = 0; m < mult_R_bin; m++) {
            xy_sig_0[m] = 0;
            xz_sig_0[m] = 0;
            yz_sig_0[m] = 0;
            xy_sig_2_re[m] = 0;
            xz_sig_2_re[m] = 0;
            yz_sig_2_re[m] = 0;
        }
        
        // Project the sigma0 and sigma2 in xy-plane.
        for (j = 0; j < particle_num; j++) {
            
            // Project in xy plane.
            rad_particle_xy = pow((pow((parts_ang[j].x - halos_ang[i].x0), 2) + pow((parts_ang[j].y - halos_ang[i].y0), 2)), 0.5);
            
            if (rad_particle_xy <= Rcut) {
                mini_xy = min(xy_coord, mult_R_bin, rad_particle_xy);
                xy_theta = atan2((parts_ang[j].y - halos_ang[i].y0), (parts_ang[j].x - halos_ang[i].x0));   // Theta of the particle respect to the center
                if (mini_xy == 0) {
                    xy_sig_0[mini_xy] += halos_ang[i].partmass / (M_PI * (pow(annulus, 2)));
                    xy_sig_2_re[mini_xy] += 2 * ((halos_ang[i].partmass * cos(2 * xy_theta)) / (M_PI * (pow(annulus, 2))));
                }
                else {
                    xy_sig_0[mini_xy] += halos_ang[i].partmass / (M_PI * (pow(((mini_xy + 1) * annulus), 2) - pow(((mini_xy) * annulus), 2)));
                    xy_sig_2_re[mini_xy] += 2 * ((halos_ang[i].partmass * cos(2 * xy_theta)) / (M_PI * (pow(((mini_xy + 1) * annulus), 2) - pow((mini_xy * annulus), 2))));
                }
            }
        }
        
        // Project the sigma0 and sigma2 in xz-plane.
        for (j = 0; j < particle_num; j++) {
            
            // Project in xz plane.
            rad_particle_xz = pow((pow((parts_ang[j].x - halos_ang[i].x0), 2) + pow((parts_ang[j].z - halos_ang[i].z0), 2)), 0.5);
            
            if (rad_particle_xz <= Rcut) {
                mini_xz = min(xz_coord, mult_R_bin, rad_particle_xz);
                xz_theta = atan2((parts_ang[j].z - halos_ang[i].z0), (parts_ang[j].x - halos_ang[i].x0));   // Theta of the particle respect to the center
                if (mini_xz == 0) {
                    xz_sig_0[mini_xz] += halos_ang[i].partmass / (M_PI * (pow(annulus, 2)));
                    xz_sig_2_re[mini_xz] += 2 * ((halos_ang[i].partmass * cos(2 * xz_theta)) / (M_PI * (pow(annulus, 2))));
                }
                else {
                    xz_sig_0[mini_xz] += halos_ang[i].partmass / (M_PI * (pow(((mini_xz + 1) * annulus), 2) - pow((mini_xz * annulus), 2)));
                    xz_sig_2_re[mini_xz] += 2 * ((halos_ang[i].partmass * cos(2 * xz_theta)) / (M_PI * (pow(((mini_xz + 1) * annulus), 2) - pow((mini_xz * annulus), 2))));
                }
            }
        }
        
        // Project the sigma0 and sigma2 in yz-plane.
        for (j = 0; j < particle_num; j++) {
            
            // Project in xy plane.
            rad_particle_yz = pow((pow((parts_ang[j].y - halos_ang[i].y0), 2) + pow((parts_ang[j].z - halos_ang[i].z0), 2)), 0.5);
            
            if (rad_particle_yz <= Rcut) {
                mini_yz = min(yz_coord, mult_R_bin, rad_particle_yz);
                yz_theta = atan2((parts_ang[j].z - halos_ang[i].z0), (parts_ang[j].y - halos_ang[i].y0));   // Theta of the particle respect to the center
                if (mini_yz == 0) {
                    yz_sig_0[mini_yz] += halos_ang[i].partmass / (M_PI * (pow(annulus, 2)));
                    yz_sig_2_re[mini_yz] += 2 * ((halos_ang[i].partmass * cos(2 * yz_theta)) / (M_PI * (pow(annulus, 2))));
                }
                else {
                    yz_sig_0[mini_yz] += halos_ang[i].partmass / (M_PI * (pow(((mini_yz + 1) * annulus), 2) - pow((mini_yz * annulus), 2)));
                    yz_sig_2_re[mini_yz] += 2 * ((halos_ang[i].partmass * cos(2 * yz_theta)) / (M_PI * (pow(((mini_yz + 1) * annulus), 2) - pow((mini_yz * annulus), 2))));
                }
            }
        }
        
        // Absolute sign for sigma 2
        for (m = 0; m < mult_R_bin; m++) {
            fprintf(halo_smd_xy, "%f\t%E\t%E\n", xy_coord[m], xy_sig_0[m], fabs(xy_sig_2_re[m]));
            fprintf(halo_smd_xz, "%f\t%E\t%E\n", xz_coord[m], xz_sig_0[m], fabs(xz_sig_2_re[m]));
            fprintf(halo_smd_yz, "%f\t%E\t%E\n", yz_coord[m], yz_sig_0[m], fabs(yz_sig_2_re[m]));
        }
        
        fclose(halo_smd_xy);
        fclose(halo_smd_xz);
        fclose(halo_smd_yz);
        fclose(halo_individual);
        
        free(parts_comov);
        free(parts_ang);
        free(xy_sig_0);
        free(xz_sig_0);
        free(yz_sig_0);
        free(xy_sig_2_re);
        free(xz_sig_2_re);
        free(yz_sig_2_re);
        free(xy_coord);
        free(xz_coord);
        free(yz_coord);
        parts_comov = NULL;
        parts_ang = NULL;
        xy_sig_0 = NULL;
        xz_sig_0 = NULL;
        yz_sig_0 = NULL;
        xy_sig_2_re = NULL;
        xz_sig_2_re = NULL;
        yz_sig_2_re = NULL;
        xy_coord = NULL;
        xz_coord = NULL;
        yz_coord = NULL;
        
        // Check progress
        if (i % 50 == 49) {
            printf("%d out of %d files completed!\t%s\n", i + 1, halo_index, label);
        }
    }
    
    free(halos_comov);
    free(halos_ang);
    free(Rcutoff);
    free(cut_halo);
    halos_comov = NULL;
    halos_ang = NULL;
    Rcutoff = NULL;
    cut_halo = NULL;
    
    fclose(halo_all);
    
    return 0;
}

