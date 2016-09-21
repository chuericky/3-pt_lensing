// Parallelized code in OpenMP
// This program calculates the wiener_filter extracts the ellipticity of the halo by using the wiener filter method, doing stacking as well.
//  mass bin (argv[1]), plane (argv[2]), pixel across R_vir (argv[3]), half number of pix across r_min (argv[4]) (for 300 pixels, 50), halo no. to be calculated (in units of 50) (For example: 0 to 49: this input should be 0; 50 to 99: should be 1) (argv[5])
// Updated: May 18, 2014


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"

int main(int argc, char *argv[])
{
    FILE *gtan_a, *gcross_a, *gtan_F, *gcross_F, *gtan_m, *dataname, *multipole, *epis_file;
    char gtan_a_str[200], gcross_a_str[200], gtan_F_str[200], gcross_F_str[200], gtan_mon_str[200], multi_str[200], file_str[200], epis_str[200], pixel_no[5], m_bin[5], plane[3], halo_na[13], NN_str[5];
    int i = 0, j = 0, k = 0, m = 0, half_pix = 0, halo_no = 0, halo_index = 0, halo_number = 0, spot_range = 0, pix_size = 0, pix_line = 0, pix_row = 0, pix_col = 0, pix_tot = 0, multi_size = 0, multi_line = 0, interval_num = 0, min_ind = 0, NN = 0, halo_est = 0, pile = 0;
    long int pixpair_count = 0;
    float r_min = 0, r_max = 0, cos_db = 0, sin_db = 0, integral = 0, bin = 0;
    double est_sum = 0;

    // Make sure the exact number of arguments are passed.
    if (argc != 6) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Input the mass bin and the projected plane.
    strcpy (m_bin, argv[1]);
    strcpy (plane, argv[2]);
    strcpy (NN_str, argv[5]);
    NN = atoi(NN_str);   // which halos to be estimated
    
    // Get the names of the files to know which halos are to be estimated.
    strcpy (file_str, data2_link);
    strcat (file_str, m_bin);
    //strcat (file_str, "_cutoff_ran.dat");
    strcat (file_str, "_cutoff.dat");    

    dataname = fopen (file_str, "r");
    
    if (!dataname) {
        printf("No halo catalog file %s found!\n", file_str);
        exit(1);
    }
    
    halo_no = -1;
    halo_index = 0;
    
    // Count the number of files.
    halo_number = count_line(spot_range, dataname, halo_no);
    
    long int *halo_name = malloc(halo_number * sizeof(long int));
        
    // Get the names of the halos.
    while(fscanf(dataname, "%ld %*f", &halo_name[halo_index]) == 1) {
        halo_index++;
    }
    
    fclose(dataname);
    
    // The pixel across R_vir (2 for diameter, 8 for FFT region)
    half_pix = 2 * 8 * atoi(argv[3]);
    sprintf(pixel_no, "%d", half_pix);
    
    // Minimum index of pixel that is ignored.
    min_ind = atoi(argv[4]);
    
    // An array to store up the epis of the halos. For stacking, the second array is needed.
    float *epis = malloc(halo_index * sizeof(float));
    float *stack_est = malloc(halo_index * sizeof(float));
    float *inter = malloc(halo_index * sizeof(float));
    float *rmin = malloc(halo_index * sizeof(float));
    float *rmax = malloc(halo_index * sizeof(float));
    
    for (k = 0; k < halo_index; k++) {
        epis[k] = 0;
        stack_est[k] = 0;
        inter[k] = 0;
        rmin[k] = 0;
        rmax[k] = 0;
    }
    
    // Open the monopole tangential shear file and subtract it away
    strcpy (gtan_mon_str, data_link);
    strcat (gtan_mon_str, "shear/");
    strcat (gtan_mon_str, m_bin);
    strcat (gtan_mon_str, "/");
    strcat (gtan_mon_str, plane);
    strcat (gtan_mon_str, "/stack_g_tan_mon.dat");

    gtan_m = fopen (gtan_mon_str, "r");
    
    if (!gtan_m) {
        printf("No monopole tan shear file %s found!\n", gtan_mon_str);
        exit(1);
    }
    
    pix_size = -1;
    pix_line = 0;
    
    // Count the number of rows of pixels, columns as well.
    pix_row = count_line(spot_range, gtan_m, pix_size);
    pix_col = pix_row;
    pix_tot = pix_row * pix_col;
    
    float **gtan_m_mat = malloc(pix_row * sizeof(float *));
    
    for (i = 0; i < pix_row; i++) {
        gtan_m_mat[i] = malloc(pix_col * sizeof(float));
        
        for (m = 0; m < pix_col; m++) {
            fscanf(gtan_m, "%f", &gtan_m_mat[i][m]);
        }
    }
    
    fclose(gtan_m);
    
    // Calculation is easier to calculate the estimator in 1D array.
    float *gtan_m_arr = malloc(pix_tot * sizeof(float));
    
    j = 0;
    
    for (i = 0; i < pix_row; i++) {
        for (m = 0; m < pix_col; m++) {
            gtan_m_arr[j] = gtan_m_mat[i][m];
            j++;
        }
    }
    
    free(gtan_m_mat);
    gtan_m_mat = NULL;

    pile = 50;
    
    if ((NN + 1) * pile > halo_index - 1) {
	halo_est = halo_index - 1;
    }
    else {
	halo_est = (NN + 1) * pile;
    }
    
    strcpy (gtan_a_str, data_link);
    strcat (gtan_a_str, "shear/");
    strcat (gtan_a_str, m_bin);
    strcat (gtan_a_str, "/");
    strcat (gtan_a_str, plane);
    strcat (gtan_a_str, "/t_a/12345678910.dat");
     
    strcpy (gcross_a_str, data_link);
    strcat (gcross_a_str, "shear/");
    strcat (gcross_a_str, m_bin);
    strcat (gcross_a_str, "/");
    strcat (gcross_a_str, plane);
    strcat (gcross_a_str, "/c_a/12345678910.dat");
     
    gtan_a = fopen (gtan_a_str, "r");
    gcross_a = fopen (gcross_a_str, "r");
    
    if (!gtan_a) {
        printf("No analytical tan shear file %s found!\n", gtan_a_str);
        exit(1);
    }
    if (!gcross_a) {
        printf("No analytical cross shear file %s found!\n", gcross_a_str);
        exit(1);
    }
    pix_size = -1;
    pix_line = 0;
    
    // Count the number of rows of pixels, columns as well.
    pix_row = count_line(spot_range, gtan_a, pix_size);
    pix_col = pix_row;
    pix_tot = pix_row * pix_col;
    
    // Input the gamma matrix profile
    float **gtan_a_mat = malloc(pix_row * sizeof(float *));
    float **gcross_a_mat = malloc(pix_row * sizeof(float *));
    
    for (i = 0; i < pix_row; i++) {
        gtan_a_mat[i] = malloc(pix_col * sizeof(float));
        gcross_a_mat[i] = malloc(pix_col * sizeof(float));
        
        for (m = 0; m < pix_col; m++) {
            fscanf(gtan_a, "%f", &gtan_a_mat[i][m]);
        }
        for (m = 0; m < pix_col; m++) {
            fscanf(gcross_a, "%f", &gcross_a_mat[i][m]);
        }
    }

    fclose(gtan_a);
    fclose(gcross_a);
   
    // Calculation is easier to calculate the estimator in 1D array.
    float *gtan_a_arr = malloc(pix_tot * sizeof(float));
    float *gcross_a_arr = malloc(pix_tot * sizeof(float));
    j = 0;

    for (i = 0; i < pix_row; i++) {
        for (m = 0; m < pix_col; m++) {
            gtan_a_arr[j] = gtan_a_mat[i][m];
            gcross_a_arr[j] = gcross_a_mat[i][m];
            j++;
        }
    }
    
    free(gtan_a_mat);
    free(gcross_a_mat);
    gtan_a_mat = NULL;
    gcross_a_mat = NULL;
    
    for (k = NN * pile; k < halo_est; k++) {   //k = NN * 50; k < halo_est
    
        // Open the FFT and analytical shear files (4 files are to be opened)
        sprintf (halo_na, "%ld", halo_name[k]);
    
        strcpy (gtan_F_str, data_link);
        strcat (gtan_F_str, "shear/");
        strcat (gtan_F_str, m_bin);
        strcat (gtan_F_str, "/");
        strcat (gtan_F_str, plane);
        strcat (gtan_F_str, "/t_F/");
        strcat (gtan_F_str, halo_na);
        strcat (gtan_F_str, ".dat");
    
        strcpy (gcross_F_str, data_link);
        strcat (gcross_F_str, "shear/");
        strcat (gcross_F_str, m_bin);
        strcat (gcross_F_str, "/");
        strcat (gcross_F_str, plane);
        strcat (gcross_F_str, "/c_F/");
        strcat (gcross_F_str, halo_na);
        strcat (gcross_F_str, ".dat");
        
        gtan_F = fopen (gtan_F_str, "r");
        gcross_F = fopen (gcross_F_str, "r");
        
        if (!gtan_F) {
            printf("No FFT tan shear file %s found!\n", gtan_F_str);
            exit(1);
        }
        if (!gcross_F) {
            printf("No FFT cross shear file %s found!\n", gcross_F_str);
            exit(1);
        }

        // Input the gamma matrix profile
        float **gtan_F_mat = malloc(pix_row * sizeof(float *));
        float **gcross_F_mat = malloc(pix_row * sizeof(float *));
        
        for (i = 0; i < pix_row; i++) {
            gtan_F_mat[i] = malloc(pix_col * sizeof(float));
            gcross_F_mat[i] = malloc(pix_col * sizeof(float));
            
            for (m = 0; m < pix_col; m++) {
                fscanf(gtan_F, "%f", &gtan_F_mat[i][m]);
            }
            for (m = 0; m < pix_col; m++) {
                fscanf(gcross_F, "%f", &gcross_F_mat[i][m]);
            }
        }
    
        fclose(gtan_F);
        fclose(gcross_F);
        
        float *gtan_F_arr = malloc(pix_tot * sizeof(float));
        float *gtan_F_arr2 = malloc(pix_tot * sizeof(float));
        float *gcross_F_arr = malloc(pix_tot * sizeof(float));
        
        j = 0;
        
        for (i = 0; i < pix_row; i++) {
            for (m = 0; m < pix_col; m++) {
                gtan_F_arr2[j] = gtan_F_mat[i][m];
                gcross_F_arr[j] = gcross_F_mat[i][m];
                j++;
            }
        }
        
        for (j = 0; j < pix_tot; j++) {
            gtan_F_arr[j] = gtan_F_arr2[j] - gtan_m_arr[j];
        }
        
        free(gtan_F_mat);
        free(gcross_F_mat);
        gtan_F_mat = NULL;
        gcross_F_mat = NULL;
        
        
        // x, y grids, the angles and radius to the center of the grids.
        float **x = malloc(pix_row * sizeof(float *));
        float **y = malloc(pix_row * sizeof(float *));
        float *db_ang = malloc(pix_tot * sizeof(float));
        float *r = malloc(pix_tot * sizeof(float));
        
        for (i = 0; i < pix_row; i++) {
            x[i] = malloc(pix_col * sizeof(float));
            y[i] = malloc(pix_col * sizeof(float));
        }
        
        for (i = 0; i < pix_row; i++) {
            for (m = 0; m < pix_col; m++) {
                x[i][m] = (-pix_row / 2 + m + 0.5) * delta_x;
                y[m][i] = (-pix_col / 2 + m + 0.5) * delta_y;
            }
        }
        
        j = 0;
        
        // Computing the angle squared of each grid.
        for (i = 0; i < pix_row; i++) {
            for (m = 0; m < pix_col; m++) {
                db_ang[j] = 2 * atan2(y[i][m], x[i][m]);
                r[j] = pow((pow(x[i][m], 2) + pow((y[i][m]), 2)), 0.5);
                j++;
            }
        }
        
        free(x);
        free(y);
        x = NULL;
        y = NULL;
        
        /////////////////////// Estimator /////////////////////////
        
        pixpair_count = 0;
        est_sum = 0;
        r_min = (min_ind - 0.5) * delta_x;
        r_max = (pix_row / 2 - 0.5) * delta_x;

        // Looping over all pixel pairs.
        #pragma omp parallel for default(shared) private(i,m,sin_db,cos_db) schedule(dynamic) reduction(+: est_sum, pixpair_count)
        for (i = 0; i < pix_tot; i++) {
            if (r[i] <= r_max && r[i] >= r_min) {
                for (m = 0; m < pix_tot; m++) { 
                    if (r[m] <= r_max && r[m] >= r_min) {
                        pixpair_count++;
                        // Not counting the same pair twice as it contributes fake signal.
                        if (m != i) {			   
                            sin_db = sin(db_ang[i] - db_ang[m]);
                            cos_db = cos(db_ang[i] - db_ang[m]);
                            est_sum += (gtan_F_arr[i] * gtan_F_arr[m] * gtan_a_arr[i] * gtan_a_arr[m] * cos_db + gtan_F_arr[i] * gcross_F_arr[m] * gtan_a_arr[i] * gcross_a_arr[m] * sin_db - gcross_F_arr[i] * gtan_F_arr[m] * gcross_a_arr[i] * gtan_a_arr[m] * sin_db + gcross_F_arr[i] * gcross_F_arr[m] * gcross_a_arr[i] * gcross_a_arr[m] * cos_db);
                        }
                    }
                }
            }
        }

        // Integral of shears.
        strcpy (multi_str, data_link);
        strcat (multi_str, "multi_exp/");
        strcat (multi_str, m_bin);
        strcat (multi_str, "/");
        strcat (multi_str, plane);
        strcat (multi_str, "/int_");
        strcat (multi_str, halo_na);
        strcat (multi_str, ".dat");
        
        multipole = fopen(multi_str, "r");
        
        if (!multipole) {
            printf("No multipole file %s found!\n", multi_str);
            exit(1);
        }
        
        multi_size = -1;
        multi_line = 0;
        
        float *r_bin = malloc(pix_row / 2 * sizeof(float));
        float *sig0 = malloc(pix_row / 2 * sizeof(float));
        float *int_int = malloc(pix_row / 2 * sizeof(float));
        float *ext_int = malloc(pix_row / 2 * sizeof(float));
        float *slopesig0 = malloc(pix_row / 2 * sizeof(float));
        float *gtan_int = malloc(pix_row / 2 * sizeof(float));
        float *gcross_int = malloc(pix_row / 2 * sizeof(float));
        float *est_int = malloc(pix_row / 2 * sizeof(float));
        
        for (m = 0; m < pix_row / 2; m++) {
            fscanf(multipole, "%f %f %*f %f %f %f", &r_bin[m], &sig0[m], &int_int[m], &ext_int[m], &slopesig0[m]);
            gtan_int[m] = (-slopesig0[m] + 3 * int_int[m] / pow(r_bin[m], 4) + ext_int[m]) / crit_den;
            gcross_int[m] = (3 * int_int[m] / pow(r_bin[m], 4) - ext_int[m]) / crit_den;
            est_int[m] = (pow(gtan_int[m], 2) + pow(gcross_int[m], 2)) * r_bin[m];
        }

        // Do the trapezoidal rule; N, number of panels.
        bin = r_bin[2] - r_bin[1];
        interval_num = (int)((r_max - r_min) / bin);        

        // x and f(x), with N + 1 intervals.
        float *integ = malloc((interval_num + 1) * sizeof(float));
        
        for (m = 0; m < interval_num + 1; m++) {
            integ[m] = est_int[m + min_ind - 1];
        }
        
        // Do the trapezoidal integral here.
        integral = 0;
        integral = trap_int(r_min, r_max, interval_num, integ);
        
        // Store up epis in the epis array.
        epis[k] = pow((pow((pow(r_max, 2) - pow(r_min, 2)), 2) * (est_sum / pixpair_count) / pow(1.837680E+09, 2)), 0.5);
        // Store up the stacked values.
        stack_est[k] = est_sum / pixpair_count;
        rmin[k] = r_min;
        rmax[k] = r_max;
        inter[k] = integral;
        
        fclose(multipole);
        
        free(integ);
        free(r);
        free(db_ang);
        free(gtan_F_arr);
        free(gtan_F_arr2);
        free(gcross_F_arr);
        free(r_bin);
        free(sig0);
        free(slopesig0);
        free(int_int);
        free(ext_int);
        free(gtan_int);
        free(gcross_int);
        integ = NULL;
        r = NULL;
        db_ang = NULL;
        gtan_F_arr = NULL;
        gtan_F_arr2 = NULL;
        gcross_F_arr = NULL;
        r_bin = NULL;
        sig0 = NULL;
        slopesig0 = NULL;
        int_int = NULL;
        ext_int = NULL;
        gtan_int = NULL;
        gcross_int = NULL;
    }
    
    strcpy (epis_str, data2_link);
    strcat (epis_str, "est_result/");
    strcat (epis_str, m_bin);
    strcat (epis_str, "/");
    strcat (epis_str, plane);
    strcat (epis_str, "_");
    strcat (epis_str, NN_str);
    //strcat (epis_str, "_test_estimate.dat");
    strcat (epis_str, "_est.dat");

    epis_file = fopen (epis_str, "w");
    
    if (!epis_file) {
        printf("Cannot open epison file!\n");
        exit(1);
    }
    
    for (k = NN * pile; k < halo_est; k++) {
	fprintf(epis_file, "%ld\t%f\t%E\t%E\t%E\t%E\n", halo_name[k], epis[k], inter[k], rmin[k], rmax[k], stack_est[k]);
    }
    
    
    free(epis);
    free(stack_est);
    free(rmin);
    free(rmax);
    free(inter);
    free(halo_name);
    free(gtan_m_arr);
    free(gtan_a_arr);
    free(gcross_a_arr);
    gtan_a_arr = NULL;
    gcross_a_arr = NULL;
    epis = NULL;
    stack_est = NULL;
    rmin = NULL;
    rmax = NULL;
    inter = NULL;
    halo_name = NULL;
    gtan_m_arr = NULL;
    
    fclose(epis_file);
    
    return 0;
}
