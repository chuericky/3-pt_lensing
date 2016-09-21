// This program selects the interested halos and write the relevant information in a text file. All distance in Mpc /h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define z_red 0.5       // Redshift of the lens
#define delta_x 0.1     // size of the grid
#define delta_y 0.1
#define delta_z 0.1


// A struct to copy the file names.
typedef struct {
    long name_of_file;
} fname;

typedef struct {
    long int filesname;
    float x0, y0, z0, Rvir, Mvir;
    int partnum;
} filescheck;


int main(int argc, char *argv[])
{
    FILE *file_range, *file_work, *outputfile;
    int spot_range, j, indices, k, checkfile = 0, m = 0, range_index = 0, range_count = 0, progress = 1;
    int name[10];
    char label[range_index][11];
    
    // Input the range of the halos
    file_range = fopen("/Users/rickyccy/Documents/Research_weak_lensing/Data/multipole_model/analysis_sigma0_sigma2/halo_catalog_0.2_11.9.dat", "r");
    
    if (!file_range) {
        printf("Open file failed!\n");
        exit(1);
    }
    
    // Count the number of files.
    do {
        spot_range = fgetc(file_range);
        if(spot_range == '\n')
            range_count++;
    } while (spot_range != EOF);
    
    // last line doesn't end with a new line!
    // but there has to be a line at least before the last line
    if(spot_range != '\n' && range_count != 0)
        range_count++;
    
    rewind(file_range);
    
    fname *fnames = malloc(range_count * sizeof *fnames);
    
    // Get the names of the files.
    while(fscanf(file_range, "%ld", &fnames[range_index].name_of_file) == 1)
        range_index++;

    // Convert the file names from long int to strings.
    for (j = 0; j < range_index; j++) {
        sprintf(label[j], "%ld", fnames[j].name_of_file);
        if (label[j][4] != label[j - 1][4]) {
            name[m] = j;
            m++;
        }
    }
    
    name[m] = range_index;
    
    fclose(file_range);
    
    outputfile = fopen("/Users/rickyccy/Documents/Research_weak_lensing/Data/shear_estimator_test/halo_catalog_0.2_11.9.dat", "w");

    // Open the files for 
    for (j = 0; j < m; j++) {
        char str[150], index[2], index_1[2];
        int spot_file, individual_count = 0;
        strcpy (str, "/Users/rickyccy/Documents/Research_weak_lensing/Data/iso_4160");
        sprintf(index, "%d", j);
        strcat (str, index);
        strcat (str, "000001_4160");
        sprintf(index_1, "%d", j+1);
        strcat (str, index_1);
        strcat (str, "000000.txt");
        printf("%s\n", str);
        file_work = fopen(str, "r");
        
        if (!file_work) {
            printf("Open file failed!\n");
            exit(1);
        }
        
        do {
            spot_file = fgetc(file_work);
            if(spot_file == '\n')
                individual_count++;
        } while (spot_file != EOF);
        
        // last line doesn't end with a new line!
        // but there has to be a line at least before the last line
        if(spot_file != '\n' && individual_count != 0)
            individual_count++;
        
        rewind(file_work);
        
        filescheck *fch = malloc(individual_count * sizeof *fch);
        
        indices = 0;
        
        // Get the names of the files.
        while(fscanf(file_work, "%ld %*f %*f %*f %f %f %f %*f %*f %*f %d %*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f", &fch[indices].filesname, &fch[indices].x0, &fch[indices].y0, &fch[indices].z0, &fch[indices].partnum, &fch[indices].Mvir, &fch[indices].Rvir) == 7) {

            if (fnames[checkfile].name_of_file == fch[indices].filesname) {
                fprintf(outputfile, "%ld\t%.4f\t%.4f\t%.4f\t%d\t%E\t%.4f\n", fch[indices].filesname, fch[indices].x0, fch[indices].y0, fch[indices].z0, fch[indices].partnum, fch[indices].Mvir, fch[indices].Rvir);
                checkfile++;
            }
            indices++;
        }
        
        printf("%d out of %d files completed!\n", progress, m);
        printf("%d\n", checkfile);
        progress++;
        
        fclose(file_work);
    }
    
    fclose(outputfile);
    return 0;
}
