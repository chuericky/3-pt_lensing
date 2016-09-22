// Updated: Feb 9, 2014

#include "tool.h"

// A function to count number of lines in a file.
int count_line(int spot, FILE *countfile, int line) {
    do {
        spot = fgetc(countfile);
        if(spot == '\n')
            line++;
    } while (spot != EOF);
    
    // last line doesn't end with a new line!
    // but there has to be a line at least before the last line
    if(spot != '\n' && line != 0)
        line++;
    
    if (line == 0)
        line++;
    
    rewind(countfile);
    
    return line;
}

// A function to find the maximum value of the difference of the array elements to the subtraction.
int max(float arr[], int size, float subtract) {
    float maximum = fabs(arr[0] - subtract);
    int j, count;
    for (j = 0; j < size; j++) {
        if (maximum < fabs(arr[j] - subtract)) {
            maximum = fabs(arr[j] - subtract);
            count = j;
        }
    }
    return count;
}

// A function to find the minimum value of the difference of the array elements to the subtraction.
int min(float arr[], int size, float subtract) {
    float minimum = fabs(arr[0] - subtract);
    int j, count;
    for (j = 0; j < size; j++) {
        if (minimum > fabs(arr[j] - subtract)) {
            minimum = fabs(arr[j] - subtract);
            count = j;
        }
    }
    return count;
}

// A function to find the maximum value of the index.
float max_index(float arr[], int size) {
    float maximum = arr[0];
    int j, count;
    for (j = 0; j < size; j++) {
        if (maximum < arr[j])
            maximum = arr[j];
    }
    return maximum;
}

// A function to find the minimum value of the index.
float min_index(float arr[], int size) {
    float minimum = arr[0];
    int j, count;
    for (j = 0; j < size; j++) {
        if (minimum > arr[j])
            minimum = arr[j];
    }
    return minimum;
}

// A function for trapezoidal rule. Interval_num should be less than the length of f by 1.
float trap_int(float a, float b, int interval_num, float f[]) {
    float sum = 0, result = 0;
    int i;
    for (i = 1; i < interval_num; i++) {
        sum += 2 * f[i];
    }
    sum += f[0];
    sum += f[interval_num];
    result = ((b - a) * sum) / (2. * interval_num);
    return result;
}

// A function to skip certain lines for a file.
void skip_line(int spot, FILE *readfile, int line, int skipped_line) {
    do {
        spot = fgetc(readfile);
        if(spot == '\n')
            line++;
        if (line == skipped_line)
            break;
    } while (spot != EOF);
    
    return;
}


