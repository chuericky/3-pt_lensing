CC = gcc
CFLAGS = -O3
HOME_DIR = /Users/rickyccy/Documents/Research_weak_lensing/Data/final_check_Feb_25_2014/
LIBS = -lm

all: sur_den_project multipole_expansion estimator rotation halo_select


sur_den_project: sur_den_project.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o sur_den_project sur_den_project.c tool.c $(LIBS) -fopenmp

rotation: rotation.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o rotation rotation.c tool.c $(LIBS)

multipole_expansion: multipole_expansion.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o multipole_expansion multipole_expansion.c tool.c $(LIBS)

estimator: estimator.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o estimator estimator.c tool.c $(LIBS) -fopenmp

halo_select:	halo_select.c Makefile
	$(CC) $(CFLAGS) -o halo_select halo_select.c $(LIBS)


clean:
	rm *.o
