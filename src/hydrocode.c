/* 
 * This is a implementation of 1-D or 2-D hydrocodes in Cartesian geometry
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <dirent.h>

#include "./include/var_struc.h"
#include "./include/tools.h"
#include "./include/file_io.h"
#include "./include/meshing.h"
#include "./include/finite_volume.h"



/*The global primitive variable and configuration array.
 */

double config[N_CONF];



int main(int argc, char *argv[])
{
	printf("\n%s\n", argv[2]);

	for (int i = 0; i < N_CONF; i++)
			config[i] = 1.0/0.0;

	// Set dimension.
	config[0] = (double)atoi(argv[3]);

	char * scheme;
	config[9] = (double)strtol(argv[4], &scheme, 10);	
	if (* scheme == '_')
		scheme++;
	else
		{
			printf("No order!\n");
			exit(2);
		}

	struct flu_var FV = flu_conf_load(argv[1]);

	struct mesh_var mv= mesh_load(argv[1], argv[5]);
	

	printf("Output: %s.\n", argv[2]);


	file_write_TEC(FV, mv, argv[2], 0.0, (int)config[0]);	

	
	Euler_scheme(&FV, mv, scheme);

	file_write_TEC(FV, mv, argv[2], config[1], (int)config[0]);

	if((int)config[0] > 1)
		file_write_VTK_3D(FV, mv, argv[2]);
	
	return 0;	
}
