#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#define CUTOFF 4.5
#define DIFF 4


void error (char *msg)
{
fprintf(stderr, "%s %s\n", msg, strerror(errno));
exit(0);
}

void UPDATE (char *str_point, char *ATOM, int *ATOM_NO, int *RESIDUE_NO, double *X, double *Y, double *Z)
{
size_t len = strlen (str_point);
unsigned int a = 0, b = 0;
char atom_no_array[5], residue_no_array[5], x_array[15], y_array[15], z_array[15];

 	for (a = 0; a < len; a++)       
	{
		if (a <= 5)
		{
	        ATOM[b] = str_point[a];
		b++;
		}
						
		if (a == 5)	
		{
		ATOM[b] = '\0';		
		b = 0;
		}

	       	if ((a >= 6) && (a <= 10))
		{
	        atom_no_array[b] = str_point[a];
		b++;
		}

		if (a == 10)
		{
		atom_no_array[b] = '\0';
		b = 0;
		}
	
	        if ((a >= 22) && (a <= 25))
		{
	        residue_no_array[b] = str_point[a];
		b++;
		}

		if (a == 25)
		{
		residue_no_array[b] = '\0';
		b = 0;
		}
	
	        if ((a >= 30) && (a <= 37))
		{
	        x_array[b] = str_point[a];
		b++;	
		}

		if (a == 37)
		{	
		x_array[b] = '\0';
		b = 0;
		}

	        if ((a >= 38) && (a <= 45))
		{
	        y_array[b] = str_point[a];
		b++;
		}

		if (a == 45)
		{
		y_array[b] = '\0';	
		b = 0;
		}
	
	 	if ((a >= 46) && (a <= 53))
		{
	       	z_array[b] = str_point[a];
		b++;
		}

		if (a == 53)
		{
		z_array[b] = '\0';
		b = 0;
		break;
		}
	}
*ATOM_NO = atoi (atom_no_array);	
*RESIDUE_NO = atoi (residue_no_array);
*X = atof (x_array);
*Y = atof (y_array);
*Z = atof (z_array);
}

int main (int argc, char *argv[])
{
char err_msg[100];
char strg[300];			/* Maximum no. of characters which can be read */
char atom[7];                   /* This are the first 6 characters in pdb */

int  atom_no = 0, residue_no = 0, last_element;
double atom_X = 0.0, atom_Y = 0.0, atom_Z = 0.0;
unsigned int chain = 0, line = 0, k = 0;
double diff_x, diff_y, diff_z, dist_calc;

int *atomno_arr = malloc (1*sizeof(int));
int *residueno_arr = malloc (1*sizeof(int));
unsigned int *MAP = malloc (1*sizeof(unsigned int));
double *X_arr = malloc (1*sizeof(double));
double *Y_arr = malloc (1*sizeof(double));
double *Z_arr = malloc (1*sizeof(double));

        if (argc < 2)
        {
        sprintf (err_msg,"Minimum 1 argument needed\n");
        error (err_msg);
        }
	
	if ((!atomno_arr) || (!residueno_arr) || (!X_arr) || (!Y_arr) || (!Z_arr))
	{
	sprintf (err_msg,"malloc error\n");
        error (err_msg);
	}
	
FILE *pdb = fopen (argv[1],"r");
        if (pdb == NULL)
        {
        sprintf (err_msg,"CANNOT OPEN %s\n",argv[1]);
        error (err_msg);
        }

	MAP[k] = line;
	while (fscanf(pdb, "%299[^\n]\n", strg) == 1)
	{
	UPDATE ( strg, atom, &atom_no, &residue_no, &atom_X, &atom_Y, &atom_Z);

		if (strstr(atom,"ATOM"))
		{
		atomno_arr = realloc (atomno_arr, (line+1)*sizeof (int));
		residueno_arr = realloc (residueno_arr, (line+1)*sizeof (int));
		X_arr = realloc (X_arr, (line+1)*sizeof (double));
		Y_arr = realloc (Y_arr, (line+1)*sizeof (double));
		Z_arr = realloc (Z_arr, (line+1)*sizeof (double));

		atomno_arr[line] = atom_no;
		residueno_arr[line] = residue_no;
		X_arr[line] = atom_X;		
		Y_arr[line] = atom_Y;
		Z_arr[line] = atom_Z;	

			if (residueno_arr[MAP[k]] != residueno_arr[line])
			{	
			k++;	
			MAP = realloc (MAP, (k+1)*sizeof (unsigned int));
			MAP[k] = line;
			}

		line++;
		}

		else if (strstr(atom,"TER"))
		{
		chain ++;
		unsigned int i, j, l, Q, m, c_alpha = 0;
		char AAname[20];
		char CAname[20];
		sprintf (AAname,"AA_%u.cont",chain);
		sprintf (CAname,"CA_%u.cont",chain);
		FILE *AA_file = fopen (AAname,"w");
		FILE *CA_file = fopen (CAname,"w");

                unsigned int **CONT_WEIGHT = (unsigned int **) malloc ((k - DIFF + 1) * sizeof (unsigned int *)); /* This creates an array of pointers */
                        if (!CONT_WEIGHT)
                                printf ("mem not allocated to rows\n");

                for (i = 0; i <= k - DIFF; i++)
                {
                CONT_WEIGHT[i] = (unsigned int *) malloc ((k - (i + DIFF) + 1) * sizeof (unsigned int));

                        if (!CONT_WEIGHT[i])
                                printf ("mem not allocated to columns");

                        for (j = i + DIFF; j <= k; j++)
                        {
                        CONT_WEIGHT[i][j - (i + DIFF)] = 0;
                        }
                }
		l = 0;
		Q = 0;
			for (i = MAP[0]; i <= MAP[k-DIFF]; i++)
			{
				if (i >= MAP[l+1])
					l++;
			
				m = l + DIFF;	
				for (j = MAP[l+DIFF]; j < line; j++)
				{
					if (j >= MAP[m+1])
						m++;	

				diff_x = X_arr[j] - X_arr[i];
        	                diff_y = Y_arr[j] - Y_arr[i];
                	        diff_z = Z_arr[j] - Z_arr[i];
                        	dist_calc = sqrt( diff_x * diff_x +
                                	          diff_y * diff_y +
                                            	  diff_z * diff_z );	
					if (dist_calc <= CUTOFF)	
					{	
					fprintf (AA_file, "%d %d\n",atomno_arr[i], atomno_arr[j]); /* This prints the atomic contacts */
/*					printf ("%d %d\n",residueno_arr[i], residueno_arr[j]);*/ /*This maps atoms to residue nos */
						if (CONT_WEIGHT[l][m - (l + DIFF)] == 0)
							c_alpha++;

					CONT_WEIGHT[l][m - (l + DIFF)] ++;
					Q++;
					}
				}	
			}

                for (i = 0; i <= k - DIFF; i++)
		{
                        for (j = i + DIFF; j <= k; j++)
                                if (CONT_WEIGHT[i][j - (i + DIFF)] > 0)
					fprintf (CA_file, "%d %d %lf\n",residueno_arr[MAP[i]], residueno_arr[MAP[j]], (double) (CONT_WEIGHT[i][j - (i + DIFF)]*c_alpha)/Q);
		free (CONT_WEIGHT[i]);
		}

		line = 0, k = 0;
		MAP[k] = line;	
                free (CONT_WEIGHT);
		fclose (AA_file);
		fclose (CA_file);
		}

		else if (strstr(atom,"END"))
			break;
	}
fclose (pdb);
free (atomno_arr);
free (residueno_arr);
free (X_arr);
free (Y_arr);
free (Z_arr);
free (MAP);
return 0;
}
