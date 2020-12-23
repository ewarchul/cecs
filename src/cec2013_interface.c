#include "cec2013_interface.h"
#include <string.h>

void cec2013_func(double *x, double *f, int nx, int mx,int func_num)
{

	int cf_num=10,i,ret;
	if (ini_flag==1)
	{
		if ((n_flag!=nx)||(func_flag!=func_num))
		{
			ini_flag=0;
		}
	}

	if (ini_flag==0)
	{
		FILE *fpt;
		char FileName[256];
		free(M);
		free(OShift);
		free(y);
		free(z);
		free(x_bound);
		y=(double *)malloc(sizeof(double)  *  nx);
		z=(double *)malloc(sizeof(double)  *  nx);
		x_bound=(double *)malloc(sizeof(double)  *  nx);
		for (i=0; i<nx; i++)
			x_bound[i]=100.0;
		
		sprintf(FileName, "%s/M_D%d.txt", extdata, nx);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
		    perror("Cannot open input file for reading");
		}

		M=(double*)malloc(cf_num*nx*nx*sizeof(double));
		if (M==NULL)
			perror("There is insufficient memory available!");
		for (i=0; i<cf_num*nx*nx; i++)
		{
				ret = fscanf(fpt,"%lf",&M[i]);
				if (ret != 1)
				{
				    perror("Error reading from the input file");
				}
		}
		fclose(fpt);
		

		sprintf(FileName, "%s/shift_data.txt", extdata);
		fpt=fopen(FileName,"r");
		if (fpt==NULL)
		{
			perror("Cannot open input file for reading");
		}
		OShift=(double *)malloc(nx*cf_num*sizeof(double));
		if (OShift==NULL)
			perror("There is insufficient memory available!");
		for(i=0;i<cf_num*nx;i++)
		{
				ret = fscanf(fpt,"%lf",&OShift[i]);
				if (ret != 1)
				{
				    perror("Error reading from the input file");
				}
		}
		fclose(fpt);

		n_flag=nx;
		func_flag=func_num;
		ini_flag=1;
	}

	for (i = 0; i < mx; i++)
	{
		switch(func_num)
		{
		case 1:	
			cec2013_sphere_func(&x[i*nx],&f[i],nx,OShift,M,0);
			f[i]+=-1400.0;
			break;
		case 2:	
			cec2013_ellips_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-1300.0;
			break;
		case 3:	
			cec2013_bent_cigar_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-1200.0;
			break;
		case 4:	
			cec2013_discus_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-1100.0;
			break;
		case 5:
			cec2013_dif_powers_func(&x[i*nx],&f[i],nx,OShift,M,0);
			f[i]+=-1000.0;
			break;
		case 6:
			cec2013_rosenbrock_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-900.0;
			break;
		case 7:	
			cec2013_schaffer_F7_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-800.0;
			break;
		case 8:	
			cec2013_ackley_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-700.0;
			break;
		case 9:	
			cec2013_weierstrass_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-600.0;
			break;
		case 10:	
			cec2013_griewank_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-500.0;
			break;
		case 11:	
			cec2013_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,0);
			f[i]+=-400.0;
			break;
		case 12:	
			cec2013_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-300.0;
			break;
		case 13:	
			cec2013_step_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=-200.0;
			break;
		case 14:	
			cec2013_schwefel_func(&x[i*nx],&f[i],nx,OShift,M,0);
			f[i]+=-100.0;
			break;
		case 15:	
			cec2013_schwefel_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=100.0;
			break;
		case 16:	
			cec2013_katsuura_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=200.0;
			break;
		case 17:	
			cec2013_bi_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,0);
			f[i]+=300.0;
			break;
		case 18:	
			cec2013_bi_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=400.0;
			break;
		case 19:	
			cec2013_grie_rosen_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=500.0;
			break;
		case 20:	
			cec2013_escaffer6_func(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=600.0;
			break;
		case 21:	
			cec2013_cf01(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=700.0;
			break;
		case 22:	
			cec2013_cf02(&x[i*nx],&f[i],nx,OShift,M,0);
			f[i]+=800.0;
			break;
		case 23:	
			cec2013_cf03(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=900.0;
			break;
		case 24:	
			cec2013_cf04(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=1000.0;
			break;
		case 25:	
			cec2013_cf05(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=1100.0;
			break;
		case 26:
			cec2013_cf06(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=1200.0;
			break;
		case 27:
			cec2013_cf07(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=1300.0;
			break;
		case 28:
			cec2013_cf08(&x[i*nx],&f[i],nx,OShift,M,1);
			f[i]+=1400.0;
			break;
		}
		
	}


}