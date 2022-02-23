#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
						
#define RANDMAX 32767	
int H[1472][2176];
int RAW_mat[46][68];
int Exchan_mat[2176][1472];

float rand49()
{ 
	static int Num = 0;
	double number;
	int i;
	i = rand() % (RANDMAX + 1);
	number = (double)i / ((unsigned)(RANDMAX + 1));
	Num++;
	if (Num >= RANDMAX) {
		srand((unsigned)(time(NULL) % RAND_MAX));
		Num = 0;
	}
	return (float)number;
}

double Normal() 
{
	static int iset = 0;
	static double qset;
	double vx, vy, r, temp;
	if (iset == 0)
	{
		do
		{
			vx = 2.0 * rand49() - 1.0; 
			vy = 2.0 * rand49() - 1.0; 
			r = vx * vx + vy * vy;
		} while (r >= 1.0 || r == 0);
		temp = sqrt(-2.0 * log(r) / r);
		qset = vy * temp;
		iset = 1;
		return (vx * temp);
	}
	else
	{
		iset = 0;
		return qset;
	}
}

double sign(double x)
{
	if (x >= 0)	return 1.0;
	else	return -1.0;
}

int main (int argc,char *argv[]) // xinYing:we can also define the value in the below equation...because the io stream is not always work
{
	int SNR_top = atoi(argv[1]);
	int inter_num = atoi(argv[2]);
	int Erroer_thrshold = atoi(argv[3]);
	double attenuation = atof(argv[4]);
	printf("-------------------------------------------------------------------------------------------\n");
	printf("|SNR from 1.0 to %.2f\titeration : %d\t Error threshold : %d\t Attenuation : %f\t  |\n", SNR_top*0.25+0.75, inter_num, Erroer_thrshold, attenuation);
    int Mes[704];
    int Enc[2176], Dec[2176], mult[1472];
	
	//LLR matrix
	double **L;
	L = (double **)malloc(1472 * sizeof(double *));
	int i;
    for (i=0; i<1472; i++) 

         L[i] = (double *)malloc(2176 * sizeof(double)); 

	//read file
	FILE* Hmtx;
	int j, k, l;

	Hmtx = fopen("H_matrix.txt", "r");
	for (i = 0; i < 46; i++)	// read H matrix
		for (j = 0; j < 68; j++)
			fscanf(Hmtx, "%d", &RAW_mat[i][j]);
	fclose(Hmtx);
	for (i = 0; i < 46; i++) {	// lifting
		for (j = 0; j < 68; j++) {
			if (RAW_mat[i][j] == -1)
				for (k = 0; k < 32; k++) {
					for (l = 0; l < 32; l++) {
						H[i * 32 + k][j * 32 + l] = 0;
					}
				}
			else
				for (k = 0; k < 32; k++) {
					for (l = 0; l < 32; l++) {
						if (l == ((k + RAW_mat[i][j]) % 32))
							H[i * 32 + k][j * 32 + l] = 1;
						else
							H[i * 32 + k][j * 32 + l] = 0;
					}
				}
		}
	}
    
	//init Exchan_mat
	for (i=0; i<2176; i++)
		for (j=0; j<1472; j++)
			Exchan_mat[i][j]=H[j][i];

	printf("-------------------------------------------------------------------------------------------\n");
	printf("|SNR\t\t| Frame_error\t| Total_frame\t| Bit_error\t| BER\t\t| FER\t  |\n");
	printf("-------------------------------------------------------------------------------------------\n");
	double SNR_start = 1.0;
	double step = 0.25;
	int iteration, snr, error,F_er ,error_tmp;
	
	int num;
    int frame;
	int Decr_stop = 1;
	int P_bit[2176 - 704];
	int temp[4 * 32];
	

	double sig[2176], prior[2176];
	double smallest, smallest2, SIG_sign;
	double R = 704.0 / 2176.0;    //coderate
	double var, sigma, SNR, BER, FER;
    for (snr=0; snr<SNR_top; snr++) { 
        SNR = SNR_start + snr*step;
        error = 0;
        error_tmp=0;
        frame=0;
        F_er =0;
		
        var = 1.0/(2.0*R*pow(10,SNR/10));
        sigma=sqrt(var); 


		do{
            Decr_stop=1;
			
		
            for (i=0;i<704;i++) { 
                Mes[i] = rand()%2;
            }
            
            //LDPC Encoder
			for (i=0; i<4*32; i++) {
				temp[i] = 0;
				for (j=0; j<704; j++)
					temp[i] += H[i][j]*Mes[j];
				temp[i] = temp[i]%2;
			}
			for (i=0; i<32; i++) {
				P_bit[i] = 0;
				for (j=0; j<4; j++)
					P_bit[i] += temp[j*32+i];
				P_bit[i] = P_bit[i]%2;
			}
			for (i=0; i<32; i++)
				P_bit[32+i] = (temp[i] + P_bit[(i+1)%32])%2;
			for (i=0; i<32; i++)
				P_bit[2*32+i] = (temp[32+i] + P_bit[i] + P_bit[32+i])%2;
			for (i=0; i<32; i++)
				P_bit[3*32+i] = (temp[2*32+i] + P_bit[2*32+i])%2;
			for (i=0; i<2176; i++)
				if (i<704)	Enc[i] = Mes[i];
				else if (i<704+4*32)	Enc[i] = P_bit[i-704];
				else	Enc[i] = 0;
			for (i=4*32; i<1472; i++) {
				Enc[704+i] = 0;
				for (j=0; j<704+4*32; j++)
					Enc[704+i] += H[i][j]*Enc[j];
				Enc[704+i] = Enc[704+i]%2;
			}		
			
			// BPSK modulation and AWGN channel
            for (i=0; i<2176; i++) {   
                sig[i] = sigma*Normal()-2*Enc[i]+1; // noise
            }

            for (i=0; i<2176; i++) {   // LLR value
            	prior[i] = 2*sig[i]/(sigma*sigma);
            }


			// LDPC Decoder
			for (i=0; i<1472; i++)
				for (j=0; j<2176; j++)
					if (H[i][j]==1)
						L[i][j]=0;
					
			for (iteration=0; iteration<inter_num; iteration++) {
				for (i=0; i<1472; i++) {
					SIG_sign = 1.0;
					for (j = 0; j < 2176; j++) {
						if (H[i][j] == 1) {
							L[i][j] = prior[j] - L[i][j];	
							SIG_sign *= sign(L[i][j]);//find out the sign of whole line after multiply
						}
					}

					//find the 1st and 2nd smallest element of a line
					smallest = 1000000000;   //set a big number
					smallest2 = 1000000000;
					
			    	for (j=0; j<2176; j++) {	
			    		if (H[i][j]==1) {
							if (fabs(L[i][j]) < smallest) {
								smallest2 = smallest;
								smallest = fabs(L[i][j]);
							}
							else if (fabs(L[i][j]) < smallest2) {
								smallest2 = fabs(L[i][j]);
							}
			    				
						}
					}

			    	for (j=0; j<2176; j++) {
			    		if (H[i][j]==1) {
			    			if (fabs(L[i][j])==smallest) {//each element of a line need to puduct a smallest fabs except itself
								prior[j]=L[i][j]+attenuation*sign(L[i][j])*SIG_sign*smallest2;
								L[i][j]=attenuation*sign(L[i][j])*SIG_sign*smallest2; 
							}
							else {
								prior[j]=L[i][j]+attenuation*sign(L[i][j])*SIG_sign*smallest;
								L[i][j]=attenuation*sign(L[i][j])*SIG_sign*smallest; 						
							}
						}
					}
				}
				for (j=0; j<2176; j++) {
					if (prior[j]>0)	Dec[j]=0;
					else	Dec[j]=1;
				}	
				for (i=0; i<1472; i++) {	
					mult[i]=0;
				    for (j=0; j<2176; j++)
						mult[i]+=Dec[j]*Exchan_mat[j][i];
					mult[i]=mult[i]%2;
				}
				num=0;
				for (i=0; i<1472; i++)	num+=mult[i];
				if (num==0)	iteration=inter_num;	
			}
    		for (i=0; i<704; i++) {
    			if (Mes[i]!=Dec[i])	error++;
    		}
  
            if (error_tmp != error) {
                F_er++;
				error_tmp=error;
			}
            frame++;
            
			printf("\r|%f\t| %d\t\t| %d\t\t| %d\t\t| ", SNR, F_er, frame, error);


            if (F_er==Erroer_thrshold)
                Decr_stop=0;
        } while(Decr_stop == 1);
		
        BER = (double)error/(double)(frame*704);
        FER = (double)F_er/(double)frame;
		
        printf("%.2E\t|%.2E|\n",BER,FER);
		printf("-------------------------------------------------------------------------------------------\n");
    } 
	
	free(L);

}


