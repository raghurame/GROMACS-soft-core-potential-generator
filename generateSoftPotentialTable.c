#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define R_BIN_DIST 0.01
#define N_POT_LINES 251
#define DR_FOR_DERIVATIVE 0.001

typedef struct softPotential
{
	float r, f, df, g, dg, h, dh;
} POTENTIAL;

float calculateSigma (float v, float w)
{
	float sigma;

	sigma = w / v;
	sigma = pow (sigma, 0.166);

	return sigma;
}

float calculateEpsilon (float v, float w, float sigma)
{
	float epsilon;

	epsilon = v / (4 * pow (sigma, 6));

	return epsilon;
}

float *initPotEnergy (float *repulsivePotEnergy, int arrayLength)
{
	for (int i = 0; i < arrayLength; ++i)
	{
		repulsivePotEnergy[i] = 0;
	}

	return repulsivePotEnergy;
}

float computeRepulsiveEnergy (int n, float lambda, float alphaLJ, float v, float w, float sigma, float epsilon, float r)
{
	float repulsivePotEnergy;

	repulsivePotEnergy = pow (lambda, n) * 4 * epsilon;
	repulsivePotEnergy /= pow ((alphaLJ * pow ((1 - lambda), 2) + pow ((r / sigma), 6)), 2);

	return repulsivePotEnergy;
}

float computeAttractiveEnergy (int n, float lambda, float alphaLJ, float v, float w, float sigma, float epsilon, float r)
{
	float attractivePotEnergy;

	// attractivePotEnergy = v / pow (r, 6);
	attractivePotEnergy = pow (lambda, n) * 4 * epsilon;
	attractivePotEnergy /= alphaLJ * pow ((1 - lambda), 2) + pow ((r / sigma), 6);

	return attractivePotEnergy;
}

float computeDerivativeRepulsiveEnergy (int n, float lambda, float alphaLJ, float v, float w, float sigma, float epsilon, float r, float dr)
{
	float derivative, y1, y2, y3, x1, x2, x3, dx1, dx2, dy1, dy2;

	x1 = r - dr;
	x2 = r;
	x3 = r + dr;

	y1 = computeRepulsiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, x1);
	y2 = computeRepulsiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, x2);
	y3 = computeRepulsiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, x3);

	dx1 = x2 - x1;
	dx2 = x3 - x2;

	dy1 = y2 - y1;
	dy2 = y3 - y2;

	derivative = (dy1 / dx1) + (dy2 / dx2);
	derivative /= 2;

	return derivative;
}

float computeDerivativeAttractiveEnergy (int n, float lambda, float alphaLJ, float v, float w, float sigma, float epsilon, float r, float dr)
{
	float derivative, x1, x2, x3, y1, y2, y3, dx1, dx2, dy1, dy2;

	x1 = r - dr;
	x2 = r;
	x3 = r + dr;

	y1 = computeAttractiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, x1);
	y2 = computeAttractiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, x2);
	y3 = computeAttractiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, x3);

	dx1 = x2 - x1;
	dx2 = x3 - x2;

	dy1 = y2 - y1;
	dy2 = y3 - y2;

	derivative = (dy1 / dx1) + (dy2 / dx2);
	derivative /= 2;

	return derivative;
}

POTENTIAL *computeSoftPotentialValues (POTENTIAL *potValues, int arrayLength, int n, float lambda, float alphaLJ, float v, float w, float sigma, float epsilon)
{
	float *repulsivePotEnergy;
	repulsivePotEnergy = (float *) malloc (arrayLength * sizeof (float));

	repulsivePotEnergy = initPotEnergy (repulsivePotEnergy, arrayLength);

	for (int i = 0; i < arrayLength; ++i)
	{
		if (i == 0)
		{
			potValues[i].r = 0;
		}
		else
		{
			potValues[i].r = potValues[i - 1].r + R_BIN_DIST;
		}
	}

	potValues[0].r = R_BIN_DIST / 1000;

	for (int i = 0; i < arrayLength; ++i)
	{
		repulsivePotEnergy[i] = computeRepulsiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, potValues[i].r);
	}

	for (int i = 0; i < arrayLength; ++i)
	{
		potValues[i].h = repulsivePotEnergy[i] / (4 * epsilon * pow (sigma, 12));
	}

	for (int i = 0; i < arrayLength; ++i)
	{
		potValues[i].dh = computeDerivativeRepulsiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, potValues[i].r, DR_FOR_DERIVATIVE);
	}

	for (int i = 0; i < arrayLength; ++i)
	{
		potValues[i].g = computeAttractiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, potValues[i].r);
		potValues[i].g /= (4 * epsilon * pow (sigma, 6));
		potValues[i].dg = computeDerivativeAttractiveEnergy (n, lambda, alphaLJ, v, w, sigma, epsilon, potValues[i].r, DR_FOR_DERIVATIVE);
	}

	return potValues;
}

void printPotValues (FILE *output_file, POTENTIAL *potValues, int arrayLength, float sigma, float epsilon)
{
	potValues[0].r = 0;

	for (int i = 0; i < arrayLength; ++i)
	{
		if (potValues[i].r < sigma)
		{
			fprintf(output_file, "%f\t%d\t%d\t%4E\t%4E\t%4E\t%4E\n", potValues[i].r, 0, 0, potValues[i].g, -potValues[i].dg, potValues[i].h, -potValues[i].dh);
		}
		else
		{
			fprintf(output_file, "%f\t%d\t%d\t%4E\t%4E\t%4E\t%4E\n", potValues[i].r, 0, 0, 0, 0, 0, 0);
		}
	}
}

int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		printf("\nERROR: INCORRECT ARGUMENTS PASSED\n\n REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = ./program\n {~} argv[1] = output filename\n {~} argv[2] = (soft potential variable) n (= 1 or 2)\n {~} argv[3] = (soft potential variable) lambda (= 0 to 1)\n {~} argv[4] = (soft potential variable) alphaLJ (= 0.5)\n {~} argv[5] = (soft potential variable) v\n {~} argv[6] = (soft potential variable) w\n\n");
		exit (1);
	}
	FILE *output_file;
	output_file = fopen (argv[1], "w");

	float n = atof (argv[2]), lambda = atof (argv[3]), alphaLJ = atof (argv[4]), v = atof (argv[5]), w = atof (argv[6]);

	float sigma = calculateSigma (v, w);
	float epsilon = calculateEpsilon (v, w, sigma);

	printf("sigma: %2E\nepsilon: %2E\n", sigma, epsilon);

	POTENTIAL *potValues;
	potValues = (POTENTIAL *) malloc (N_POT_LINES * sizeof (POTENTIAL));

	potValues = computeSoftPotentialValues (potValues, N_POT_LINES, n, lambda, alphaLJ, v, w, sigma, epsilon);
	printPotValues (output_file, potValues, N_POT_LINES, sigma, epsilon);

	fclose (output_file);
	return 0;
}