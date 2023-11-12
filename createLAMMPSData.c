#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct inputdata
{
	float x, y, z, vx, vy, vz;
	int sino, monomerID;
	char atomType[10], resIndex[10];
} GRO_DATA;

typedef struct outputdata
{
	int sino, atomType, molType, monomerID;
	float x, y, z, vx, vy, vz;
} LAMMPS_DATA;

typedef struct outputbonds
{
	int atom1, atom2, bondType, bondID;
} LAMMPS_BONDS;

typedef struct simulationboundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
} SIMULATION_BOUNDARY;


int getNumberOfLines (int nLines, FILE *input)
{
	char lineString[4000];
	nLines = 0;

	while (fgets (lineString, 4000, input) != NULL)
	{
		nLines++;
	}

	rewind (input);
	return nLines;
}

LAMMPS_BONDS *initOutputBonds (LAMMPS_BONDS *outputBonds, int nLines)
{
	for (int i = 0; i < nLines; ++i)
	{
		outputBonds[i].atom1 = 0;
		outputBonds[i].atom2 = 0;
		outputBonds[i].bondType = 0;
		outputBonds[i].bondID = 0;
	}

	return outputBonds;
}

LAMMPS_BONDS *generateBonds (LAMMPS_BONDS *outputBonds, LAMMPS_DATA *outputData, int nAtoms)
{
	int currentBond = 0;

	for (int i = 1; i < nAtoms; ++i)
	{
		if (
			((outputData[i].atomType == 1) && (outputData[i - 1].atomType == 2)) ||
			((outputData[i].atomType == 2) && (outputData[i - 1].atomType == 1))
			)
		{
			if ((outputData[i].monomerID == outputData[i - 1].monomerID) ||
				abs (outputData[i].monomerID - outputData[i - 1].monomerID) == 1
				)
			{
				outputBonds[currentBond].atom1 = i;
				outputBonds[currentBond].atom2 = i + 1;
				outputBonds[currentBond].bondType = 1;
				outputBonds[currentBond].bondID = currentBond + 1;
				currentBond++;
			}
		}

		if (outputData[i].atomType == 3 && outputData[i - 2].atomType == 1 ||
			outputData[i].atomType == 1 && outputData[i - 2].atomType == 3)
		{
			if (outputData[i].monomerID == outputData[i - 2].monomerID)
			{
				outputBonds[currentBond].atom1 = i + 1;
				outputBonds[currentBond].atom2 = i - 1;
				outputBonds[currentBond].bondType = 1;
				outputBonds[currentBond].bondID = currentBond + 1;
				currentBond++;
			}
		}
	}

	return outputBonds;
}

void printDatafile (LAMMPS_DATA *outputData, LAMMPS_BONDS *outputBonds, int nAtoms, int nBonds, FILE *output, SIMULATION_BOUNDARY boundaryInformation)
{
	fprintf(output, "%s\n\n", "LAMMPS data file converted from *.gro");
	fprintf(output, "%d atoms\n3 atom types\n%d bonds\n1 bond type\n\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nMasses\n\n1 13.019\n2 14.0270\n3 15.035\n\nAtoms\n\n", nAtoms, nBonds, boundaryInformation.xlo, boundaryInformation.xhi, boundaryInformation.ylo, boundaryInformation.yhi, boundaryInformation.zlo, boundaryInformation.zhi);

	for (int i = 0; i < nAtoms; ++i)
	{
		fprintf(output, "%d %d %d 0.0 %f %f %f\n", outputData[i].sino, outputData[i].molType, outputData[i].atomType, outputData[i].x, outputData[i].y, outputData[i].z);
	}

	fprintf(output, "\n\nBonds\n\n");

	for (int i = 0; i < nBonds; ++i)
	{
		fprintf(output, "%d %d %d %d\n", outputBonds[i].bondID, outputBonds[i].bondType, outputBonds[i].atom1, outputBonds[i].atom2);
	}
}

int countNBonds (LAMMPS_BONDS *outputBonds, int nLines, int nBonds)
{
	nBonds = 0;

	for (int i = 0; i < nLines; ++i)
	{
		if (outputBonds[i].atom1 > 0 && outputBonds[i].atom2 > 0)
		{
			nBonds++;
		}
	}

	return nBonds;
}

SIMULATION_BOUNDARY computeSimulationBoundary (SIMULATION_BOUNDARY boundaryInformation, LAMMPS_DATA *outputData, int nLines)
{
	boundaryInformation.xlo = outputData[0].x;
	boundaryInformation.xhi = outputData[0].x;
	boundaryInformation.ylo = outputData[0].y;
	boundaryInformation.yhi = outputData[0].y;
	boundaryInformation.zlo = outputData[0].z;
	boundaryInformation.zhi = outputData[0].z;

	for (int i = 0; i < nLines; ++i)
	{
		if (outputData[i].x < boundaryInformation.xlo)
		{
			boundaryInformation.xlo = outputData[i].x;
		}
		else if (outputData[i].x > boundaryInformation.xhi)
		{
			boundaryInformation.xhi = outputData[i].x;
		}

		if (outputData[i].y < boundaryInformation.ylo)
		{
			boundaryInformation.ylo = outputData[i].y;
		}
		else if (outputData[i].y > boundaryInformation.yhi)
		{
			boundaryInformation.yhi = outputData[i].y;
		}

		if (outputData[i].z < boundaryInformation.zlo)
		{
			boundaryInformation.zlo = outputData[i].z;
		}
		else if (outputData[i].z > boundaryInformation.zhi)
		{
			boundaryInformation.zhi = outputData[i].z;
		}
	}

	return boundaryInformation;
}

int main(int argc, char const *argv[])
{
	FILE *input, *output;
	input = fopen (argv[1], "r");
	output = fopen (argv[2], "w");

	char lineString[4000];

	int nLines = getNumberOfLines (nLines, input);

	GRO_DATA *inputData;
	inputData = (GRO_DATA *) malloc (nLines * sizeof (GRO_DATA));

	LAMMPS_DATA *outputData;
	outputData = (LAMMPS_DATA *) malloc (nLines * sizeof (LAMMPS_DATA));

	printf("Number of lines in the input file: %d\n", nLines);

	int currentLine = 0;

	while (fgets (lineString, 4000, input) != NULL)
	{
		sscanf (lineString, "%d%s %s %d %f %f %f %f %f %f\n", &inputData[currentLine].monomerID, &inputData[currentLine].resIndex, &inputData[currentLine].atomType, &inputData[currentLine].sino, &inputData[currentLine].x, &inputData[currentLine].y, &inputData[currentLine].z, &inputData[currentLine].vx, &inputData[currentLine].vy, &inputData[currentLine].vz);

		if (strstr (inputData[currentLine].atomType, "C1"))
		{
			outputData[currentLine].atomType = 1;
		}
		else if (strstr (inputData[currentLine].atomType, "C2"))
		{
			outputData[currentLine].atomType = 2;
		}
		else if (strstr (inputData[currentLine].atomType, "C3"))
		{
			outputData[currentLine].atomType = 3;
		}

		outputData[currentLine].sino = currentLine + 1;
		outputData[currentLine].x = inputData[currentLine].x * 10;
		outputData[currentLine].y = inputData[currentLine].y * 10;
		outputData[currentLine].z = inputData[currentLine].z * 10;
		outputData[currentLine].vx = inputData[currentLine].vx;
		outputData[currentLine].vy = inputData[currentLine].vy;
		outputData[currentLine].vz = inputData[currentLine].vz;
		outputData[currentLine].monomerID = inputData[currentLine].monomerID;

		if (currentLine < 7200) {
			outputData[currentLine].molType = 1; }
		else {
			outputData[currentLine].molType = 2; }

		// printf("%d %s %s %d %f %f %f %f %f %f\n", inputData[currentLine].monomerID, inputData[currentLine].resIndex, inputData[currentLine].atomType, inputData[currentLine].sino, inputData[currentLine].x, inputData[currentLine].y, inputData[currentLine].z, inputData[currentLine].vx, inputData[currentLine].vy, inputData[currentLine].vz);

		currentLine++;
	}

	int maxBonds = nLines;

	SIMULATION_BOUNDARY boundaryInformation;
	boundaryInformation = computeSimulationBoundary (boundaryInformation, outputData, nLines);

	printf("Simulation boundary:\n\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n", boundaryInformation.xlo, boundaryInformation.xhi, boundaryInformation.ylo, boundaryInformation.yhi, boundaryInformation.zlo, boundaryInformation.zhi);

	LAMMPS_BONDS *outputBonds;
	outputBonds = (LAMMPS_BONDS *) malloc (nLines * sizeof (LAMMPS_BONDS));

	outputBonds = initOutputBonds (outputBonds, nLines);
	outputBonds = generateBonds (outputBonds, outputData, nLines);
	int nBonds = countNBonds (outputBonds, nLines, nBonds);
	printf("Number of bonds: %d\n", nBonds);

/*	for (int i = 0; i < nLines; ++i)
	{//
		printf("%d %d %d %d\n", outputBonds[i].bondID, outputBonds[i].bondType, outputBonds[i].atom1, outputBonds[i].atom2);
		usleep (100000);
	}
*/

	printDatafile (outputData, outputBonds, nLines, nBonds, output, boundaryInformation);

	fclose (input);
	fclose (output);
	return 0;
}