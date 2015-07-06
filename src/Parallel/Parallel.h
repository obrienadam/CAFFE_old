#ifndef PARALLEL_H
#define PARALLEL_H

#include "Array3D.h"
#include "Matrix.h"

class Parallel
{

public:

static int nProcesses();
static int processNo();

static int sum(int number);
static int min(int number);
static int max(int number);

static double sum(double number);
static double min(double number);
static double max(double number);

static void sendMatrix(int root, int dest, Matrix& matrix);

};

#endif
