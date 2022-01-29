#include <iostream>
#include <cmath>

using namespace std;
int numOfIter = 0;
double funVal(double x) {
	return cos(x) / powl((1 - cos(x)), 2);
}

double IntegralResultFun(double x) {
	return 1.0 / 2.0 * (1.0 / tan(x / 2.0)) - 1.0 / 6.0 * (powl((1.0 / tan(x / 2.0)), 3.0));
}

double GetRealIntegralValue(double a, double b) {
	double result = IntegralResultFun(b) - IntegralResultFun(a);
	return result;
}

double testFun(double x) {
	return x;
}

double *GetEvenGrid(int numOfDots, double a, double b, double *delta) {
	*delta = (b - a) / (double)numOfDots;
	double *grid = new double[numOfDots - 1];
	grid[0] = a + *delta;
	for (int i = 1; i < numOfDots - 1; i++) {
		grid[i] = grid[i - 1] + (*delta);
	}
	return grid;
}

double GetSimpsonValue(int N, double a, double b, double(*funVal) (double)) {
	double S;
	double S1 = funVal(a) + funVal(b);
	double S2 = 0;////sum over even nodes
	double S3 = 0;////sum over odd nodes
	double delta;
	double I1 = 0;

	double *nodes = GetEvenGrid(N, a, b, &delta);

	for (int i = 1; i < N - 2; i += 2) {
		S2 += funVal(nodes[i]);
	}

	//S3 count
	for (int i = 0; i < N - 1; i += 2) {
		S3 += funVal(nodes[i]);
	}

	//calculation of the total
	S = S1 + 4 * S3 + 2 * S2;

	//we get the value of the integral
	I1 = S * (delta / 3);
	
	return I1;
}

double IntegralValSimp(double a, double b, double e, double (*funVal) (double), double factValue, double d) {
	int n = 4;//number of partition segments
	double delta;

	double I1 = 0, I2 = INFINITY;
	
	double S;
	double S1 = funVal(a) + funVal(b);
	double S2 = 0;//sum over even nodes
	double S3 = 0;//sum over odd nodes
	bool isReady = false;
	bool isFirstIteration = true;
	numOfIter = 0;

	while (abs(I1 - I2) > 15 * e) {
		numOfIter++;
		if (I2 != INFINITY) {
			I1 = I2;
			n *= 2;
		}
		I2 = 0;

		//get an array of nodes
		double *nodes = GetEvenGrid(n, a, b, &delta);
		
		//S2 count
		if (isFirstIteration) {
			for (int i = 1; i < n - 2; i += 2) {
				S2 += funVal(nodes[i]);
			}
			isFirstIteration = false;
		}
		
		//S3 count
		for (int i = 0; i < n - 1; i += 2) {
			S3 += funVal(nodes[i]);
		}

		//calculation of the total
		S = S1 + 4 * S3 + 2 * S2;

		//we get the value of the integral
        I2 = S * ( delta / 3 );

		//economically saving the calculated result
		S2 = S3 + S2;

		S3 = 0;
		if (abs(factValue - I2) < d && !isReady) {
			printf("%lf\n", I2);
			printf("%lf\n", d / factValue);
			printf("%i\n", n);
			isReady = true;
		}
	}
	printf("%i\n", n);
	return I2;
}

double ElementaryChebIntegral(double a, double b, double(*funVal) (double), int n, double *masOfNodes) {
	double coef1 = (b - a) / 2;
	double coef2 = (b + a) / 2;
	double I = 0;

	for (int i = 0; i < n; i++) {
		I += funVal(coef2 + coef1 * masOfNodes[i]);
	}

	I *= (double)2 / (double)3;
	I *= coef1;

	return I;
}

double *GetGridForCheb(int numOfDots, double a, double b, double *delta) {
	*delta = (b - a) / (double)numOfDots;
	double *grid = new double[numOfDots + 1];
	grid[0] = a;
	for (int i = 1; i < numOfDots + 1; i++) {
		grid[i] = grid[i - 1] + (*delta);
	}
	return grid;
}

double NumericChebIntegration(double a, double b, double e, double(*funVal) (double), double *masOdNodes, int n, double factValue, double d) {
	int numOfIntervals = 1;
	double delta;
	double *masOfIntervals;
	double integralSum = 0;
	double prevSum = INFINITY;
	bool isReady = false;
	while (abs(prevSum - integralSum) > e) {
		if (integralSum != 0) {
			prevSum = integralSum;
			integralSum = 0;
			numOfIntervals += 1;
		}
		masOfIntervals = GetGridForCheb(numOfIntervals, a, b, &delta);
		for (int i = 0; i < numOfIntervals; i++) {
			integralSum += ElementaryChebIntegral(masOfIntervals[i], masOfIntervals[i + 1], funVal, n, masOdNodes);
		}
		if (abs(factValue - integralSum) < d && !isReady) {
			printf("%lf\n", integralSum);
			printf("%lf\n", d / factValue);
			printf("%i\n", n);
			isReady = true;
		}
	}
	return integralSum;
}

int main(void) {
	double masOfNodes[3] = {-0.707107, 0, 0.707107};
    double a = 2;
	double b = 4;
	double e = 0.0001;
	double d = 0.1;
	double fact = GetRealIntegralValue(a, b);
	double res = NumericChebIntegration(a, b, e, testFun, masOfNodes, 3, fact, d);
	
	/*double a = 2;
	double b = 4;
	double res = ChebIntegral(a, b, testFun, 3, masOfNodes);*/
}