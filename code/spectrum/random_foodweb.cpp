// including .cpp-files
#include "/conv1/zfg663/foodwebs/code/evolution/distributions.cpp"
#include "/conv1/zfg663/foodwebs/code/evolution/species.cpp"
#include "/conv1/zfg663/foodwebs/code/evolution/stability_analysis.cpp"
#include "/conv1/zfg663/foodwebs/code/evolution/food_web.cpp"

// including templates
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

// defining print-function
void printMatrix(MatrixXd M, int N) {
	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) { cout << "C(" << i << "," << j << "): " << M(i, j) << " "; }
		cout << endl;
	}
}

bool feasibleWeb(Species S[], Producer P[]) {
	// declaring vectors and matrices
	MatrixXd R(Species::nTotal, Species::nTotal);
	VectorXd K(Species::nTotal);
	VectorXd Ssteady(Species::nTotal);
	
	// initializing R and K
	R = initializeR(R, S, P);				//cout << R << endl;
	K = initializeK(K, S, P);				//cout << K << endl;
	
	// computing steady states
	Ssteady = R.inverse() * K;
	
	if (feasible(Ssteady)) { return true; }
	else { return false; }
}

// initializing global variables from foodwebs
int Species::nTotal = 0;
int Producer::nProducer = 0;
bool FoodWeb::feasible = 0;
bool FoodWeb::stable = 0;

int main() {
// random seed
srand(27);

// parameters
int Nrep = 1e5;		// number of repetitions
int N;			// number of species
cout << "Enter size of community matrix: " << endl;
cin >> N;

// initializing normal distribution
int mu = 0;
int sig = 1;
random_device rd;
mt19937 e2(rd());
normal_distribution<> dist(mu, sig);

// opening file for saving data
ofstream file("/conv1/zfg663/foodwebs/data/spectra/FWstruct_N" + to_string(N) + ".txt");

for (int rep=0; rep<Nrep; rep++) {
	// arrays	
	MatrixXd C(N,N);
	VectorXcd E(N);
	Species S[N];
	Producer P[N];

//	for (int i = 0; i < N; i++) {
//		for (int j = 0; j < N; j++) {
//			cout << S[i].consumers[j] << " " << S[i].resources[j] << " ";
//		}
//		cout << endl;
//	}


	// assembling food web
	Producer s(0);
	P[0] = s;
	S[0] = s;
	
	for (int i=1; i<N; i++) {
		addSpecies(S, P, i);
		updateTrophicLevel(S);
	}
	
	// assigning random densities
	for (int i=0; i<Producer::nProducer; i++) {
		S[i].density = dist(e2);
		P[i].density = S[i].density;
	}
	for (int i=Producer::nProducer; i<N; i++) {
		S[i].density = dist(e2);
	}
	
	// checking if feasible
	if (feasibleWeb(S, P)) {
		// computing eigenvalues of Jacobian
		C = CommunityMatrix(C, S, P);
		EigenSolver<MatrixXd> es(C);
		E = es.eigenvalues();
	
		// saving
		for (int i=0; i<N; i++) {
			file << E(i).real() << " " << E(i).imag() << " ";
		}
		file << endl;
		//printMatrix(C, N);
	}

	// removing all species
	for (int i=N-1; i>=0; i--) {
		removeSpecies(S, P, i);
	}
}

file.close();
return 0;
}

