#ifndef CAGETTI_DENARDI_H
#define CAGETTI_DENARDI_H

#include <vector>

using namespace std;

class CagettiDeNardi {
public:
    CagettiDeNardi();

    void computeAssetGrid(double growthRate = 0.025);

    void computeIncomeInvDist(double eps = EPS, int maxIter = MAX_ITER);

    void computePolicy(double interestRate, double eps = EPS, int maxIter = MAX_ITER);

    void computePolicyEuler(double interestRate, double eps = EPS, int maxIter = MAX_ITER);

    void simulate(bool plotDistribution = false, double eps = EPS, int maxIter = MAX_ITER);

    void solveEquilibrium(double eps = EPS, int maxIter = MAX_ITER);

    void plot();

    double computeGiniCoefficient(const vector<double>& wealth, const vector<double>& dist);

    double computeParetoCoefficient(const vector<double>& wealth, const vector<double>& dist, double threshold = 1.0);

    void debug();

private:
    inline int id(int age, int type, int asset, int income, int ability);
    inline int id(int asset, int income, int ability);
    inline int id(int asset, int ability);
    inline double u(double c);
    inline double mu_c(double c);
    inline double mu_c_inverse(double u);
    inline double computeWageFromInterestRate(double interestRate);
    inline double computeModifiedInterestRate(double interestRate, double asset, double ability);
    inline double computeDerivative(double y1, double y2, double x1, double x2);
    inline double computeDerivative(double y1, double y2, double y3, double x1, double x2, double x3);
    inline double computePercentile(const vector<double>& wealth, const vector<double>& dist, double k);

    static constexpr double INF = 1e9;
    static constexpr double EPS = 1e-6;
    static constexpr int MAX_ITER = (int) 1e3;
    static constexpr int none = 0;
    static constexpr int ageGridSize = 1;       // Only young agents in this version
    static constexpr int typeGridSize = 1;      // Trick: Entrepreneurs are also waged workers!
    
    int assetGridSize;
    int incomeGridSize;
    int abilityGridSize;
    int totalGridSize;
    int totalGridSizeYoung;
    int totalGridSizeOld;
    double alpha;
    double beta;
    double sigma;
    double delta;
    double eta;
    double nu;
    double fracDefault;
    double fracCapitalConstraint;
    double replacementRate;
    double wageRate;
    double taxRate;
    double aggregateIncome;
    double aggregateAssetSupply;
    double aggregateAssetDemand;
    double eqmInterestRate;
    array<double, 2> pi;
    pair<double, double> assetBounds;
    pair<double, double> interestRateBounds;
    
    vector<double> assets;
    vector<double> incomes;
    vector<double> abilities;
    vector<vector<double>> transIncome;
    vector<vector<double>> transAbility;
    vector<double> incomeInvDist;
    vector<double> assetPolicy;
    vector<double> consumptionPolicy;
    vector<double> v;
    vector<double> expectedV;
    
    enum Age {
        young = 0,
        old = 1
    };

    enum Type {
        entrepreneur = 0,
        worker = 1,
        retiree = 1
    };
};

#endif