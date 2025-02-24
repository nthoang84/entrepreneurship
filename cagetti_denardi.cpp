#include "cagetti_denardi.h"

#include <fstream>
#include <iostream>

CagettiDeNardi::CagettiDeNardi() 
    : assetGridSize(300),
      incomeGridSize(5),         
      abilityGridSize(2),
      alpha(0.33),                      // Capital share of income
      beta(0.865),                      // Discount factor
      sigma(1.5),                       // Coefficient of relative risk aversion
      delta(0.06),                      // Depreciation rate
      eta(1.0),                         // Altruism
      nu(0.88),                         // Degree of decreasing returns to scale to entrepreneurial ability
      fracDefault(0.75),                // Fraction of working capital kept if default
      replacementRate(0.40),            // Pension and SS replacement rate
      fracCapitalConstraint(2.0)        // Max fraction of asset that can be borrowed for working capital
{
    totalGridSize = ageGridSize * typeGridSize * assetGridSize * incomeGridSize * abilityGridSize;
    totalGridSizeYoung = assetGridSize * incomeGridSize * abilityGridSize;
    totalGridSizeOld = assetGridSize * abilityGridSize;
    
    assetBounds = {0.0, 5000.0};
    interestRateBounds = make_pair(0.005, (1.0 / beta - 1.0));

    incomes = {.2468, .4473, .7654, 1.3097, 2.3742};
    transIncome = {
        {.7376, .2473, .0150, .0002, .0000},
        {.1947, .5555, .2328, .0169, .0001},
        {.0113, .2221, .5333, .2221, .0113},
        {.0001, .0169, .2328, .5555, .1947},
        {.0000, .0002, .0150, .2473, .7376}
    };

    abilities = {0.0, 0.514};
    transAbility = {
        {.964, .036},
        {.206, .794}
    };

    pi[young] = .978;
    pi[old] = .911;
}

void CagettiDeNardi::computeAssetGrid(double growthRate) {
    const double& assetMin = assetBounds.first;
    const double& assetMax = assetBounds.second;
    assets.resize(assetGridSize);
    if (fabs(growthRate) < EPS) {
        double stepSize = (assetMax - assetMin) / (assetGridSize - 1);
        for (int i = 0; i < assetGridSize; i++) {
            assets[i] = (i == 0 ? assetMin : assets[i - 1] + stepSize);
        }
        return;
    }
    for (int i = 0; i < assetGridSize; i++) {
        assets[i] = assetMin + (assetMax - assetMin) * 
                   ((pow(1 + growthRate, i) - 1) / (pow(1 + growthRate, assetGridSize - 1) - 1));
    }
}

void CagettiDeNardi::computeIncomeInvDist(double eps, int maxIter) {
    incomeInvDist.resize(incomeGridSize);
    incomeInvDist[0] = 1.0;
    int iter = 0;
    while (iter < maxIter) {
        vector<double> dist(incomeGridSize);
        for (int j = 0; j < incomeGridSize; j++) {
            for (int i = 0; i < incomeGridSize; i++) {
                dist[j] += incomeInvDist[i] * transIncome[i][j];
            }
        }
        double diff = 0;
        for (int i = 0; i < incomeGridSize; i++) {
            diff = max(diff, fabs(dist[i] - incomeInvDist[i]));
            incomeInvDist[i] = dist[i];
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
    double sum = 0;
    for (int i = 0; i < incomeGridSize; i++) {
        sum += incomeInvDist[i];
    }
    if (fabs(sum - 1.0) > EPS) {
        for (int i = 0; i < incomeGridSize; i++) {
            incomeInvDist[i] /= sum;
        }
    }
}

void CagettiDeNardi::computePolicy(double interestRate, double eps, int maxIter) {
    double wageRate = computeWageFromInterestRate(interestRate);
    assetPolicy.resize(totalGridSize);
    consumptionPolicy.resize(totalGridSize);
    v.resize(totalGridSize);
    vector<vector<double>> expectedV(ageGridSize);
    vector<vector<double>> dExpectedV(ageGridSize);
    vector<vector<double>> endogenousAssets(ageGridSize);
    expectedV[young].resize(totalGridSizeYoung);
    dExpectedV[young].resize(totalGridSizeYoung);
    endogenousAssets[young].resize(totalGridSizeYoung);

    auto updateGrid = [&]() {
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                for (int i = 0; i < assetGridSize; i++) {
                    double& discountedExpectedV = expectedV[young][id(i, j, t)];
                    discountedExpectedV = 0;
                    for (int jj = 0; jj < incomeGridSize; jj++) {
                        for (int tt = 0; tt < abilityGridSize; tt++) {
                            discountedExpectedV += beta * v[id(young, entrepreneur, i, jj, tt)] *
                                                   transIncome[j][jj] * transAbility[t][tt];
                        }
                    }
                    if (i == 1) {
                        dExpectedV[young][id(0, j, t)] = computeDerivative(
                            expectedV[young][id(0, j, t)], expectedV[young][id(1, j, t)],
                            assets[0], assets[1]
                        );
                    }
                    if (i == assetGridSize - 1) {
                        dExpectedV[young][id(i, j, t)] = computeDerivative(
                            expectedV[young][id(i - 1, j, t)], expectedV[young][id(i, j, t)],
                            assets[i - 1], assets[i]
                        );
                    }
                    if (i > 1) {
                        dExpectedV[young][id(i - 1, j, t)] = computeDerivative(
                            expectedV[young][id(i - 2, j, t)],
                            expectedV[young][id(i - 1, j, t)],
                            expectedV[young][id(i, j, t)],
                            assets[i - 2], assets[i - 1], assets[i]
                        );
                    }
                }
            }
        }
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                for (int i = 0; i < assetGridSize; i++) {
                    double maxWorkingCapital = fracCapitalConstraint * assets[i];
                    double modifiedInterestRate = interestRate;
                    if (fabs(maxWorkingCapital) > EPS) {
                        modifiedInterestRate = max(
                            modifiedInterestRate, 
                            interestRate + fracCapitalConstraint * (abilities[t] * pow(maxWorkingCapital, nu - 1.0) - interestRate - delta)
                        );
                    }
                    endogenousAssets[young][id(i, j, t)] = (mu_c_inverse(dExpectedV[young][id(i, j, t)]) +
                                                                assets[i] - wageRate * incomes[j]) /
                                                               (1 + modifiedInterestRate);
                }
            }
        }
    };

    // Initial guess where asset policy is zero everywhere
    for (int j = 0; j < incomeGridSize; j++) {
        for (int t = 0; t < abilityGridSize; t++) {
            for (int i = 0; i < assetGridSize; i++) {
                const int& index = id(young, entrepreneur, i, j, t);
                consumptionPolicy[index] = (1 + interestRate) * assets[i] + wageRate * incomes[j];
                v[index] = u(consumptionPolicy[index]) / (1 + beta);
            }
        }
    }
    updateGrid();

    int iter = 0;
    while (iter < maxIter) {
        auto vPrev = v;
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                int current_i = 1;
                for (int i = 0; i < assetGridSize; i++) {
                    while (current_i < assetGridSize - 1 && endogenousAssets[young][id(current_i, j, t)] < assets[i]) {
                        current_i++;
                    }
                    double weight = 0;
                    if (fabs(endogenousAssets[young][id(current_i, j, t)] - endogenousAssets[young][id(current_i - 1, j, t)]) > EPS) {
                        weight = (endogenousAssets[young][id(current_i, j, t)] - assets[i]) / 
                                    (endogenousAssets[young][id(current_i, j, t)] - endogenousAssets[young][id(current_i - 1, j, t)]);
                    }
                    weight = min(max(weight, 0.0), 1.0);
                    const int& index = id(young, entrepreneur, i, j, t);
                    double& nextAsset = assetPolicy[index];
                    nextAsset = weight * assets[current_i - 1] + (1 - weight) * assets[current_i];
                    nextAsset = max(nextAsset, assets[0]);
                    double tempV = weight * expectedV[young][id(current_i - 1, j, t)] + 
                                   (1 - weight) * expectedV[young][id(current_i, j, t)];
                    double maxWorkingCapital = fracCapitalConstraint * assets[i];
                    double modifiedInterestRate = interestRate;
                    if (fabs(maxWorkingCapital) > EPS) {
                        modifiedInterestRate = max(
                            modifiedInterestRate, 
                            interestRate + fracCapitalConstraint * (abilities[t] * pow(maxWorkingCapital, nu - 1.0) - interestRate - delta)
                        );
                    }
                    consumptionPolicy[index] = (1 + modifiedInterestRate) * assets[i] +
                                               wageRate * incomes[j] -
                                               assetPolicy[index];
                    v[index] = u(consumptionPolicy[index]) + tempV;
                }
            }
        }
        updateGrid();
        double diff = 0;
        for (int i = 0; i < totalGridSize; i++) {
            diff = max(diff, fabs(v[i] - vPrev[i]));
        }
        cout << "Iteration " << iter + 1 << ": diff = " << diff << '\n';
        if (diff < eps) {
            break;
        }
        iter++;
    }
}

void CagettiDeNardi::simulate(double eps, int maxIter) {
    vector<double> where(totalGridSize);
    vector<double> weight(totalGridSize);
    for (int j = 0; j < incomeGridSize; j++) {
        for (int t = 0; t < abilityGridSize; t++) {
            int current_i = 1;
            for (int i = 0; i < assetGridSize; i++) {
                const int& index = id(young, entrepreneur, i, j, t);
                double x = assetPolicy[index];
                while (current_i < assetGridSize - 1 && assets[current_i] < x) {
                    current_i++;
                }
                where[index] = current_i;
                weight[index] = (assets[current_i] - x) / (assets[current_i] - assets[current_i - 1]);
                weight[index] = min(max(weight[id(i, j)], 0.0), 1.0);
            }
        }
    }
    vector<double> dist(totalGridSize);
    dist[0] = 1.0;
    int iter = 0;
    while (iter < maxIter) {
        vector<double> tempDist(totalGridSize);
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                for (int i = 0; i < assetGridSize; i++) {
                    const int& index = id(young, entrepreneur, i, j, t);
                    double new_i = where[index];
                    double w = weight[index];
                    tempDist[id(young, entrepreneur, new_i - 1, j, t)] += w * dist[index];
                    tempDist[id(young, entrepreneur, new_i, j, t)] += (1 - w) * dist[index];
                }
            }            
        }
        vector<double> newDist(totalGridSize);
        for (int i = 0; i < assetGridSize; i++) {
            for (int jj = 0; jj < incomeGridSize; jj++) {
                for (int tt = 0; tt < abilityGridSize; tt++) {
                    for (int j = 0; j < incomeGridSize; j++) {
                        for (int t = 0; t < abilityGridSize; t++) {
                            newDist[id(young, entrepreneur, i, jj, tt)] += tempDist[id(young, entrepreneur, i, j, t)] * transIncome[j][jj] * transAbility[t][tt];
                        }
                    }
                }
            }
        }
        double diff = 0;
        for (int i = 0; i < totalGridSize; i++) {
            diff = max(diff, fabs(newDist[i] - dist[i]));
            dist[i] = newDist[i];
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
    double totalSavings = 0;
    vector<double> assetDist(assetGridSize);
    for (int i = 0; i < assetGridSize; i++) {
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                assetDist[i] += dist[id(young, entrepreneur, i, j, t)];
            }
        }
        totalSavings += assets[i] * assetDist[i];
    }
    string dataFile = "./data/assetDistribution.dat";
    ofstream dataStream(dataFile);
    for (int i = 0; i < assetGridSize; i++) {
        dataStream << assets[i] << " " << assetDist[i] * 100 << endl; 
    }
    dataStream.close();
    FILE *gnuplot = popen("gnuplot", "w");
    if (!gnuplot) {
        cerr << "Error: Unable to open gnuplot." << endl;
        return;
    }
    fprintf(gnuplot, "set terminal pdfcairo\n");
    fprintf(gnuplot, "set output './figures/assetDistribution.pdf'\n");
    fprintf(gnuplot, "set xlabel 'Asset'\n");
    fprintf(gnuplot, "set ylabel 'Percentage of agents'\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "plot '%s' with lines\n", dataFile.c_str());
    fflush(gnuplot);
    pclose(gnuplot);
}

inline int CagettiDeNardi::id(int age, int type, int asset, int income, int ability) {
    return age * (typeGridSize * assetGridSize * incomeGridSize * abilityGridSize) + 
           type * (assetGridSize * incomeGridSize * abilityGridSize) + 
           asset * (incomeGridSize * abilityGridSize) + 
           income * abilityGridSize +
           ability;
}

inline int CagettiDeNardi::id(int asset, int income, int ability) {
    return asset * (incomeGridSize * abilityGridSize) + 
           income * abilityGridSize +
           ability;
}

inline int CagettiDeNardi::id(int asset, int ability) {
    return asset * abilityGridSize + ability;
}

inline double CagettiDeNardi::u(double c) {
    return pow(c, 1 - sigma) / (1 - sigma);
}

inline double CagettiDeNardi::mu_c(double c) {
    return pow(c, -sigma);
}

inline double CagettiDeNardi::mu_c_inverse(double u) {
    double c = pow(u, -1 / sigma);
    return c;
}

inline double CagettiDeNardi::computeWageFromInterestRate(double interestRate) {
    return (1 - alpha) * pow(alpha / (interestRate + delta), alpha / (1 - alpha));
}

inline double CagettiDeNardi::computeDerivative(double y1, double y2, double x1, double x2) {
    double d = (y2 - y1) / (x2 - x1);
    return d;
}

inline double CagettiDeNardi::computeDerivative(double y1, double y2, double y3, double x1, double x2, double x3) {
    double d = ((1.0 - (x3 - x2)/(x3 - x1)) * ((y3 - y2) / (x3 - x2)) + 
                ((x3 - x2) / (x3 - x1)) * ((y2 - y1) / (x2 - x1)));
    return d;
}