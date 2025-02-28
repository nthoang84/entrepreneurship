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
      fracCapitalConstraint(1.5)        // Max fraction of asset that can be borrowed for working capital
{
    totalGridSize = ageGridSize * typeGridSize * assetGridSize * incomeGridSize * abilityGridSize;
    totalGridSizeYoung = assetGridSize * incomeGridSize * abilityGridSize;
    totalGridSizeOld = assetGridSize * abilityGridSize;
    
    assetBounds = {0.0, 3000.0};
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
    aggregateIncome = 0;
    for (int i = 0; i < incomeGridSize; i++) {
        aggregateIncome += incomeInvDist[i] * incomes[i];
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
                    double modifiedInterestRate = computeModifiedInterestRate(interestRate, assets[i], abilities[t]);
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
                    if (assets[i] < endogenousAssets[young][id(young, entrepreneur, 0, j, t)]) {
                        const int& index = id(young, entrepreneur, i, j, t);
                        double& nextAsset = assetPolicy[index];
                        nextAsset = assets[0];
                        double tempV = expectedV[young][id(young, entrepreneur, 0, j, t)];
                        double modifiedInterestRate = computeModifiedInterestRate(interestRate, assets[i], abilities[t]);
                        consumptionPolicy[index] = (1 + modifiedInterestRate) * assets[i] +
                                                   wageRate * incomes[j] -
                                                   nextAsset;
                        v[index] = u(consumptionPolicy[index]) + tempV;
                    } else {
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
                        double modifiedInterestRate = computeModifiedInterestRate(interestRate, assets[i], abilities[t]);
                        consumptionPolicy[index] = (1 + modifiedInterestRate) * assets[i] +
                                                   wageRate * incomes[j] -
                                                   nextAsset;
                        v[index] = u(consumptionPolicy[index]) + tempV;
                    }
                }
            }
        }
        updateGrid();
        double diff = 0;
        for (int i = 0; i < totalGridSize; i++) {
            diff = max(diff, fabs(v[i] - vPrev[i]));
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
}

void CagettiDeNardi::computePolicyEuler(double interestRate, double eps, int maxIter) {
    double wageRate = computeWageFromInterestRate(interestRate);
    assetPolicy.resize(totalGridSize);
    consumptionPolicy.resize(totalGridSize);
    vector<double> MU(totalGridSize);
    vector<double> expectedMU(totalGridSize);
    vector<double> endogenousAssets(totalGridSize);

    auto updateGrid = [&]() {
        for (int t = 0; t < abilityGridSize; t++) {
            for (int j = 0; j < incomeGridSize; j++) {
                for (int i = 0; i < assetGridSize; i++) {
                    const int& index = id(young, entrepreneur, i, j, t);
                    double modifiedInterestRate = computeModifiedInterestRate(interestRate, assets[i], abilities[t]);
                    consumptionPolicy[index] = (1 + modifiedInterestRate) * assets[i] + 
                                               wageRate * incomes[j] -
                                               assetPolicy[index];
                    MU[index] = mu_c(consumptionPolicy[index]);
                }
            }
        }
        for (int t = 0; t < abilityGridSize; t++) {
            for (int j = 0; j < incomeGridSize; j++) {
                for (int i = 0; i < assetGridSize; i++) {
                    const int& index = id(young, entrepreneur, i, j, t);
                    expectedMU[index] = 0;
                    double modifiedInterestRate = computeModifiedInterestRate(interestRate, assets[i], abilities[t]);
                    for (int tt = 0; tt < abilityGridSize; tt++) {
                        for (int jj = 0; jj < incomeGridSize; jj++) {
                            expectedMU[index] += beta * (1 + modifiedInterestRate) *
                                                 transIncome[j][jj] * transAbility[t][tt] * 
                                                 MU[id(young, entrepreneur, i, jj, tt)];
                        }
                    }
                    endogenousAssets[index] = (mu_c_inverse(expectedMU[index]) + 
                                               assets[i] - wageRate * incomes[j]) / 
                                              (1 + modifiedInterestRate);
                }
            }
        }
    };
    
    // Start with an initial guess where tomorrow's assets are zero everywhere
    updateGrid();

    int iter = 0;
    while (iter < maxIter) {     
        for (int t = 0; t < abilityGridSize; t++) {
            for (int j = 0; j < incomeGridSize; j++) {
                int current_i = 1;
                for (int i = 0; i < assetGridSize; i++) {
                    while (current_i < assetGridSize - 1 && endogenousAssets[id(young, entrepreneur, current_i, j, t)] < assets[i]) {
                        current_i++;
                    }
                    double weight = (endogenousAssets[id(young, entrepreneur, current_i, j, t)] - assets[i]) / 
                                    (endogenousAssets[id(young, entrepreneur, current_i, j, t)] - endogenousAssets[id(young, entrepreneur, current_i - 1, j, t)]);
                    weight = min(max(weight, 0.0), 1.0);
                    const int& index = id(young, entrepreneur, i, j, t);
                    assetPolicy[index] = weight * assets[current_i - 1] + (1 - weight) * assets[current_i];
                    assetPolicy[index] = max(assetPolicy[index], assets[0]);
                }
            }
        }  
        auto prevMU = MU;
        updateGrid();
        double diff = 0;
        for (int i = 0; i < totalGridSize; i++) {
            diff = max(diff, fabs(MU[i] - prevMU[i]));
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
}

void CagettiDeNardi::simulate(bool plotDistribution, double eps, int maxIter) {
    vector<double> where(totalGridSize);
    vector<double> weight(totalGridSize);
    for (int j = 0; j < incomeGridSize; j++) {
        for (int t = 0; t < abilityGridSize; t++) {
            int current_i = 1;
            for (int i = 0; i < assetGridSize; i++) {
                const int& index = id(young, entrepreneur, i, j, t);
                while (current_i < assetGridSize - 1 && assets[current_i] < assetPolicy[index]) {
                    current_i++;
                }
                where[index] = current_i;
                weight[index] = (assets[current_i] - assetPolicy[index]) / (assets[current_i] - assets[current_i - 1]);
                weight[index] = min(max(weight[index], 0.0), 1.0);
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
        double sum = 0;
        for (int i = 0; i < assetGridSize; i++) {
            for (int j = 0; j < incomeGridSize; j++) {
                for (int t = 0; t < abilityGridSize; t++) {
                    for (int jj = 0; jj < incomeGridSize; jj++) {
                        for (int tt = 0; tt < abilityGridSize; tt++) {
                            if(tempDist[id(young, entrepreneur, i, j, t)] > 0){
                                newDist[id(young, entrepreneur, i, jj, tt)] += 
                                    tempDist[id(young, entrepreneur, i, j, t)] * 
                                    transIncome[j][jj] * transAbility[t][tt];
                            }
                        }
                    }
                }
            }
        }
        double diff = 0;
        for (int i = 0; i < assetGridSize; i++) {
            for (int j = 0; j < incomeGridSize; j++) {
                for (int t = 0; t < abilityGridSize; t++) {
                    const int& index = id(young, entrepreneur, i, j, t);
                    diff = max(diff, fabs(newDist[index] - dist[index]));
                    dist[index] = newDist[index];
                }
            }
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
    aggregateAssetSupply = 0;
    vector<double> assetDist(assetGridSize);
    for (int i = 0; i < assetGridSize; i++) {
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                assetDist[i] += dist[id(young, entrepreneur, i, j, t)];
            }
        }
    }
    for (int i = 0; i < assetGridSize; i++) {
        aggregateAssetSupply += assets[i] * assetDist[i];
    }
    if (plotDistribution) {
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
}

void CagettiDeNardi::solveEquilibrium(double eps, int maxIter) {
    double interestRate = 0.04;
    int iter = 0;
    while (iter < maxIter) {
        double wageRate = computeWageFromInterestRate(interestRate);
        aggregateAssetDemand = aggregateIncome * pow(alpha / (interestRate + delta), 1.0 / (1.0 - alpha));
        computePolicy(interestRate);
        simulate();
        double diff = (aggregateAssetSupply - aggregateAssetDemand) /
                      ((aggregateAssetSupply + aggregateAssetDemand) / 2);
        cout << "Iteration " << iter + 1 << ": r = " << interestRate << ", diff = " << diff << '\n';
        if (fabs(diff) < eps) {
            break;
        }
        interestRate = alpha * pow(aggregateIncome, 1 - alpha) / pow((aggregateAssetSupply + aggregateAssetDemand) / 2, 1 - alpha) - delta;
        interestRate = max(interestRate, interestRateBounds.first);
        iter++;
    }
    eqmInterestRate = interestRate;
    simulate(true);
    cout << ">> Equilibrium: r = " << eqmInterestRate << ", K(r) = " << aggregateAssetSupply << '\n';
}

void CagettiDeNardi::plot() {
    double interestRateStep = 0.005;
    vector<double> interestRate;
    for (double currentRate = interestRateBounds.first; currentRate < interestRateBounds.second; ) {
        interestRate.push_back(currentRate);
        currentRate += interestRateStep;
        if (fabs(currentRate - interestRateBounds.second) < EPS) {
            interestRate.push_back(currentRate);
        }
    }
    vector<double> demand(interestRate.size());
    vector<double> supply(interestRate.size());
    for (int i = 0; i < interestRate.size(); i++) {
        double wageRate = computeWageFromInterestRate(interestRate[i]);
        demand[i] = aggregateIncome * pow(alpha / (interestRate[i] + delta), 1.0 / (1.0 - alpha));
        computePolicy(interestRate[i]);
        simulate();
        supply[i] = aggregateAssetSupply;
    }
    FILE *gnuplot = popen("gnuplot", "w");
    if (!gnuplot) {
        cerr << "Error: Unable to open gnuplot." << endl;
        return;
    }
    fprintf(gnuplot, "set terminal pdfcairo\n");
    fprintf(gnuplot, "set output './figures/assetSupplyDemand.pdf'\n");
    fprintf(gnuplot, "set xlabel 'Aggregate wealth'\n");
    fprintf(gnuplot, "set ylabel 'Interest rate'\n");
    fprintf(gnuplot, "plot '-' with lines title 'Demand', '-' with lines title 'Supply'\n");
    for (size_t i = 0; i < interestRate.size(); ++i) {
        fprintf(gnuplot, "%lf %lf\n", demand[i], interestRate[i]);
    }
    fprintf(gnuplot, "e\n");
    for (size_t i = 0; i < interestRate.size(); ++i) {
        fprintf(gnuplot, "%lf %lf\n", supply[i], interestRate[i]);
    }
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
    pclose(gnuplot);
}

void CagettiDeNardi::debug() {
    {
        string dataFile = "./data/assetPolicy.dat";
        ofstream dataStream(dataFile);    
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                dataStream << ">> INCOME = " << incomes[j] << "\tABILITY = " << abilities[t] << endl;
                for (int i = 0; i < assetGridSize; i++) {
                    dataStream << fixed << setprecision(6)
                        << setw(15) << left << assets[i]
                        << setw(15) << left << assetPolicy[id(young, entrepreneur, i, j, t)]
                        << endl;
                }
                dataStream << endl;
            }
        }
        dataStream.close();
    }
    {
        double interestRate = 0.04;
        string dataFile = "./data/returns.dat";
        ofstream dataStream(dataFile);    
        for (int i = 0; i < assetGridSize; i++) {
            for (int t = 0; t < abilityGridSize; t++) {
                double modifiedInterestRate = computeModifiedInterestRate(interestRate, assets[i], abilities[t]);
                dataStream << fixed << setprecision(6)
                    << setw(15) << left << assets[i]
                    << setw(15) << left << abilities[t]
                    << setw(15) << left << modifiedInterestRate
                    << endl;
            }
        }
        dataStream.close();
    }
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
    return pow(u, -1 / sigma);
}

inline double CagettiDeNardi::computeWageFromInterestRate(double interestRate) {
    return (1 - alpha) * pow(alpha / (interestRate + delta), alpha / (1 - alpha));
}

inline double CagettiDeNardi::computeModifiedInterestRate(double interestRate, double asset, double ability) {
    double maxWorkingCapital = fracCapitalConstraint * asset;
    double modifiedInterestRate = interestRate;
    if (fabs(maxWorkingCapital) > EPS) {
        modifiedInterestRate = max(
            modifiedInterestRate, 
            interestRate + fracCapitalConstraint * (ability * pow(maxWorkingCapital, nu - 1.0) - interestRate - delta)
        );
    }
    return modifiedInterestRate;
}

inline double CagettiDeNardi::computeDerivative(double y1, double y2, double x1, double x2) {
    return (y2 - y1) / (x2 - x1);
}

inline double CagettiDeNardi::computeDerivative(double y1, double y2, double y3, double x1, double x2, double x3) {
    return (1.0 - (x3 - x2) / (x3 - x1)) * ((y3 - y2) / (x3 - x2)) + 
           ((x3 - x2) / (x3 - x1)) * ((y2 - y1) / (x2 - x1));
}