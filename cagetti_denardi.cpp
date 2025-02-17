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
      replacementRate(0.40)             // Pension and SS replacement rate
{
    totalGridSize = ageGridSize * typeGridSize * assetGridSize * incomeGridSize * abilityGridSize;
    totalGridSizeYoung = assetGridSize * incomeGridSize * abilityGridSize;
    totalGridSizeOld = assetGridSize * abilityGridSize;
    
    assetBounds = {0.0, 300.0};
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
    // TODO: Maybe aggregate labor/income later
}

void CagettiDeNardi::computePolicy(double interestRate, double eps, int maxIter) {
    double wageRate = computeWageFromInterestRate(interestRate);
    assetPolicy.resize(totalGridSize);
    consumptionPolicy.resize(totalGridSize);
    v.resize(totalGridSize);
    V.resize(ageGridSize);
    V[young].resize(totalGridSizeYoung);
    V[old].resize(totalGridSizeOld);
    vector<vector<double>> expectedV(ageGridSize);
    expectedV[young].resize(totalGridSizeYoung);
    expectedV[old].resize(totalGridSizeOld);

    auto find = [&](double retainedCapital) {
        int where = -1, low = 0, high = assetGridSize - 1;
        while (low <= high) {
            int mid = (low + high) >> 1;
            if (assets[mid] <= retainedCapital) {
                where = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return where;
    };

    auto update = [&]() {
        for (int j = 0; j < incomeGridSize; j++) {
            for (int t = 0; t < abilityGridSize; t++) {
                for (int i = 0; i < assetGridSize; i++) {
                    // Update value functions
                    V[young][id(i, j, t)] = max(v[id(young, entrepreneur, i, j, t)], v[id(young, worker, i, j, t)]);
                    if (j == none) {
                        V[old][id(i, t)] = max(v[id(old, entrepreneur, i, none, t)], v[id(old, retiree, i, none, none)]);
                    }
                    // Update expected value functions
                    expectedV[young][id(i, j, t)] = 0;
                    for (int jj = 0; jj < incomeGridSize; jj++) {
                        for (int tt = 0; tt < abilityGridSize; tt++) {
                            expectedV[young][id(i, j, t)] += V[young][id(i, jj, tt)] * transIncome[j][jj] * transAbility[t][tt];
                        }
                    }
                    if (j == none) {
                        expectedV[old][id(i, t)] = 0;
                        for (int tt = 0; tt < abilityGridSize; tt++) {
                            expectedV[old][id(i, t)] += V[old][id(i, tt)] * transAbility[t][tt];
                        }
                    }
                }
            }
        }
    };

    auto vPrev = v;
    int iter = 0;
    while (iter < maxIter) {
        for (int j = 0; j < incomeGridSize; j++) {
            for (int i = 0; i < assetGridSize; i++) {
                for (int t = 0; t < abilityGridSize; t++) {
                    double currentAsset = assets[i];
                    double currentIncome = incomes[j];
                    double theta = abilities[t];
                    {
                        // Young entrepreneur
                        const int index = id(young, entrepreneur, i, j, t);
                        double& bestValue = v[index];
                        double& bestAsset = assetPolicy[index];
                        double& bestConsumption = consumptionPolicy[index];
                        for (int k = assetGridSize - 1; k >= 0; k--) {
                            double investment = assets[k];
                            if (investment < currentAsset) {
                                break;
                            }
                            bestValue = -INF;
                            for (int ii = 0; ii < assetGridSize; ii++) {
                                double nextAsset = assets[ii];
                                double consumption = (1 - delta) * investment + 
                                                     theta * pow(investment, nu) -
                                                     (1 + interestRate) * (investment - currentAsset) -
                                                     nextAsset;
                                double value = u(consumption) + 
                                               beta * pi[young] * expectedV[young][id(ii, j, t)] +
                                               beta * (1 - pi[young]) * expectedV[old][id(ii, t)];
                                if (value < v[id(young, worker, find(fracDefault * investment), j, t)]) {
                                    continue;
                                }
                                if (value > bestValue) {
                                    bestValue = value;
                                    bestAsset = nextAsset;
                                    bestConsumption = consumption;
                                }
                            }
                        }
                    }
                    {
                        // Young worker
                        const int index = id(young, worker, i, j, t);
                        double& bestValue = v[index];
                        double& bestAsset = assetPolicy[index];
                        double& bestConsumption = consumptionPolicy[index];
                        bestValue = -INF;
                        for (int ii = 0; ii < assetGridSize; ii++) {
                            double nextAsset = assets[ii];
                            double consumption = (1 + interestRate) * currentAsset + 
                                                 (1 - taxRate) * wageRate * currentIncome -
                                                 nextAsset;
                            double value = u(consumption) + 
                                           beta * pi[young] * expectedV[young][id(ii, j, t)] +
                                           beta * (1 - pi[young]) * v[id(old, retiree, nextAsset, none, none)];
                            if (value > bestValue) {
                                bestValue = value;
                                bestAsset = nextAsset;
                                bestConsumption = consumption;
                            }
                        }
                    }
                    if (j == none) {
                        // Old entrepreneur
                        const int index = id(old, entrepreneur, i, none, t);
                        double& bestValue = v[index];
                        double& bestAsset = assetPolicy[index];
                        double& bestConsumption = consumptionPolicy[index];
                        for (int k = assetGridSize - 1; k >= 0; k--) {
                            double investment = assets[k];
                            if (investment < currentAsset) {
                                break;
                            }
                            bestValue = -INF;
                            for (int ii = 0; ii < assetGridSize; ii++) {
                                double nextAsset = assets[ii];
                                double consumption = (1 - delta) * investment + 
                                                     theta * pow(investment, nu) -
                                                     (1 + interestRate) * (investment - currentAsset) -
                                                     nextAsset;
                                double value = u(consumption) + 
                                               beta * pi[old] * expectedV[old][id(ii, t)] +
                                               eta * beta * (1 - pi[old]) * expectedV[young][id(ii, j, t)];         // TODO: expectedV must be adjusted because nextIncome is unconditional
                                if (value < v[id(old, retiree, find(fracDefault * investment), none, none)]) {
                                    continue;
                                }
                                if (value > bestValue) {
                                    bestValue = value;
                                    bestAsset = nextAsset;
                                    bestConsumption = consumption;
                                }
                            }
                        }
                    }
                    if (j == none && t == none) {
                        // Old retiree
                        const int index = id(old, retiree, i, none, none);
                        double& bestValue = v[index];
                        double& bestAsset = assetPolicy[index];
                        double& bestConsumption = consumptionPolicy[index];
                        bestValue = -INF;
                        for (int ii = 0; ii < assetGridSize; ii++) {
                            double nextAsset = assets[ii];
                            double consumption = (1 + interestRate) * currentAsset - nextAsset;     // TODO: Social security payments can be added later
                            double value = u(consumption) + 
                                           beta * pi[old] * v[id(old, retiree, ii, none, none)] +
                                           eta * beta * (1 - pi[old]) * expectedV[young][id(ii, j, t)];     // TODO: expectedV must be adjusted because nextIncome is unconditional
                            if (value > bestValue) {
                                bestValue = value;
                                bestAsset = nextAsset;
                                bestConsumption = consumption;
                            }
                        }
                    }
                }
            }
        }
        double diff = 0;
        for (int i = 0; i < totalGridSize; i++) {
            diff = max(diff, fabs(v[i] - vPrev[i]));
        }
        if (diff < eps) {
            break;
        }
        update();
        vPrev = v;
        iter++;
        cout << "Iteration " << iter << ": diff = " << diff << '\n';
    }
}

void CagettiDeNardi::simulate(double eps, int maxIter) {
    vector<double> where(totalGridSize);
    vector<double> weight(totalGridSize);
    for (int j = 0; j < incomeGridSize; j++) {
        for (int t = 0; t < abilityGridSize; t++) {
            for (int age = 0; age < ageGridSize; age++) {
                for (int type = 0; type < typeGridSize; type++) {
                    if (age == old && j != none) {
                        continue;
                    }
                    if (age == old && type == retiree && t != none) {
                        continue;
                    }
                    int current_i = 0;
                    for (int i = 0; i < assetGridSize; i++) {
                        const int index = id(age, type, i, j, t);
                        double nextAsset = assetPolicy[index];
                        while (current_i < assetGridSize && assets[current_i] < nextAsset) {
                            current_i++;
                        }
                        current_i = min(current_i, assetGridSize - 1);
                        where[index] = current_i;
                        weight[index] = 0; // TODO: Interpolation will be added later
                    }                    
                }
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
                    {
                        int index = id(young, entrepreneur, i, j, t);
                        int next_i = where[index];
                        tempDist[id(young, entrepreneur, next_i, j, t)] += dist[index];
                        // TODO: I am unsure where to update the next distribution, added later
                    }
                }
            }
        }

        vector<double> newDist(totalGridSize);
        // TODO: Update with transition probabilities

        double diff = 0;
        for (int i = 0; i < totalGridSize; i++) {
            diff = max(diff, fabs(dist[i] - newDist[i]));
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
                assetDist[i] += dist[id(young, worker, i, j, t)];
                if (j == none) {
                    assetDist[i] += dist[id(old, entrepreneur, i, none, t)];
                }
                if (j == none && t == none) {
                    assetDist[i] += dist[id(old, retiree, i, none, none)];
                }
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
           asset * (incomeGridSize * ability) + 
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