#include "cagetti_denardi.h"

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
    totalGridSize = 2 * 2 * assetGridSize * incomeGridSize * abilityGridSize;
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

void CagettiDeNardi::computePolicy(double interestRate, double eps, int maxIter) {
    double wageRate = computeWageFromInterestRate(interestRate);
    assetPolicy.resize(totalGridSize);
    consumptionPolicy.resize(totalGridSize);
    v.resize(totalGridSize);
    V.resize(2);
    V[young].resize(totalGridSizeYoung);
    V[old].resize(totalGridSizeOld);
    vector<vector<double>> expectedV(2);
    expectedV[young].resize(totalGridSizeYoung);
    expectedV[old].resize(totalGridSizeOld);

    auto find = [&](int i) {    // TODO: This needs to be filled in later with precise search procedure
        return i;
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
                        double& bestValue = v[id(young, entrepreneur, i, j, t)];
                        double& bestAsset = assetPolicy[id(young, entrepreneur, i, j, t)];
                        double& bestConsumption = consumptionPolicy[id(young, entrepreneur, i, j, t)];
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
                        double& bestValue = v[id(young, worker, i, j, t)];
                        double& bestAsset = assetPolicy[id(young, worker, i, j, t)];
                        double& bestConsumption = consumptionPolicy[id(young, worker, i, j, t)];
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
                        double& bestValue = v[id(old, entrepreneur, i, none, t)];
                        double& bestAsset = assetPolicy[id(old, entrepreneur, i, none, t)];
                        double& bestConsumption = consumptionPolicy[id(old, entrepreneur, i, none, t)];
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
                        double& bestValue = v[id(old, retiree, i, none, none)];
                        double& bestAsset = assetPolicy[id(old, retiree, i, none, none)];
                        double& bestConsumption = consumptionPolicy[id(young, worker, i, none, none)];
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
                    V[young][id(i, j, t)] = max(v[id(young, entrepreneur, i, j, t)], v[id(young, worker, i, j, t)]);
                    V[old][id(i, t)] = max(v[id(old, entrepreneur, i, none, t)], v[id(old, retiree, i, none, none)]);
                    // TODO: Update expectedV[.]
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
        vPrev = v;
        iter++;
    }
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