//
// Created by andre on 7/9/22.
//

#ifndef RESTRICTEDGENETICALGORITHM__UTILS_H_
#define RESTRICTEDGENETICALGORITHM__UTILS_H_

#include <random>
#include <vector>
#include <bitset>
#include <cmath>
#include <sstream>
#include <list>
#include <map>
#include <algorithm>
#include "Individual.h"
#include "lib/RandomGA.h"

/// Calculate exterior penalization without parameter R for problem problemIndex
/// \param x
/// \param problemIndex
/// \return
double restrictionsValue(const std::vector<double> &x, const size_t &problemIndex){
  // Convertimos todas las restricciones del tipo h a g
  double result = 0.0;
  double eps = 1e-4;
  if(problemIndex == 1){
    double g1Squared;
    double mult = 1.0;
    for(auto &x_i: x){
      mult *= x_i;
    }
    g1Squared = 0.75 - mult;
    g1Squared = 0.0 < g1Squared ? g1Squared * g1Squared : 0.0;
    double g2Squared;
    double sum = 0.0;
    for(auto &x_i: x){
      sum += x_i;
    }
    g2Squared = sum - (7.5*20.0);
    g2Squared = 0.0 < g2Squared ? g2Squared * g2Squared : 0.0;
    result = g1Squared + g2Squared;
  } else if (problemIndex == 2){
    double g1Squared = -x[3] + x[2] - 0.55;
    g1Squared = 0.0 < g1Squared ? g1Squared * g1Squared : 0.0;
    double g2Squared = -x[3] + x[2] - 0.55;
    g2Squared = 0.0 < g2Squared ? g2Squared * g2Squared : 0.0;
    double g3Squared = 1000.0*sin(-x[2]-0.25) + 1000.0*sin(-x[3]-0.25) + 894.8 - x[0];
    g3Squared = std::abs(g3Squared) - eps;
    g3Squared = 0.0 < g3Squared ? g3Squared * g3Squared : 0.0;
    double g4Squared = 1000.0*sin(x[2]-0.25) + 1000.0*sin(x[2]-x[3]-0.25) + 894.8 - x[1];
    g4Squared = std::abs(g4Squared) - eps;
    g4Squared = 0.0 < g4Squared ? g4Squared * g4Squared : 0.0;
    double g5Squared = 1000.0*sin(x[3]-0.25) + 1000.0*sin(x[3]-x[2]-0.25) + 1294.8;
    g5Squared = std::abs(g5Squared) - eps;
    g5Squared = 0.0 < g5Squared ? g5Squared * g5Squared : 0.0;
    result = g1Squared + g2Squared + g3Squared + g4Squared + g5Squared;
  } else if (problemIndex == 3){
    double g1Squared = 0.0;
    for(auto &x_i: x){
      g1Squared += x_i*x_i;
    }
    g1Squared -= 10;
    g1Squared = std::abs(g1Squared) - eps;
    g1Squared = 0.0 < g1Squared ? g1Squared * g1Squared : 0.0;
    double g2Squared = x[1]*x[2] - 5.0*x[3]*x[4];
    g2Squared = std::abs(g2Squared) - eps;
    g2Squared = 0.0 < g2Squared ? g2Squared * g2Squared : 0.0;
    double g3Squared = x[0]*x[0]*x[0] + x[1]*x[1] + 1.0;
    g3Squared = std::abs(g3Squared) - eps;
    g3Squared = 0.0 < g3Squared ? g3Squared * g3Squared : 0.0;
    result = g1Squared + g2Squared + g3Squared;
  } else {
    throw std::runtime_error("Número de problema no existe");
  }

  return result;
}

double restrictionsValueNormalized(const std::vector<double> &x, const size_t &problemIndex){

  // Mayor y menor valor de f(x) y maximo de valores de restricciones
  static double maxRestrictionValues = -std::numeric_limits<double>::max();
  static double minRestrictionValues = std::numeric_limits<double>::max();

  maxRestrictionValues = maxRestrictionValues < restrictionsValue(x, problemIndex) ? restrictionsValue(x, problemIndex) : maxRestrictionValues;
  minRestrictionValues = restrictionsValue(x, problemIndex) < minRestrictionValues ? restrictionsValue(x, problemIndex) : minRestrictionValues;
  double restrictionNormalized;

  if(std::abs(maxRestrictionValues - minRestrictionValues) < 1e-8){
    restrictionNormalized = restrictionsValue(x, problemIndex) / maxRestrictionValues;
  } else {
    restrictionNormalized = (restrictionsValue(x, problemIndex) - minRestrictionValues)/(maxRestrictionValues - minRestrictionValues);
  }

  return restrictionNormalized;
}

bool isFeasible(const std::vector<double> &x, const size_t &problemIndex){
  // Convertimos todas las restricciones del tipo h a g
  bool result = true;
  double eps = 1e-4;
  if(problemIndex == 1){
    double g1;
    double mult = 1.0;
    for(auto &x_i: x){
      mult *= x_i;
    }
    g1 = 0.75 - mult;
    if(g1 > 0.0){
      result = false;
    }

    double g2;
    double sum = 0.0;
    for(auto &x_i: x){
      sum += x_i;
    }
    g2 = sum - (7.5*20.0);
    if(g1 > 0.0){
      result = false;
    }

  } else if (problemIndex == 2){
    double g1 = -x[3] + x[2] - 0.55;
    if(g1 > 0.0){
      result = false;
    }

    double g2 = -x[3] + x[2] - 0.55;
    if(g2 > 0.0){
      result = false;
    }

    double g3 = 1000.0*sin(-x[2]-0.25) + 1000.0*sin(-x[3]-0.25) + 894.8 - x[0];
    if(std::abs(g3) - eps > 0.0){
      result = false;
    }

    double g4 = 1000.0*sin(x[2]-0.25) + 1000.0*sin(x[2]-x[3]-0.25) + 894.8 - x[1];
    if(std::abs(g4) - eps > 0.0){
      result = false;
    }

    double g5 = 1000.0*sin(x[3]-0.25) + 1000.0*sin(x[3]-x[2]-0.25) + 1294.8;
    if(std::abs(g5) - eps > 0.0){
      result = false;
    }

  } else if (problemIndex == 3){
    double g1 = 0.0;
    for(auto &x_i: x){
      g1 += x_i*x_i;
    }
    g1 -= 10;
    if(std::abs(g1) - eps > 0.0){
      result = false;
    }

    double g2 = x[1]*x[2] - 5.0*x[3]*x[4];
    if(std::abs(g2) - eps > 0.0){
      result = false;
    }

    double g3 = x[0]*x[0]*x[0] + x[1]*x[1] + 1.0;
    if(std::abs(g3) - eps > 0.0){
      result = false;
    }
  } else {
    throw std::runtime_error("Número de problema no existe");
  }

  return result;
}

/// Exactly what it is
/// \param x
/// \return problem 1 function value
double function1(const std::vector<double> &x){
  double sumCos4th = 0.0;
  double prodCos2th = 1.0;
  double sumSquared = 0.0;
  for (auto &x_i: x) {
    sumCos4th += cos(x_i)*cos(x_i)*cos(x_i)*cos(x_i);
  }

  for (auto &x_i: x) {
    prodCos2th *= cos(x_i)*cos(x_i);
  }

  double i = 1.0;
  for (auto &x_i: x) {
    sumSquared += i*x_i*x_i;
     i += 1.0;
  }

  double result = std::abs((sumCos4th - 2.0*prodCos2th)/(sqrt(sumSquared)));

  return result;
}

/// Exactly what it is
/// \param x
/// \return problem 2 function value
double function2(const std::vector<double> &x){
  return (3.0*x[0]) + (0.000001*x[0]*x[0]*x[0]) + (2.0*x[1]) + ((0.000002/3.0)*x[1]*x[1]*x[1]);
}

/// Exactly what it is
/// \param x
/// \return problem 3 function value
double function3(const std::vector<double> &x){
  return exp(x[0]*x[1]*x[2]*x[3]*x[4]);
}

/// Calculate function value for problemIndex
/// \param x
/// \param problemIndex
/// \return
double function(const std::vector<double> &x, const size_t &problemIndex) {
  double result;
  switch (problemIndex) {
    case 1:
      result = function1(x);
      break;
    case 2:
      result = function2(x);
      break;
    case 3:
      result = function3(x);
      break;
    default:
      throw std::runtime_error("No existe el problema seleccionado");
  }

  return result;
}

double functionNormalized(const std::vector<double> &x, const size_t &problemIndex){
  static double maxFx = -std::numeric_limits<double>::max();
  static double minFx = std::numeric_limits<double>::max();

  double functionNormalized;

  maxFx = maxFx < function(x, problemIndex) ? function(x, problemIndex) : maxFx;
  minFx = function(x, problemIndex) < minFx ? function(x, problemIndex) : minFx;
  if(std::abs(maxFx - minFx) < 1e-8){
    functionNormalized = function(x, problemIndex) / maxFx;
  } else {
    functionNormalized = (function(x, problemIndex) - minFx)/(maxFx - minFx);
  }

  return functionNormalized;
}

/// Calculate function with exterior penalised function
/// \param x
/// \param sigmoidParameters
/// \param problemIndex
/// \param generationIndex
/// \return
double functionPenalizedNormalized(const std::vector<double> &x, const std::vector<double> &sigmoidParameters, const size_t &problemIndex, const size_t &generationIndex){
  double result;
  double sigmoidValue = (sigmoidParameters[2]/(1 + exp(sigmoidParameters[0]*(-(round((double )generationIndex) - sigmoidParameters[1])))));

  switch (problemIndex) {
    case 1:
      result = functionNormalized(x, problemIndex) - sigmoidValue * restrictionsValueNormalized(x, problemIndex);
      break;
    case 2:
      result = functionNormalized(x, problemIndex) + sigmoidValue * restrictionsValueNormalized(x, problemIndex);
      break;
    case 3:
      result = functionNormalized(x, problemIndex) + sigmoidValue * restrictionsValueNormalized(x, problemIndex);
      break;
    default:
      throw std::runtime_error("No existe el problema seleccionado");
  }

  return result;
}

/// Calcular el promedio de un vector de pares usando el segundo
/// \param x
/// \return avg(x)
double avgFitnessPopulation(const std::vector<individual> &population) {
  double sum = 0.0;
  for (auto &individuo: population) {
    sum += (double) individuo.fitness;
  }
  return sum / (double) population.size();
}

/// Inicializa el arreglo de población
/// \param population: arreglo de población, asume que tiene el tamaño de la población
/// \param dim: dimension of chromosome string
void initalizePopulation(std::vector<individual> &population, const std::vector<std::pair<double, double>> &constraint) {

  // Inicializar cromosomas
  for (auto &individuo: population) {
    individuo.x = std::vector<double>(constraint.size(),0.0);
    for(size_t indexVariable = 0; indexVariable < individuo.x.size(); indexVariable++){
      // Las variables deben estar dentro del rango
      individuo.x[indexVariable] = rndreal(constraint[indexVariable].first, constraint[indexVariable].second);
    }
  }
}

/// Calcula la aptitud de cada individuo
/// \param population: arreglo de individuos
void fitnessCalculation(std::vector<individual> &population, const std::vector<double> &sigmoidParameters, const size_t &problemIndex, const size_t &generationIndex) {
  static double minFxPenalized = std::numeric_limits<double>::max();
  static double maxFxPenalized = -std::numeric_limits<double>::max();

  // Cálculo de Fx con penalidades
  for (auto &individuo: population) {
    individuo.Fx = function(individuo.x, problemIndex);
    individuo.fitness = functionPenalizedNormalized(individuo.x, sigmoidParameters, problemIndex, generationIndex);
    individuo.penaltiesValueNormalized = restrictionsValueNormalized(individuo.x, problemIndex);
    // Si las penalidades son cercanas a cero en epsilon, entonces está en la zona factible
    if(isFeasible(individuo.x, problemIndex)) {
      individuo.isFeasible = true;
    } else {
      individuo.isFeasible = false;
    }
  }

  double minFx = std::min_element(population.begin(), population.end(),
                                  [](const individual& i, const individual& j) {
                                    return i.fitness < j.fitness;
                                  })->fitness;

  double maxFx = std::max_element(population.begin(), population.end(),
                                  [](const individual& i, const individual& j) {
                                    return i.fitness < j.fitness;
                                  })->fitness;

  maxFxPenalized = maxFxPenalized < maxFx ? maxFx : maxFxPenalized;
  minFxPenalized = minFxPenalized > minFx ? minFx : minFxPenalized;

  if(std::abs(maxFxPenalized - minFxPenalized) < 1e-8){
    for (auto &individuo: population) {
      individuo.fitness /= maxFxPenalized;
    }
  } else {
    for (auto &individuo: population) {
      individuo.fitness = (individuo.fitness - minFxPenalized) / (maxFxPenalized - minFxPenalized);
    }
  }

  if(problemIndex == 2 || problemIndex == 3){
    // Ajuste de aptitud para minimizar
    for (auto &individuo: population) {
      individuo.fitness = 1.0 - individuo.fitness;
    }
  }
}

/// Selección por torneo binario inspirado en deb
/// \param population: arreglo de individuos
/// \return arreglo de pares de padre
std::vector<std::pair<size_t, size_t>> selection(const std::vector<individual> &population) {
  // Poner índices en un arreglo
  std::vector<size_t> populationIndexes(population.size());
  for (size_t i = 0; i < population.size(); ++i) {
    populationIndexes[i] = i;
  }

  // Inicializar dispositivo aleatorio con la misma semilla
  unsigned seed = 123654789;
  std::mt19937 g(seed);

  // Selección por torneo binario determinístico
  std::vector<std::pair<size_t, size_t>> selected;
  std::vector<size_t> preselected;

  for (size_t selectionRound = 0; selectionRound < 2; ++selectionRound) {
    // Barajear
    std::shuffle(populationIndexes.begin(), populationIndexes.end(), g);
    for (size_t popIndex = 0; popIndex < population.size(); popIndex += 2) {
      // Si uno es factible y otro no, entonces gana el factible
      if(population[populationIndexes[popIndex]].isFeasible != population[populationIndexes[popIndex + 1]].isFeasible){
        if(population[populationIndexes[popIndex]].isFeasible){
          preselected.push_back(populationIndexes[popIndex]);
        } else {
          preselected.push_back(populationIndexes[popIndex + 1]);
        }
        // Si ambos son factibles, entonces gana el que tiene mejor aptitud
      } else if (population[populationIndexes[popIndex]].isFeasible  && population[populationIndexes[popIndex + 1]].isFeasible){
        if (population[populationIndexes[popIndex + 1]].fitness < population[populationIndexes[popIndex]].fitness)
          preselected.push_back(populationIndexes[popIndex]);
        else
          preselected.push_back(populationIndexes[popIndex + 1]);
        // Si ninguno es factible, entonces gana quien tenga menor valor de violaciones
      } else {
        if (population[populationIndexes[popIndex + 1]].penaltiesValueNormalized < population[populationIndexes[popIndex]].penaltiesValueNormalized)
          preselected.push_back(populationIndexes[popIndex + 1]);
        else
          preselected.push_back(populationIndexes[popIndex]);
      }
    }
  }

  // Los padres son pares ordenados como: i-ésimo seleccionado de la ronda uno está emparejado con el i-ésimo de la segunda ronda.
  /*
   * (preselected[0],...,preselected[(n/2)-1]) son los seleccionados de la primera ronda.
   * (preselected[n/2],...,preselected[n]) son los seleccionados de la primera ronda.
   * n = tamaño de la población
   */
  for (size_t selectedIndex = 0; selectedIndex < (population.size()/2); ++selectedIndex) {
    selected.emplace_back(preselected[selectedIndex], preselected[(population.size() / 2) + selectedIndex]);
  }

  return selected;
}

///// Cruza PMX
///// \param selected: arreglo de par de índices individuos en la población seleccionados para cruza
///// \param oldPopulation: población de padres
///// \param newPopulation: población de hijos
///// \param pc: probabilidad de cruza
///// \param file: archivo
void SimpleCrossover(const std::vector<std::pair<size_t, size_t>> &selected, const std::vector<individual> &oldPopulation,
                  std::vector<individual> &newPopulation, const double pc, std::ofstream &file) {
  size_t numberCrossover = 0;
  size_t newPopIndex = 0;
  for(size_t indexSelect = 0; indexSelect < selected.size(); ++indexSelect, newPopIndex +=2) {
    auto selectPair = selected[indexSelect];
    // Inicializar los hijos
    newPopulation[newPopIndex].x = oldPopulation[selectPair.first].x;
    newPopulation[newPopIndex + 1].x = oldPopulation[selectPair.second].x;
    if (flip(pc)) {
      ++numberCrossover;

      // Indices de sección
      size_t index1 = rnd(0, oldPopulation[newPopIndex].x.size());

      // Cruza
      for(size_t indexCrossover = index1; indexCrossover < newPopulation[newPopIndex].x.size(); indexCrossover++){
        // Hijo 1 toma genes del padre 2
        newPopulation[newPopIndex].x[indexCrossover] = oldPopulation[newPopIndex + 1].x[indexCrossover];
        // Hijo 2 toma genes del padre 1
        newPopulation[newPopIndex + 1].x[indexCrossover] = oldPopulation[newPopIndex].x[indexCrossover];
      }
    }
  }

  file << "Número de cruzas: " << numberCrossover << '\n';
}

///// Muta por inserción con probabilidad pm
///// \param newPopulation: población de hijos a mutar
///// \param pm: probabilidad de mutación
/// \param file: archivo donde guardar resultados
void uniformMutation(std::vector<individual> &newPopulation, const std::vector<std::pair<double, double>> &constraint, const double pm, std::ofstream &file) {
  size_t numberMutations = 0;
  for (auto &individuo: newPopulation) {
    if(flip(pm)){
      ++numberMutations;
      // Indices de sección
      size_t index1 = rnd(0, individuo.x.size() - 1);

      individuo.x[index1] = rndreal(constraint[index1].first, constraint[index1].second);
    }
  }

  file << "Número de mutaciones totales: " << numberMutations << '\n';
}

/// Tomar el más apto y ponerlo en el primer hijo
/// \param newPopulation: población de hijos
/// \param bestIndividual: individuo más apto de la población de padres
void elitism(std::vector<individual> &newPopulation, individual &bestGenerationIndividual) {
  newPopulation[0] = bestGenerationIndividual;
}

#endif //RESTRICTEDGENETICALGORITHM__UTILS_H_
