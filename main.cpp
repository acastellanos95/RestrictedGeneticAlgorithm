#include <iostream>
#include <limits>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "Utils.h"

int main(int argc, char *argv[]) {
  double pm, pc;
  size_t problemIndex;
  size_t Gmax, populationSize;
  std::string reportFilename;
  if (argc == 1) {
//    std::cout << "Introduzca nombre de archivo que contiene las matrices: ";
//    std::cin >> matrixFilename;
    std::cout << "Introduzca problema (1,2,3): ";
    std::cin >> problemIndex;
    std::cout << "Introduzca la probabilidad de mutación: ";
    std::cin >> pm;
    std::cout << "Introduzca la probabilidad de cruza: ";
    std::cin >> pc;
    std::cout << "Introduzca el número máximo de generaciones: ";
    std::cin >> Gmax;
    std::cout << "Introduzca el tamaño máximo de población: ";
    std::cin >> populationSize;
    std::cout << "Introduzca la semilla: ";
    std::cin >> Rseed;
    std::cout << "Nombrar archivo con reporte: ";
    std::cin >> reportFilename;
  } /*else if (argc == 8) {
    pm = std::strtod(argv[1], nullptr);
    pc = std::strtod(argv[2], nullptr);
    pi = std::strtod(argv[3], nullptr);
    Gmax = std::strtoul(argv[4], nullptr, 0);
    populationSize = std::strtoul(argv[5], nullptr, 0);
    reportFilename = std::string(argv[6]);
    Rseed = std::strtof(argv[7], nullptr);
  }*/else {
    throw std::runtime_error("programa no admite entradas en línea de comandos");
  }

  // Restricciones del tipo a_i <= x_i <= b_i de la función por problema
  std::vector<std::pair<double, double>> constraints;
  switch (problemIndex) {
    case 1:
      constraints = std::vector<std::pair<double, double>>(20,{0.0,10.0});
      break;
    case 2:
      constraints = std::vector<std::pair<double, double>> {{0.0,1200.0},{0.0,1200.0},{-0.55,0.55},{-0.55,0.55}};
      break;
    case 3:
      constraints = std::vector<std::pair<double, double>> {{-2.3,2.3},{-2.3,2.3},{-3.2,3.2},{-3.2,3.2},{-3.2,3.2}};
      break;
    default:
      throw std::runtime_error("No se eligó un problema");
      break;
  }

  /* variables de sigmoide */
  /* Que tan empinado, entre más alto el valor más empinado */
  double c1 = 0.5;
  /* Alrededor de que punto queremos que empiece a empinar, elegimos que en la primera tercera parte de la ejecución */
  double c2 = 0.3*round((double) Gmax);
  /* Altura de sigmoide, en este caso queremos que empiece en 1 y avance hasta 2 */
  double c3 = 1.0;

  std::vector<double> sigmoidParameters{c1,c2,c3};



  // Inicializa población y generador de números aleatorios
  std::ofstream file(reportFilename);
  randomize();
  std::vector<individual> oldPopulation(populationSize);
  initalizePopulation(oldPopulation, constraints);
  std::vector<individual> newPopulation(populationSize);

  std::cout << "Mejor solución conocida para función 2: " << function({679.9453, 1026.067, 0.1188764, -0.3962336}, 2) << "\n";
  std::cout << "Mejor solución conocida para función 3: " << function({-1.717143, 1.595709, 1.827247, -0.7636413, -0.763645}, 3) << "\n";

  // Reporte Inicial
  file << std::fixed;
  file << std::setprecision(5);
  file << "probabilidad de mutación: " << pm << ", probabilidad de cruza: " << pc << ", problema: " << problemIndex << ", número de generaciones: "
       << Gmax << ", tamaño de población: " << populationSize << ", semilla: " << Rseed << '\n';

  // Rondas aleatorias para inicializar máximos y mínimos
  std::vector<double> x(constraints.size(), 0.0);
  for(size_t indexRounds = 0; indexRounds < 5000; ++indexRounds){
    for(size_t indexChromosome = 0; indexChromosome < x.size(); ++indexChromosome){
      x[indexChromosome] = rndreal(constraints[indexChromosome].first, constraints[indexChromosome].second);
    }
    restrictionsValueNormalized(x, problemIndex);
  }

  // Generaciones y variables globales para generación
  // Empezamos con un individuo factible
  std::pair<size_t, individual> bestIndividual;
  // Gráfica de convergencia
  std::vector<std::pair<size_t, double>> convergence;
  for (size_t generationIndex = 1; generationIndex <= Gmax; ++generationIndex) {
    file << "--------------------- Generación " << generationIndex << " ---------------------\n";

    // Cálculo de aptitud
    fitnessCalculation(oldPopulation, sigmoidParameters, problemIndex, generationIndex);

    // Selección
    auto selected = selection(oldPopulation);

    /* Cálculo de estadísticas de aptitud de la generación y el mejor individuo de la generación */
    // El mejor individuo es aquel que es factible y con la mayor aptitud
    std::vector<individual> feasiblePopulation;

    std::copy_if(oldPopulation.begin(), oldPopulation.end(), std::back_inserter(feasiblePopulation), [](const individual &i){ return i.isFeasible; });

    individual generationBestFeasibleIndividual;

    if(feasiblePopulation.empty()){
      // Si no hubo factibles es que hubo un ajuste de restricciones
      generationBestFeasibleIndividual = bestIndividual.second;
    } else {
      generationBestFeasibleIndividual = *std::max_element(feasiblePopulation.begin(), feasiblePopulation.end(),
                                                           [](const individual &i, const individual &j) {
                                                             return i.fitness < j.fitness;
                                                           });
    }

    // Estadísticas de aptitud de los factibles
    std::stringstream ssFeasible;
    if(feasiblePopulation.size() > 0){
      double avgfitnessFeasible = avgFitnessPopulation(feasiblePopulation);
      double minfitnessFeasible = std::min_element(feasiblePopulation.begin(), feasiblePopulation.end(),
                                                   [](const individual& i, const individual& j) {
                                                     return i.fitness < j.fitness;
                                                   })->fitness;
      double maxfitnessFeasible = std::max_element(feasiblePopulation.begin(), feasiblePopulation.end(),
                                                   [](const individual& i, const individual& j) {
                                                     return i.fitness < j.fitness;
                                                   })->fitness;
      ssFeasible << "Media de aptitud de población factible: " << avgfitnessFeasible << '\n';
      ssFeasible << "Aptitud máxima de población factible: " << maxfitnessFeasible << '\n';
      ssFeasible << "Aptitud mínima de población factible: " << minfitnessFeasible << '\n';
    } else {
      ssFeasible << "No hubo factibles en esta población\n";
    }

    // Estadísticas de aptitud de todos
    double avgfitness = avgFitnessPopulation(oldPopulation);
    double minfitness = std::min_element(oldPopulation.begin(), oldPopulation.end(),
                                         [](const individual& i, const individual& j) {
                                           return i.fitness < j.fitness;
                                         })->fitness;
    double maxfitness = std::max_element(oldPopulation.begin(), oldPopulation.end(),
                                         [](const individual& i, const individual& j) {
                                           return i.fitness < j.fitness;
                                         })->fitness;

    // Revisar si el mejor de la generación es el mejor de todas las generaciones
    // Actualizar mejor individuo global para comparar en caso de que maxF haya cambiado
    std::vector<individual> tmp {bestIndividual.second};
    fitnessCalculation(tmp, sigmoidParameters, problemIndex, generationIndex);
    bestIndividual.second = tmp[0];
    if (bestIndividual.second.fitness < generationBestFeasibleIndividual.fitness && generationBestFeasibleIndividual.isFeasible){
      bestIndividual.second = generationBestFeasibleIndividual;
      bestIndividual.first = generationIndex;
    }

    std::stringstream ssGeneration;
    ssGeneration << "(";
    for(auto &xi: generationBestFeasibleIndividual.x){
      ssGeneration << xi << " ";
    }
    ssGeneration << ")";

    std::stringstream ssBestOfAll;
    ssBestOfAll << "(";
    for(auto &xi: bestIndividual.second.x){
      ssBestOfAll << xi << " ";
    }
    ssBestOfAll << ")";

    file << "Media de aptitud de población: " << avgfitness << '\n';
    file << "Aptitud máxima de población: " << maxfitness << '\n';
    file << "Aptitud mínima de población: " << minfitness << '\n';
    file << ssFeasible.rdbuf();
    file << "Mejor individuo de la generación con valor x: " << ssGeneration.rdbuf()
         << ", fitness: " << generationBestFeasibleIndividual.fitness << ", es factible?: " << std::boolalpha << generationBestFeasibleIndividual.isFeasible <<'\n';
    file << "Mejor individuo global con valor x: " << ssBestOfAll.rdbuf() <<", fitness: "
         << bestIndividual.second.fitness << ", f(x): " << bestIndividual.second.Fx << ", generación: " << bestIndividual.first << '\n';
    // Gráfica de convergencia
    convergence.emplace_back(generationIndex, generationBestFeasibleIndividual.fitness);

    // Cruza
    SimpleCrossover(selected, oldPopulation, newPopulation, pc, file);

    // Mutación
    uniformMutation(newPopulation, constraints, pm, file);

    // Elitismo
    elitism(newPopulation, bestIndividual.second);

    auto tmpPopulation = oldPopulation;
    oldPopulation = newPopulation;
    newPopulation = tmpPopulation;
  }

  // Reporte final
  file << "f(x): " << bestIndividual.second.Fx << '\n';
  file << "$(";
  for(size_t i = 0; i < bestIndividual.second.x.size(); ++i){
    file << bestIndividual.second.x[i] << ",)"[i == (bestIndividual.second.x.size() - 1)];
  }
  file << "$ & $" << bestIndividual.second.fitness << "$\\\\\n";
//   Punto 5 convergencia de la mediana
//  std::ofstream fileConvergence("convergenceTai30Mediana.dat");
//  for (auto &pair: convergence) {
//    fileConvergence << pair.first << ", " << std::setprecision(10) << std::to_string(pair.second) << '\n';
//  }
//  fileConvergence.close();
  file.close();
  return 0;
}

