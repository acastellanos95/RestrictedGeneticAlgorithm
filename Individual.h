//
// Created by andre on 7/9/22.
//

#ifndef RESTRICTEDGENETICALGORITHM__INDIVIDUAL_H_
#define RESTRICTEDGENETICALGORITHM__INDIVIDUAL_H_

#include <vector>
#include <string>

/**
 * Estructura para representar un individuo
 */
struct individual {
  std::vector<double> x; /* Cromosoma es el vector x */
  double fitness; /* Valor de aptitud */
  double Fx;
  double penaltiesValueNormalized;
  bool isFeasible; /* Es factible? */
};

#endif //RESTRICTEDGENETICALGORITHM__INDIVIDUAL_H_
