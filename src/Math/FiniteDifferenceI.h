#include "FiniteDifference.h"

template <typename T>
T fd2ndDeriv2ndOrderCentral(const T& phiMinus1,
                            const T& phi,
                            const T& phiPlus1,
                            double deltaMinus1,
                            double deltaPlus1)
{
  double alpha(deltaPlus1/deltaMinus1);

  return (alpha*phiMinus1 - (1. + alpha)*phi + phiPlus1)/(0.5*alpha*(alpha + 1.)*deltaMinus1*deltaMinus1);
}

template <typename T>
T fd1stDeriv2ndOrderCentral(const T& phiMinus1,
                            const T& phi,
                            const T& phiPlus1,
                            double deltaMinus1,
                            double deltaPlus1)
{

}

template <typename T>
T fd1stDeriv2ndOrderBackwards(const T& phiMinus2,
                              const T& phiMinus1,
                              const T& phi,
                              double deltaMinus2,
                              double deltaMinus1)
{

}

template <typename T>
T fd1stDeriv2ndOrderForwards(const T& phi,
                             const T& phiPlus1,
                             const T& phiPlus2,
                             double deltaPlus1,
                             double deltaPlus2)
{

}
