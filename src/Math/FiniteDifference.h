#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

template <typename T>
T fd2ndDeriv2ndOrderCentral(const T& phiMinus1,
                            const T& phi,
                            const T& phiPlus1,
                            double deltaMinus1,
                            double deltaPlus1);

template <typename T>
T fd1stDeriv2ndOrderCentral(const T& phiMinus1,
                            const T& phi,
                            const T& phiPlus1,
                            double deltaMinus1,
                            double deltaPlus1);

template <typename T>
T fd1stDeriv2ndOrderBackwards(const T& phiMinus2,
                              const T& phiMinus1,
                              const T& phi,
                              double deltaMinus2,
                              double deltaMinus1);

template <typename T>
T fd1stDeriv2ndOrderForwards(const T& phi,
                             const T& phiPlus1,
                             const T& phiPlus2,
                             double deltaPlus1,
                             double deltaPlus2);

#include "FiniteDifferenceI.h"

#endif
