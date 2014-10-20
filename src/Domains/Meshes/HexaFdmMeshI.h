#include "HexaFdmMesh.h"

template <class STATE_TYPE>
HexaFdmMesh<STATE_TYPE>::HexaFdmMesh()
{

}

template <class STATE_TYPE>
HexaFdmMesh<STATE_TYPE>::HexaFdmMesh(int nI, int nJ, int nK)
{
  
}

template <class STATE_TYPE>
void HexaFdmMesh<STATE_TYPE>::initialize(Input &input)
{

    StructuredMesh<STATE_TYPE>::initialize(input);

}

template <class STATE_TYPE>
void HexaFdmMesh<STATE_TYPE>::computeTimeDerivatives(STATE_TYPE* timeDerivatives)
{

}
