#ifndef __METHODS_H_INCLUDED__
#define __METHODS_H_INCLUDED__

#include "linear-system.h"

#include <assert.h>

namespace methods {
void Beautify(SLAE* slae /* out */);

void MakeUpperTriangle(SLAE* slae /* out */);

}; // namespace methods

#endif