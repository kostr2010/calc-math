#include "methods.h"

#include <assert.h>

namespace methods {
namespace {
void BeautifyPartial(SLAE* slae /* out */, size_t offset) {
    const auto pos = slae->FindBiggestElement(offset);

    slae->ExchangeRow(offset, pos.first);
    slae->ExchangeCol(offset, pos.second);

    if (++offset < slae->GetSize()) {
        BeautifyPartial(slae, offset);
    }

    return;
}

void MakeUpperTrianglePartial(SLAE* slae /* out */, size_t offset) {
    assert(offset + 1 < slae->GetSize());

    const auto main_elem = slae->A_.at(offset, offset);

    assert(main_elem != 0);

    for (size_t row = offset + 1; row < slae->GetSize(); row++) {
        const auto multiplier = slae->A_.at(row, offset) / main_elem;

        slae->AddRow(row, offset, 1, -1 * multiplier);
    }

    if (++offset + 1 < slae->GetSize()) {
        MakeUpperTrianglePartial(slae, offset);
    }

    return;
}
}; // namespace

void Beautify(SLAE* slae /* out */) {
    BeautifyPartial(slae, 0);

    return;
}

void MakeUpperTriangle(SLAE* slae /* out */) {
    MakeUpperTrianglePartial(slae, 0);

    return;
}

}; // namespace methods
