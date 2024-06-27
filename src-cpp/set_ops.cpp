//
// Created by Achille Nazaret on 11/20/23.
//

#include "set_ops.h"

void union_with_single_element(const FlatSet &a, int b, FlatSet &out) {
    bool inserted = false;
    for (auto p: a) {
        if (p > b && !inserted) {
            out.insert(b);
            inserted = true;
            out.insert(p);
        } else {
            out.insert(p);
        }
        out.insert(p);
    }
    if (!inserted) { out.insert(b); }
}
