#include "NeightNode.h"

void NeightNode::add(int id, int ind) {
    nei[id].push_back(ind);
}

std::vector<int> NeightNode::getNei(int id) const {
    auto it = nei.find(id);
    if (it != nei.end()) {
        return it->second;
    } else {
        return {};
    }
}
