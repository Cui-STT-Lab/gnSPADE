#ifndef NEIGHTNODE_H
#define NEIGHTNODE_H

#include <unordered_map>
#include <vector>

class NeightNode {
public:
    std::unordered_map<int, std::vector<int>> nei;

    // Add a neighbor
    void add(int id, int ind);

    // Get neighbors for a given ID
    std::vector<int> getNei(int id) const;
};

#endif // NEIGHTNODE_H

