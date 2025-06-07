#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace Rcpp;

// Define the NeightNode class
class NeightNode {
public:
    std::unordered_map<int, std::vector<int>> nei;

    // Add a neighbor
    void add(int id, int ind) {
        nei[id].push_back(ind);
    }

    // Get neighbors for a given ID
    std::vector<int> getNei(int id) const {
        auto it = nei.find(id);
        if (it != nei.end()) {
            return it->second;
        } else {
            return {};
        }
    }
};

// [[Rcpp::export]]
List computNeiIndex(List pseData, IntegerVector num_neighbors_foreach_word, List word_neighbors) {
    std::vector<NeightNode> sentNeighArr;

    // Loop over pseData (equivalent to Java's outer loop)
    for (int i = 0; i < pseData.size(); i++) {
        NeightNode nn;
        std::unordered_set<int> idSet; // To avoid processing duplicate IDs

        // Get the current pseData element (a vector of integers)
        IntegerVector currentData = pseData[i];

        // Process each ID in currentData
        for (int j = 0; j < currentData.size(); j++) {
            int id = currentData[j];

            // Skip if already processed
            if (idSet.find(id) != idSet.end()) {
                continue;
            }

            idSet.insert(id);

            // Get the number of neighbors for this ID
            int neiNum = num_neighbors_foreach_word[id];
            if (neiNum < 1) {
                continue;
            }

            // Get neighbors for this ID
            IntegerVector neighbors = word_neighbors[id];

            // Process each neighbor
            for (int n = 0; n < neighbors.size(); n++) {
                int neiId = neighbors[n];

                // Check if neiId exists in currentData
                for (int ind = 0; ind < currentData.size(); ind++) {
                    if (currentData[ind] == neiId) {
                        nn.add(id, ind);
                    }
                }
            }
        }

        // Add the NeightNode to sentNeighArr
        sentNeighArr.push_back(nn);
    }

    // Convert sentNeighArr to an Rcpp List
    List result(sentNeighArr.size());
    for (size_t i = 0; i < sentNeighArr.size(); i++) {
        List node;
        for (const auto& pair : sentNeighArr[i].nei) {
            node[std::to_string(pair.first)] = pair.second;
        }
        result[i] = node;
    }

    return result;
}
