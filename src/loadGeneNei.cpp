#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
List load_wordid_neighcnt(List neiPath, CharacterVector vocab) {
    // Step 1: Create a mapping for vocab words to their indices
    std::unordered_map<std::string, int> vocab_map;
    for (int i = 0; i < vocab.size(); i++) {
        vocab_map[as<std::string>(vocab[i])] = i;
    }

    // Step 2: Initialize the word_neighbors array
    int V = vocab.size();
    std::vector<std::vector<int>> word_neighbors(V); // V rows, each with an empty vector

    // Step 3: Populate word_neighbors from neiPath
    CharacterVector names = neiPath.names(); // Extract the names
    for (int i = 0; i < neiPath.size(); i++) {
        std::string word = as<std::string>(names[i]); // Access the i-th name
        if (vocab_map.find(word) != vocab_map.end()) {
            int wordIndex = vocab_map[word];

            CharacterVector neighbors = neiPath[i];
            for (int j = 0; j < neighbors.size(); j++) {
                std::string neighbor = as<std::string>(neighbors[j]);
                if (vocab_map.find(neighbor) != vocab_map.end()) {
                    word_neighbors[wordIndex].push_back(vocab_map[neighbor]);
                }
            }
        }
    }



    // Step 4: Convert word_neighbors to an Rcpp List for R compatibility
    List result(V);
    for (int i = 0; i < V; i++) {
        result[i] = wrap(word_neighbors[i]);
    }

    return result;
}
