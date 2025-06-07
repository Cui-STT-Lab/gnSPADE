#ifndef __LDAbase__INCLUDED__
#define __LDAbase__INCLUDED__
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_meta.h"
#include "NeightNode.h"
using namespace Eigen;
using namespace Rcpp;
using namespace std;

class LDAbase : virtual public keyATMmeta
{
  // Base function for the Weighted LDA models
  // This inherits keyATMmeta class.

  public:
    // Constructor
    LDAbase(List model_) :
      keyATMmeta(model_) {};

    // Variables
    MatrixXd n_kv;
    VectorXd n_k;
    VectorXd n_k_noWeight;
    //std::vector<std::vector<int>> Z_vector;
    // Functions
    // In LDA, we do not need to read and initialize X
    virtual void read_data_common() override final;
    virtual void initialize_common() override final;
    virtual void parameters_store(int r_index) override final;
    virtual int sample_z(VectorXd &alpha, int z, int s, int w, int doc_id, RowVectorXd &vocab_weights_doc, MatrixXd &Phi) override final;
    VectorXd computeNeiProb(int doc_id, int word_name, int num_topics) {
      VectorXd neiProb = VectorXd::Zero(num_topics);

      // Retrieve the NeightNode for the current document
      if (doc_id >= sentNeighArr.size()) {
          Rcpp::Rcerr << "Invalid document ID: " << doc_id << std::endl;
          return neiProb;
      }

      const NeightNode& docNeightNode = sentNeighArr[doc_id];
      std::vector<int> neighbors = docNeightNode.getNei(word_name);
      //IntegerVector doc_topics = Z[doc_id];
      const std::vector<int>& doc_topics = Z_vector[doc_id];
      if (!neighbors.empty()) {
          for (int neiIdx : neighbors) {
              int neighbor_topic = doc_topics[neiIdx]; // Assuming Z[doc_id] holds topic assignments
              neiProb[neighbor_topic] += 1.0;
          }
          neiProb /= neighbors.size(); // Normalize probabilities
      }

      return neiProb;
    };
    //virtual void generateSentNeighArr(int num_doc, List word_neighbors) override final;

    // Convert Z from Rcpp::List to std::vector<std::vector<int>>
    std::vector<std::vector<int>> convertZ() {
        std::vector<std::vector<int>> Z_vector(Z.size());
        for (size_t i = 0; i < Z.size(); ++i) {
            Rcpp::IntegerVector doc_topics = Z[i];
            Z_vector[i] = std::vector<int>(doc_topics.begin(), doc_topics.end());
        }
        return Z_vector;
    };

};

#endif

