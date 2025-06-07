#include "LDA_base.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void LDAbase::read_data_common()
{
  // Read data
  Weight_Mat = Rcpp::as<Eigen::MatrixXd>(model["Weight_Mat"]);
  W = model["W"]; Z = model["Z"];
  vocab = model["vocab"];
  regular_k = model["no_keyword_topics"];
  model_fit = model["model_fit"];
  //sentNeighArr = model["sentNeighArr"];
  try {
    // Attempt to load word_neighbors from model
    word_neighbors = model["word_neighbors"];

    // Print success message
    Rcpp::Rcout << "Successfully loaded word_neighbors." << std::endl;
  } catch (std::exception& e) {
      // Catch any exceptions and print an error message
      Rcpp::Rcerr << "Error loading word_neighbors: " << e.what() << std::endl;
      throw; // Rethrow the exception if needed
  }
  
  keyword_k = 0;
  num_topics = regular_k;

  // document-related constants
  num_vocab = vocab.size();
  num_doc = W.size();
  // alpha -> specific function

  // Options
  options_list = model["options"];
  use_weights = options_list["use_weights"];
  slice_A = options_list["slice_shape"];
  store_theta = options_list["store_theta"];
  store_phi = options_list["store_phi"];
  //Rcpp::Rcout << "store_theta: " << store_theta << std::endl;
  //Rcpp::Rcout << "store_phi: " << store_phi << std::endl;

  thinning = options_list["thinning"];
  llk_per = options_list["llk_per"];
  verbose = options_list["verbose"];
  weights_type = Rcpp::as<std::string>(options_list["weights_type"]);
  estimate_phi = options_list["estimate_phi"];

  // Priors
  priors_list = model["priors"];
  beta = priors_list["beta"];
  lambda = priors_list["lambda"];

  eta_1 = priors_list["eta_1"];
  eta_2 = priors_list["eta_2"];
  eta_1_regular = priors_list["eta_1_regular"];
  eta_2_regular = priors_list["eta_2_regular"];

  if (estimate_phi == 0) {
    Phi = Rcpp::as<Eigen::MatrixXd>(priors_list["phi"]);
  }

  // Stored values
  stored_values = model["stored_values"];

  // Slice Sampling
  // this is used except cov models
  model_settings = model["model_settings"];
  min_v = model_settings["slice_min"];
  min_v = shrinkp(min_v);

  max_v = model_settings["slice_max"];
  max_v = shrinkp(max_v);
}


void LDAbase::initialize_common()
{
  // Prior values are set in `LDAbase::read_data_common()`

  // Slice sampling initialization
  max_shrink_time = 200;

  //
  // Vocabulary weights
  //
  vocab_weights = VectorXd::Constant(num_vocab, 1.0);

  int z, w;
  int doc_len;
  IntegerVector doc_z, doc_w;
  RowVectorXd vocab_weights_doc;


  // Construct vocab weights
  for (int doc_id = 0; doc_id < num_doc; ++doc_id) {
    doc_w = W[doc_id];
    doc_len = doc_w.size();
    doc_each_len.push_back(doc_len);

    for (int w_position = 0; w_position < doc_len; ++w_position) {
      w = doc_w[w_position];
      vocab_weights(w) += 1.0;
    }
  }
  total_words = (int)vocab_weights.sum();


  //
  // Construct data matrices
  //
  n_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_dk_noWeight = MatrixXd::Zero(num_doc, num_topics);
  n_kv_noWeight = MatrixXd::Zero(num_topics, num_vocab);
  n_k = VectorXd::Zero(num_topics);

  total_words_weighted = 0.0;
  double temp;
  for(int doc_id = 0; doc_id < num_doc; ++doc_id){
    doc_z = Z[doc_id], doc_w = W[doc_id];
    doc_len = doc_each_len[doc_id];
    vocab_weights_doc = Weight_Mat.row(doc_id);
    for(int w_position = 0; w_position < doc_len; ++w_position){
      z = doc_z[w_position], w = doc_w[w_position];

      n_kv(z, w) += vocab_weights_doc(w);
      n_k(z) += vocab_weights_doc(w);
      n_dk(doc_id, z) += vocab_weights_doc(w);
      n_dk_noWeight(doc_id, z) += 1.0;
      n_kv_noWeight(z, w) += 1.0;
    }

    temp = n_dk.row(doc_id).sum();
    doc_each_len_weighted.push_back(temp);
    total_words_weighted += temp;
  }

  // Use during the iteration
  z_prob_vec = VectorXd::Zero(num_topics);
  generateSentNeighArr(num_doc, word_neighbors);
  Z_vector = convertZ();
  /*
  try {
    // Main processing logic
    void generateSentNeighArr(num_doc, word_neighbors);
  } catch (std::exception& e) {
      Rcpp::Rcerr << "Error: " << e.what() << std::endl;
      std::vector<NeightNode>().swap(sentNeighArr); // Clear allocated memory
      throw;
  }
  */

}


void LDAbase::parameters_store(int r_index)
{
  if (store_theta)
    store_theta_iter(r_index);

  // Store phi
  if (store_phi){
    //Rcpp::Rcout << "store_phi is implemented." << std::endl;
    store_phi_iter(r_index);
  }
}


int LDAbase::sample_z(VectorXd &alpha, int z, int s, int w, int doc_id, RowVectorXd &vocab_weights_doc, MatrixXd &Phi)
{
  int new_z = -1;
  double numerator, denominator, sum, phi_kw;

  // Remove data
  n_kv(z, w) -= vocab_weights_doc(w);
  n_k(z) -= vocab_weights_doc(w);
  n_dk(doc_id, z) -= vocab_weights_doc(w);
  n_dk_noWeight(doc_id, z) -= 1;
  n_kv_noWeight(z, w) -= 1;

  VectorXd neiProb = computeNeiProb(doc_id, w, num_topics);
  //VectorXd neiProb = VectorXd::Zero(num_topics);
  if (estimate_phi == 1) {
    for (int k = 0; k < num_topics; ++k) {
      numerator = (beta + n_kv(k, w)) * (n_dk(doc_id, k) + alpha(k));
      denominator = (num_vocab * beta + n_k(k));
      z_prob_vec(k) = numerator / denominator;

      if (neiProb.sum() > 0) {
        z_prob_vec(k) *= exp(lambda*neiProb(k));
      }
    }
  } else {
    for (int k = 0; k < num_topics; ++k) {
      phi_kw = exp(Phi(k, w));
      z_prob_vec(k) = phi_kw * (n_dk(doc_id, k) + alpha(k));

      if (neiProb.sum() > 0) {
        z_prob_vec(k) *= exp(lambda*neiProb(k));
      }
    }
  }

  // Normalize and sample
  sum = z_prob_vec.sum();
  new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics);

  // Add data back
  n_kv(new_z, w) += vocab_weights_doc(w);
  n_k(new_z) += vocab_weights_doc(w);
  n_dk(doc_id, new_z) += vocab_weights_doc(w);
  n_dk_noWeight(doc_id, new_z) += 1;
  n_kv_noWeight(new_z, w) += 1;

  return new_z;
}


/*
VectorXd LDAbase::computeNeiProb(int doc_id, int word_name, int num_topics) {
    VectorXd neiProb = VectorXd::Zero(num_topics);

    // Retrieve the NeightNode for the current document
    if (doc_id >= sentNeighArr.size()) {
        Rcpp::Rcerr << "Invalid document ID: " << doc_id << std::endl;
        return neiProb;
    }

    const NeightNode& docNeightNode = sentNeighArr[doc_id];
    std::vector<int> neighbors = docNeightNode.getNei(word_name);
    IntegerVector doc_topics = Z[doc_id];

    if (!neighbors.empty()) {
        for (int neiIdx : neighbors) {
            int neighbor_topic = doc_topics[neiIdx]; // Assuming Z[doc_id] holds topic assignments
            neiProb[neighbor_topic] += 1.0;
        }
        neiProb /= neighbors.size(); // Normalize probabilities
    }

    return neiProb;
}
*/

/*
void LDAbase::generateSentNeighArr(int num_doc, List word_neighbors){
  sentNeighArr.clear(); // Clear any existing data

  for (int i = 0; i < num_doc; i++) {
    NeightNode nn;
    std::unordered_set<int> idSet; // To avoid processing duplicate IDs

    // Get the current pseData element (a vector of integers)
    IntegerVector currentData = W[i];

    // Process each ID in currentData
    for (int j = 0; j < currentData.size(); j++) {
      int id = currentData[j];

      // Skip if already processed
      if (idSet.find(id) != idSet.end()) {
        continue;
      }

      idSet.insert(id);

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
  
}
*/