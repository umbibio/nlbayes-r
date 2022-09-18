#ifndef RCPPMOD_MODELORNOR
#define RCPPMOD_MODELORNOR

#include <Rcpp.h>
#include "ModelORNOR.h"


class ModelORNOR {
public:
    ModelORNOR();
    ~ModelORNOR();
    void set_config(unsigned int = 3, double = .3, double = 2., double = 2., double = 0., double = 0.);
    void load_data(Rcpp::ListOf<Rcpp::IntegerVector>, Rcpp::RObject, Rcpp::RObject);
    void load_network(Rcpp::ListOf<Rcpp::IntegerVector>);
    void load_evidence(Rcpp::RObject);
    void load_active_tf_set(Rcpp::RObject);
    void print_network();
    void build_model();
    int sample_posterior(unsigned int, unsigned int, double);
    void burn_stats();
    void set_verbosity(unsigned int);
    double get_gr();
    Rcpp::NumericVector get_all_gr();
    Rcpp::DataFrame get_posterior(Rcpp::String);
private:
    unsigned int verbosity = 0;
    unsigned int n_graphs = 3;
    double s_leniency = 0.3, t_alpha = 2., t_beta = 2., zy = 0., zn = 0.;
    nlb::network_t interaction_network;
    nlb::evidence_dict_t evidence;
    nlb::prior_active_tf_set_t active_tf_set;
    nlb::ModelORNOR* wrapped;
};

#endif
