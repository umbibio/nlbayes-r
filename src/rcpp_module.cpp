// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// rcpp_module.cpp: Rcpp R/C++ interface class library -- Rcpp Module examples
//
// Copyright (C) 2010 - 2012  Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include "rcpp_module.h"
// #include <boost/format.hpp>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


ModelORNOR::ModelORNOR() {
    Rcpp::Rcout << std::endl << "This is NLBayes with model OR-NOR ..." << std::endl;
}
ModelORNOR::~ModelORNOR() {
    delete this->wrapped;
}

void ModelORNOR::set_config(unsigned int n_graphs, double s_leniency, double t_alpha, double t_beta, double zy, double zn) {
    this->n_graphs = n_graphs;
    this->s_leniency = s_leniency;
    this->t_alpha = t_alpha;
    this->t_beta = t_beta;
    this->zy = zy;
    this->zn = zn;
    Rcpp::Rcout << "Model configured ..." << std::endl;
};

void ModelORNOR::load_network(Rcpp::ListOf<Rcpp::IntegerVector> input_network) {

    if (this->interaction_network.size() > 0)
        this->interaction_network.clear();
    
    Rcpp::CharacterVector src_names = input_network.names();
    for(int i = 0; i < input_network.size(); i++) {

        Rcpp::IntegerVector   trg_list  = input_network[i];
        Rcpp::CharacterVector trg_names = input_network[i].names();

        for(int j = 0; j < trg_list.size(); j++) {

            nlb::src_trg_pair_t src_trg_pair = nlb::src_trg_pair_t(src_names[i], trg_names[j]);
            nlb::network_edge_t edge         = nlb::network_edge_t(src_trg_pair, trg_list[j]);

            this->interaction_network.push_back(edge);
        }
    }

    // this->interaction_network = interaction_network;
    Rcpp::Rcout << "Loaded network with " << this->interaction_network.size() << " edges ..." << std::endl;

};

void ModelORNOR::print_network() {
    unsigned int n_edges = this->interaction_network.size();
    int mor;
    std::string src, trg;
    nlb::src_trg_pair_t src_trg_pair;
    for(unsigned int i = 0; i < n_edges; i++) {
        tie(src_trg_pair, mor) = this->interaction_network[i];
        tie(src, trg) = src_trg_pair;
        // Rcpp::Rcout << "| " << boost::str(boost::format("% 4d") % i) << " | " << src << " ---> " << trg << " | " << boost::str(boost::format("% 2d") % mor) << " | ";
        Rcpp::Rcout << "| " << i << " | " << src << " ---> " << trg << " | " << mor << " | ";
        if (this->active_tf_set.count(src)) Rcpp::Rcout << "*";
        Rcpp::Rcout << std::endl;
    }

}

void ModelORNOR::load_evidence(Rcpp::RObject input_obj) {

    nlb::evidence_dict_t evidence;

    int dex, n_upreg = 0, n_dwreg = 0;
    std::string trg;

    Rcpp::CharacterVector src_names;
    Rcpp::CharacterVector trg_names;

    Rcpp::IntegerVector input_evidence;

    if(Rcpp::is<Rcpp::IntegerVector>(input_obj) || Rcpp::is<Rcpp::NumericVector>(input_obj)) {
        input_evidence = Rcpp::as<Rcpp::IntegerVector>(input_obj);

        trg_names = input_evidence.names();
        for(int i = 0; i < input_evidence.size(); i++) {
            trg = trg_names[i];
            dex = input_evidence[i];
            evidence[trg] = dex;

            if (dex > 0) n_upreg += 1;
            if (dex < 0) n_dwreg += 1;
        }

        Rcpp::Rcout << "Loaded evidence: " << evidence.size() << ". ";
        Rcpp::Rcout << n_upreg << " up, " << n_dwreg << " down " << std::endl;
    }
    this->evidence = evidence;

};

void ModelORNOR::load_active_tf_set(Rcpp::RObject input_obj) {


    nlb::prior_active_tf_set_t active_tf_set;
    std::string src;

    Rcpp::StringVector input_active_tf_set;

    if(Rcpp::is<Rcpp::StringVector>(input_obj)) {
        input_active_tf_set = Rcpp::as<Rcpp::StringVector>(input_obj);
        Rcpp::Rcout << "active tf set size: " << input_active_tf_set.size();
        for(int i = 0; i < input_active_tf_set.size(); i++) {
            src = input_active_tf_set[i];
            active_tf_set.insert(src);
        }
    }
    this->active_tf_set = active_tf_set;
};

void ModelORNOR::build_model() {
    double s_leniency = this->s_leniency;
    const double SPRIOR [3 * 3] = {1.0 - s_leniency, 0.9 * s_leniency, 0.1 * s_leniency,
                                   0.5 * s_leniency, 1.0 - s_leniency, 0.5 * s_leniency,
                                   0.1 * s_leniency, 0.9 * s_leniency, 1.0 - s_leniency};

    double zy_value = this->zy == 0. ? 1. : this->zy;
    double zn_value;
    if(this->evidence.size() > 0 && this->zn == 0.) {
        unsigned int n_edges = this->interaction_network.size();
        unsigned int n_edges_deg = 0;

        int mor;
        std::string src, trg;
        nlb::src_trg_pair_t src_trg_pair;
        for(auto edge: this->interaction_network) {
            tie(src_trg_pair, mor) = edge;
            tie(src, trg) = src_trg_pair;
            if(evidence.find(trg) != evidence.end() && evidence[trg] != 0)
                n_edges_deg += 1;
        }
        zn_value = (double) n_edges_deg / n_edges / 10;
        // zn_value = std::pow((double) n_edges_deg / n_edges, 1.5);
    } else if (this->zn == 0.) {
        zn_value = 1.;
    } else {
        zn_value = this->zn;
    }

    if (this->verbosity >= 2) {
        Rcpp::Rcout << "\tT alpha   : " << this->t_alpha << std::endl;
        Rcpp::Rcout << "\tT beta    : " << this->t_beta << std::endl;
        Rcpp::Rcout << "\tS leniency: " << s_leniency << std::endl;
        Rcpp::Rcout << "\tZY value  : " << zy_value << std::endl;
        Rcpp::Rcout << "\tZN value  : " << zn_value << std::endl;
        Rcpp::Rcout << "\t# Graphs  : " << this->n_graphs << std::endl;
    }

    this->wrapped = new nlb::ModelORNOR(this->interaction_network, this->evidence, this->active_tf_set, SPRIOR,
                                        this->t_alpha, this->t_beta, zy_value, zn_value, this->n_graphs);

}

void ModelORNOR::load_data( Rcpp::ListOf<Rcpp::IntegerVector> input_network, Rcpp::RObject input_evidence, Rcpp::RObject input_active_tf_set) {
    this->load_network(input_network);
    this->load_evidence(input_evidence);
    this->load_active_tf_set(input_active_tf_set);
    this->build_model();
};

int ModelORNOR::sample_posterior(unsigned int N, unsigned int dN, double gr_level) {

    Rcpp::Rcout << std::endl;
    double gr = INFINITY;
    Progress p(N, true);

    int status = 0;
    unsigned int n = 0;
    for (
        ;
        n < N && gr > gr_level;
        dN=std::min(dN, N-n), n+=dN, gr=this->wrapped->get_max_gelman_rubin()
    ) {
        if (Progress::check_abort()) {
            status = -1;
            break;
        }
        this->wrapped->sample_n(dN);
        p.increment(dN);
    }

    bool converged = gr <= gr_level;

    if (converged) {
        p.increment(N - n);
        Rcpp::Rcout << "Converged after " << this->wrapped->total_sampled << " samples. ";
    } else if (status == 0) {
        status = 1;
        Rcpp::Rcout << "Drawed " << this->wrapped->total_sampled << " samples so far. ";
    } else if (status == -1) {
        Rcpp::Rcout << "Process interrupted. " << std::endl;
        Rcpp::Rcout << "Drawed " << this->wrapped->total_sampled << " samples so far. ";
    }
    Rcpp::Rcout << "Max Gelman-Rubin statistic is " << gr;
    Rcpp::Rcout << " (target" << (converged ? " was " : " is ") << gr_level << ")" << std::endl;

    return status;
}

void ModelORNOR::burn_stats() {
    this->wrapped->burn_stats();
    this->wrapped->total_sampled = 0;
}

void ModelORNOR::set_verbosity(unsigned int v) {
    this->verbosity = v;
}

double ModelORNOR::get_gr() {
    return this->wrapped->get_max_gelman_rubin();
}

Rcpp::NumericVector ModelORNOR::get_all_gr() {
    nlb::gelman_rubin_vector_t gelman_rubin_vector = this->wrapped->get_gelman_rubin();

    std::vector<std::string> names;
    std::vector<double> values;
    names.reserve(gelman_rubin_vector.size());
    values.reserve(gelman_rubin_vector.size());
    for (auto node: gelman_rubin_vector) {
        names.push_back(node.first);
        values.push_back(node.second);
    }

    Rcpp::NumericVector output_vector(values.begin(), values.end());
    Rcpp::CharacterVector names_vector(names.begin(), names.end());

    output_vector.names() = names_vector;
    return output_vector;
}

Rcpp::DataFrame ModelORNOR::get_posterior(Rcpp::String input_var_name) {
    nlb::posterior_vector_t result = this->wrapped->get_posterior(input_var_name);
    std::string name, uid;
    std::vector<double> means, sds;
    std::vector<std::string> col_names, col_uids;
    std::vector<double> col_V1_means, col_V2_means, col_V3_means;
    std::vector<double> col_V1_sds, col_V2_sds, col_V3_sds;
    unsigned int n_stats = std::get<2>(result[0]).size();

    for (auto row: result) {
        tie(name, uid, means, sds) = row;
        col_names.push_back(name);
        col_uids.push_back(uid);

        col_V1_means.push_back(means[0]);
        if (n_stats > 1) col_V2_means.push_back(means[1]);
        if (n_stats > 2) col_V3_means.push_back(means[2]);

        col_V1_sds.push_back(sds[0]);
        if (n_stats > 1) col_V2_sds.push_back(sds[1]);
        if (n_stats > 2) col_V3_sds.push_back(sds[2]);
    }

    Rcpp::DataFrame out;
    switch(n_stats) {
        case 1:
            out = Rcpp::DataFrame::create(
                Rcpp::_["name"] = col_names,
                Rcpp::_["uid"] = col_uids,
                Rcpp::_["V1.mean"] = col_V1_means,
                Rcpp::_["V1.sd"] = col_V1_sds
            );
            break;
        case 2:
            out = Rcpp::DataFrame::create(
                Rcpp::_["name"] = col_names,
                Rcpp::_["uid"] = col_uids,
                Rcpp::_["V1.mean"] = col_V1_means,
                Rcpp::_["V2.mean"] = col_V2_means,
                Rcpp::_["V1.sd"] = col_V1_sds,
                Rcpp::_["V2.sd"] = col_V2_sds
            );
            break;
        default:
            out = Rcpp::DataFrame::create(
                Rcpp::_["name"] = col_names,
                Rcpp::_["uid"] = col_uids,
                Rcpp::_["V1.mean"] = col_V1_means,
                Rcpp::_["V2.mean"] = col_V2_means,
                Rcpp::_["V3.mean"] = col_V3_means,
                Rcpp::_["V1.sd"] = col_V1_sds,
                Rcpp::_["V2.sd"] = col_V2_sds,
                Rcpp::_["V3.sd"] = col_V3_sds
            );
    }
    
    return out;
}


RCPP_MODULE(nlbInference){
    
    Rcpp::class_<ModelORNOR>("ModelORNOR")
    .constructor()
    .method("set.config", &ModelORNOR::set_config, "Configure the model parameters")
    .method("load.network", &ModelORNOR::load_network, "Load data and build the model")
    .method("load.evidence", &ModelORNOR::load_evidence, "Load data and build the model")
    .method("load.active.tf.set", &ModelORNOR::load_active_tf_set, "Load data and build the model")
    .method("print.network", &ModelORNOR::print_network, "Print the network edges")
    .method("load.data", &ModelORNOR::load_data, "Load data and build the model")
    .method("build.model", &ModelORNOR::build_model, "build the model")
    .method("sample.posterior", &ModelORNOR::sample_posterior, "Sample from model posterior distribution")
    .method("burn.stats", &ModelORNOR::burn_stats, "Reset statistics for all random variables (chains burning)")
    .method("set.verbosity", &ModelORNOR::set_verbosity, "Set the verbosity")
    .method("get.gr", &ModelORNOR::get_gr, "Get the maximum Gelman-Rubin statistic across all random variables in the model")
    .method("get.all.gr", &ModelORNOR::get_all_gr, "Get a list of GR stats for all random variables")
    .method("get.posterior", &ModelORNOR::get_posterior, "Get statistics (mean and sd) for a variable class")
    ;
}

