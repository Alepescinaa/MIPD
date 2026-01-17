#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List process_patient(
    int i,
    NumericVector time,
    NumericVector hazard12,
    NumericVector hazard13,
    NumericVector h12,
    NumericVector lp12,
    NumericVector lp13,
    NumericMatrix c1,
    NumericMatrix c2,
    DataFrame patient_data,
    int m,
    double onset_rate
) {
  int n = time.size();
  NumericVector S1(n), S2(n), p(n), P(n);

  double exp_lp12 = std::exp(lp12[i-1]);
  double exp_lp13 = std::exp(lp13[i-1]);


  for (int j = 0; j < n; ++j) {
    S1[j] = std::exp(-hazard12[j] * exp_lp12 - hazard13[j] * exp_lp13);
    S2[j] = std::exp(-hazard13[j] * exp_lp13);
    p[j] = S1[j] * h12[j] / S2[j];
  }

  P[0] = 0;
  for (int j = 1; j < n; ++j) {
    double incr = p[j - 1] * (time[j] - time[j - 1]);
    P[j] = std::min(1.0, P[j - 1] + incr);
  }

  NumericVector da(m, NA_REAL);
  NumericVector ds(m, NA_REAL);

  NumericVector age = patient_data["age"];
  IntegerVector onset = patient_data["onset"];
  IntegerVector visits = patient_data["visits"];
  LogicalVector dead_vec = patient_data["dead"];
  NumericVector death_time_vec = patient_data["death_time"];
  NumericVector onset_age_vec = patient_data["onset_age"];

  bool has_onset = false;
  int onset_index = -1;
  for (int j = 0; j < onset.size(); ++j) {
    if (onset[j] == 1) {
      has_onset = true;
      onset_index = j;
      break;
    }
  }
  // Rcpp::Rcout << "onset_index = " << onset_index << std::endl;

  if (has_onset && onset_index > 0) {
    double a = age[onset_index - 1];
    double b = age[onset_index];


    std::vector<int> idx;
    for (int j = 0; j < n; ++j) {
      if (time[j] > a && time[j] < b) {
        idx.push_back(j);
      }
    }

    int nF = idx.size();
    if (nF > 0) {
      double minP = P[idx[0]];
      double maxP = P[idx[nF - 1]];

      if (minP != maxP) {
        NumericVector q(nF);
        for (int j = 0; j < nF; ++j)
          q[j] = (P[idx[j]] - minP) / (maxP - minP);

        for (int j = 0; j < m; ++j) {
          double x = c1(i-1, j);
          double selected_time = a;
          for (int k = 0; k < nF; ++k) {
            if (q[k] <= x) selected_time = time[idx[k]];
          }
          da[j] = selected_time;
        }
      } else {
        double mean_ab = (a + b) / 2.0;
        std::fill(da.begin(), da.end(), mean_ab);
      }
    } else {
      double mean_ab = (a + b) / 2.0;
      std::fill(da.begin(), da.end(), mean_ab);
    }

    std::fill(ds.begin(), ds.end(), 1);
  } else {
    int n_visits = visits.size();
    double a = age[n_visits - 1];
    double b = death_time_vec[0];

    // Rcpp::Rcout << "a = " << a << std::endl;
    // Rcpp::Rcout << "b = " << b << std::endl;

    std::vector<int> idx;
    for (int j = 0; j < n; ++j) {
      if (time[j] > a && time[j] < b) {
        idx.push_back(j);
      }
    }

    int nF = idx.size();
    double pD;
    NumericVector q(nF);

    if (nF > 0) {
      double minP = P[idx[0]];
      double maxP = P[idx[nF - 1]];

      double S1t = S1[idx[nF - 1]];
      double S2t = S2[idx[nF - 1]];
      double S2h23P = S2t * std::exp(lp13[i-1] ) * (maxP - minP);
      pD = S2h23P / (S1t + S2h23P);

      // Rcpp::Rcout << "S2t = " << S2t << std::endl;
      // Rcpp::Rcout << "S1t = " << S1 << std::endl;
      // Rcpp::Rcout << "pD = " << pD << std::endl;

      if (minP != maxP) {
        for (int j = 0; j < nF; ++j)
          q[j] = (P[idx[j]] - minP) / (maxP - minP);
      } else {
        std::fill(q.begin(), q.end(), 0.0);
      }
    } else {
      pD = (a < b) ? onset_rate : 0.0;
    }

    for (int j = 0; j < m; ++j) {
      if (c2(i-1, j) <= pD) {
        ds[j] = 1;
        if (nF > 0) {
          double x = c1(i-1, j);
          double selected_time = a;
          for (int k = 0; k < nF; ++k) {
            if (q[k] <= x) selected_time = time[idx[k]];
          }
          da[j] = selected_time;
        } else {
          da[j] = (a + b) / 2.0;
        }
      } else {
        ds[j] = 0;
        da[j] = onset_age_vec[0];
      }
    }
  }

  return List::create(
    Named("age") = da,
    Named("status") = ds
  );
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
