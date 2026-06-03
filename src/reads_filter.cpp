#include "umi_dist.h"
#include "reads_filter.h"
#include "isoform_core.h"



int iso2_dis(std::string a,std::string b,
             const std::string split,
             const std::string sep){
  return longcellsrc::iso2_dis_core(a, b, split, sep);
}
// [[Rcpp::export]]
DataFrame isos_dis(const std::vector<std::string> isoforms,const int thresh,
             const std::string split,const std::string sep){
  std::vector<longcellsrc::IsoDistanceRow> rows = longcellsrc::isos_dis_core(isoforms, thresh, split, sep);
  std::vector<int> node1;
  std::vector<int> node2;
  std::vector<int> dis;
  node1.reserve(rows.size());
  node2.reserve(rows.size());
  dis.reserve(rows.size());
  for (const auto& row : rows) {
    node1.push_back(row.node1 + 1);
    node2.push_back(row.node2 + 1);
    dis.push_back(row.dis);
  }
  return DataFrame::create(Named("node1") = node1, Named("node2") = node2, Named("dis") = dis);
}

// [[Rcpp::export]]
NumericVector size_filter_cpp(NumericVector size,double ratio){
  int n = size.size();
  double left = 0, right = sum(size);
  int id = 0;
  double diff = 100;
  for(int i = 0;i < n;i++){
    left += size[i];
    right -= size[i];
    double thresh = left-right*ratio/(1-ratio);
    //Rcout << left<< "-" << right*ratio/(1-ratio) << endl;
    if(diff > abs(thresh)){
      diff = abs(thresh);
      id = i;
    }
    if(thresh >= 0){
      break;
    }
  }
  NumericVector weight(n,0.0);
  weight[size > size[id]] = 1;
  weight[size == size[id]] = (sum(size <= size[id])-double(id)-1.0)/sum(size == size[id]);
  return(weight);
}

// Parse a single isoform string into numeric vector of sites
static inline void parse_sites(const char* p, std::vector<double>& v) {
  v.clear();
  double cur = 0.0;
  bool in_num = false;
  bool seen_dot = false;
  double frac = 0.1;

  for (; *p; ++p) {
    char c = *p;
    if (c == 32) continue; // skip spaces

    if ((c >= 48 && c <= 57)) { // digit
      if (!in_num) { in_num = true; cur = 0.0; seen_dot = false; frac = 0.1; }
      if (!seen_dot) {
        cur = cur * 10.0 + (c - 48);
      } else {
        cur += (c - 48) * frac;
        frac *= 0.1;
      }
    } else if (c == 46) { // dot
      if (!in_num) { in_num = true; cur = 0.0; }
      if (!seen_dot) { seen_dot = true; frac = 0.1; }
    } else {
      if (in_num) { v.push_back(cur); in_num = false; }
      // any non-numeric char (e.g., sep or split) is a delimiter
    }
  }
  if (in_num) v.push_back(cur);
}

// [[Rcpp::export]]
Rcpp::NumericVector isos_len_cpp(Rcpp::CharacterVector isos) {
  const int n = isos.size();
  Rcpp::NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    if (isos[i] == NA_STRING) { out[i] = NA_REAL; continue; }
    out[i] = longcellsrc::iso_length_core(Rcpp::as<std::string>(isos[i]));
  }
  return out;
}
