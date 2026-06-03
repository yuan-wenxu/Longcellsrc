#pragma once

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "isoform_core.h"

namespace longcellsrc {

// ============================================================
// Data structures
// ============================================================

struct CigarResult {
  std::vector<std::string> marks;
  std::vector<int> counts;
};

using IsoformMatrix = std::vector<std::pair<int, int>>;

struct ReadsResult {
  std::vector<int> ids;
  std::vector<std::string> isoforms;
  std::vector<int> isoends;
  std::vector<std::string> polyAs;
};

struct SpliceSiteTableResult {
  std::vector<int> ids;
  std::vector<std::string> starts;
  std::vector<std::vector<int>> mid;
  std::vector<std::string> ends;
  std::vector<std::string> splice_site_names;
};

struct ShareNeighborRow {
  std::string node1;
  std::string node2;
  int share;
};

struct BarcodeMatchResult {
  std::vector<int> ids;
  std::vector<std::string> barcodes;
  std::vector<int> positions;
  std::vector<int> edits;
  std::vector<std::string> umis;
  std::vector<std::string> adapters;
};

// ============================================================
// Normal quantile
// ============================================================

double qnorm_core(double p, double mu, double sigma);

// ============================================================
// reads_extraction
// ============================================================

CigarResult cigar_process_core(const std::string& cigar);
bool cigar_check_core(const std::vector<std::string>& marks,
                      const std::vector<int>& counts);
int seq_end_core(int start, const std::vector<std::string>& marks,
                 const std::vector<int>& counts);

IsoformMatrix extract_isoform_core(int read_start,
    const std::vector<std::string>& marks, const std::vector<int>& counts,
    int refer_start, int refer_end, int flank);

IsoformMatrix chop_out_bound_core(const IsoformMatrix& isoform,
    const IsoformMatrix& annotation, const std::string& strand,
    bool& chopped);

IsoformMatrix splicesite_correct_core(const IsoformMatrix& isoform,
    const IsoformMatrix& annotation, const std::string& strand, int bin);

int isoform_end_core(const IsoformMatrix& isoform, const std::string& strand);

std::string isoform_to_string_core(const IsoformMatrix& isoform,
                                   const std::string& sep = "|");

ReadsResult extract_reads_core(const std::vector<std::string>& seqs,
    const std::vector<std::string>& cigars,
    const std::vector<int>& positions, const IsoformMatrix& annotation,
    const std::string& strand, int toolkit, int end_flank,
    int splice_site_bin);

// ============================================================
// reads_filter
// ============================================================

std::vector<double> size_filter_core(const std::vector<double>& size,
                                     double ratio = 0.1);

// ============================================================
// splice_site_correct
// ============================================================

SpliceSiteTableResult splice_site_table_core(
    const std::vector<std::string>& isoforms, const std::string& split,
    const std::string& sep, int splice_site_thresh);

std::vector<std::vector<bool>> matrix_xor_core(
    const std::vector<std::vector<int>>& mat);

// ============================================================
// exon_corres
// ============================================================

int bin_sum_core(const std::vector<std::pair<int, int>>& bins);

std::string bin2exonid_core(const std::string& bin_str, int status,
    const std::vector<int>& start, const std::vector<int>& end,
    const std::vector<std::string>& exon_id, int mid_bias, int end_bias,
    int end_overlap, const std::string& nonsense_label,
    const std::string& split, const std::string& sep);

std::map<std::string, std::string> bins2exonids_core(
    const std::vector<std::string>& bins, const std::vector<int>& status,
    const std::vector<int>& start, const std::vector<int>& end,
    const std::vector<std::string>& exon_id, int mid_bias, int end_bias,
    int end_overlap, const std::string& nonsense_label,
    const std::string& split, const std::string& sep);

std::map<std::string, std::string> isos2exonids_index_core(
    const std::vector<std::string>& isoforms,
    const std::vector<int>& start, const std::vector<int>& end,
    const std::vector<std::string>& exon_id, int mid_bias, int end_bias,
    int end_overlap, const std::string& nonsense_label,
    const std::string& split, const std::string& sep);

// ============================================================
// umi_dist
// ============================================================

std::pair<double, double> iso2_mid_diff_core(const std::string& a,
    const std::string& b, int end_bias,
    const std::string& split = "|", const std::string& sep = ",");

int needle_core(const std::string& A, const std::string& B,
    int match_score = 1, int mismatch_score = -1, int gap_score = -1);

std::vector<std::vector<int>> index_core(
    const std::vector<std::string>& data,
    const std::vector<std::string>& uniq);

std::vector<std::vector<int>> umi_graph_table_core(
    const std::vector<std::string>& umi,
    const std::vector<std::string>& isoform,
    const std::vector<int>& count, int sim_thresh, int iso_thresh,
    const std::string& split, const std::string& sep);

// ============================================================
// umi_cluster
// ============================================================

std::vector<ShareNeighborRow> share_neighbor_core(
    const std::vector<std::string>& index,
    const std::vector<std::vector<std::string>>& neighbor,
    const std::vector<int>& count);

// ============================================================
// bc
// ============================================================

std::vector<int> pos_filter_core(const std::vector<double>& start,
    const std::vector<int>& edit, std::string& warning_msg);

BarcodeMatchResult barcode_match_core(
    const std::vector<std::string>& seq,
    const std::vector<std::string>& barcodes, double mu, double sigma,
    double sigma_start, int k, int batch, int top, double cos_thresh,
    double alpha, int edit_thresh, int UMI_len, int flank);

}  // namespace longcellsrc

// ============================================================
// Re-declarations of pure-C++ helpers from headers that
// otherwise pull in Rcpp.  In the global namespace to match
// their original definitions.
// ============================================================

// edit.h
std::vector<std::string> reverse(std::vector<std::string> s);
std::set<std::string> kmer(std::vector<std::string> s, int k, int step);
std::vector<std::string> kmer(std::string s, int k, int step);
int editDist(std::string word1, std::string word2);
std::pair<int, int> traceback(std::vector<std::vector<int>> dp);
std::pair<int, int> minEditDist(std::string seq, std::string barcode);

// tag_extraction.h
int baseCount(std::string seq, char base);
std::string reverseComplement(const std::string& sequence);
bool polyADetect(const std::string& seq, int bin = 20, int count = 15,
                 char base = 'A');
size_t polyARm(std::string seq, const int polyA_len = 10);
int strSlideSearch(std::string seq, const std::string adapter,
                   const int window = 12, const int step = 3,
                   const bool first = true);
std::vector<std::string> strSubset(std::string str, const int window,
                                   const int step);
std::string replicate(std::string mode, int times);

// bc.h
double mean(std::vector<double> num);
double var(std::vector<double> num);
double update_sigma(std::vector<double> num, double sigma_start);
double update_prob(std::vector<double> m, double n);
std::vector<int> kmer_include(std::string seq, std::set<std::string> dic);
std::vector<std::vector<int>> barcodes_cos_vec(
    std::vector<std::string> barcodes, std::set<std::string> dic);
double cos_sim(std::vector<int> a, std::vector<int> b);
bool cos_sim_comp(std::pair<int, double> cos1, std::pair<int, double> cos2);
std::vector<std::string> barcode_cand_cos(
    std::string seq, std::vector<std::string> barcodes,
    std::set<std::string> dic, std::vector<std::vector<int>> index, int top = 8,
    double thresh = 0.25);
std::vector<std::vector<std::string>> NeighborExtract(
    std::vector<std::string> reads, std::vector<int> id,
    std::vector<int> pos, const int UMI_len, const int flank,
    const int bar_len);

// splice_site_correct.h
std::vector<std::string> splice_site_cpp(const std::string& input,
                                         const std::string& delimiters);
std::map<std::string, int> isoform_count(std::vector<std::string> isoform,
                                         const std::string sep);
std::unordered_map<std::string, int> splice_site_count_cpp(
    std::vector<std::string> isoform, const std::string split,
    const std::string sep);
std::string getSubString(const std::string& strValue,
                         const std::string& startChar,
                         const std::string& endChar);
bool isin(const std::string element, const std::vector<std::string>& vec);

// exon_corres.h
std::pair<int, int> exonstr2bin(std::string exon,
                                const std::string& delimiters);
std::vector<std::pair<int, int>> isostr2bins(std::string isoform,
                                             const std::string& delimiters);
std::string paste(std::vector<std::string> s, std::string sep = "|");
std::vector<std::string> exon_status(std::vector<std::string> exons,
                                     std::string split = "|");

// umi_dist.h
int iso_len(std::vector<int> sites);
int bin2_intersect(int a_start, int a_end, int b_start, int b_end);
std::vector<int> sites_chop(std::vector<int> sites, int start, int end);
int iso2_mid_dist(std::string a, std::string b, const std::string split = "|",
                  const std::string sep = ",");
int minEditDis(std::string seq1, std::string seq2, int k = 10);
int pair2id(int x, int y);
std::vector<std::string> vec_extract(std::vector<std::string> data,
                                     std::vector<int> index);
std::vector<std::vector<int>> umi_needle(std::vector<std::string> umi,
                                         std::vector<int> count,
                                         std::vector<int> umi_id1,
                                         std::vector<int> umi_id2,
                                         int thresh = 5);
std::vector<std::vector<int>> umi_edit(std::vector<std::string> umi,
                                       std::vector<int> umi_id1,
                                       std::vector<int> umi_id2, int thresh = 2,
                                       int k = 10);
void initialize(std::vector<std::string> umi);
std::vector<std::string> flatten(
    std::vector<std::vector<std::string>> const& vec);
std::vector<std::string> str_split(std::string s, std::string split);

// umi_cluster.h
int selectSum(std::map<std::string, int> count, std::vector<std::string> id);
