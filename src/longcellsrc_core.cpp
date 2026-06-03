#include "longcellsrc_core.h"

#include <algorithm>
#include <cmath>
#include <regex>
#include <sstream>
#include <stdexcept>

#include "edit.h"

// ============================================================
// Internal helpers (not declared in header)
// ============================================================
namespace {

std::vector<std::string> split_by_chars(const std::string& input,
                                        const std::string& delimiters) {
  std::vector<std::string> tokens;
  size_t start = 0, pos;
  while (start < input.size()) {
    pos = input.find_first_of(delimiters, start);
    if (pos == std::string::npos) {
      if (start < input.size()) tokens.push_back(input.substr(start));
      break;
    }
    if (pos > start) tokens.push_back(input.substr(start, pos - start));
    start = pos + 1;
  }
  return tokens;
}

std::string get_sub_string(const std::string& strValue,
                           const std::string& startChar,
                           const std::string& endChar) {
  size_t s = strValue.find(startChar);
  size_t e = strValue.rfind(endChar);
  if (s != std::string::npos && e != std::string::npos && s != e)
    return strValue.substr(s + 1, e - s - 1);
  return "";
}

double bc_mean(const std::vector<double>& num) {
  double sum = 0;
  for (auto v : num) sum += v;
  return sum / num.size();
}

}  // anonymous namespace

// ============================================================
// Global-namespace function definitions
// (declared in longcellsrc_core.h)
// ============================================================

// --- tag_extraction.h ---

int baseCount(std::string seq, char base) {
  int count = 0;
  for (char c : seq) if (c == base) count++;
  return count;
}

std::string reverseComplement(const std::string& sequence) {
  std::string comp;
  comp.reserve(sequence.length());
  for (char base : sequence) {
    switch (base) {
      case 'A': comp += 'T'; break;
      case 'T': comp += 'A'; break;
      case 'G': comp += 'C'; break;
      case 'C': comp += 'G'; break;
      default: comp += base;
    }
  }
  std::reverse(comp.begin(), comp.end());
  return comp;
}

bool polyADetect(const std::string& seq, int bin, int count, char base) {
  int n = static_cast<int>(seq.size());
  if (bin <= 0 || count <= 0 || n < bin) return false;
  int bc = 0;
  for (int i = 0; i < bin; ++i) bc += (seq[i] == base);
  if (bc >= count) return true;
  for (int i = bin; i < n; ++i) {
    bc += (seq[i] == base);
    bc -= (seq[i - bin] == base);
    if (bc >= count) return true;
  }
  return false;
}

size_t polyARm(std::string seq, const int polyA_len) {
  std::string polyA(polyA_len, 'A');
  return seq.rfind(polyA);
}

int strSlideSearch(std::string seq, const std::string adapter,
                   const int window, const int step, const bool first) {
  size_t pos = first ? seq.find(adapter) : seq.rfind(adapter);
  if (pos != std::string::npos) return static_cast<int>(pos);
  auto sub = strSubset(adapter, window, step);
  std::vector<int> pv;
  for (const auto& s : sub) {
    size_t p = first ? seq.find(s) : seq.rfind(s);
    if (p != std::string::npos) pv.push_back(static_cast<int>(p));
  }
  if (pv.empty() || pv.size() < sub.size() / 2) return -1;
  long long sum = 0;
  for (int p : pv) sum += p;
  return static_cast<int>(sum / pv.size());
}

std::vector<std::string> strSubset(std::string str, const int window,
                                   const int step) {
  std::vector<std::string> vec;
  int len = static_cast<int>(str.size());
  if (window >= len) {
    vec.push_back(str);
  } else {
    int s = 0, e = window;
    while (e < len) {
      vec.push_back(str.substr(s, window));
      s += step; e += step;
    }
    if (e != len) vec.push_back(str.substr(len - window + 1, window));
  }
  return vec;
}

std::string replicate(std::string mode, int times) {
  std::string out;
  for (int i = 0; i < times; i++) out += mode;
  return out;
}

// --- splice_site_correct.h ---

std::vector<std::string> splice_site_cpp(const std::string& input,
                                         const std::string& delimiters) {
  return split_by_chars(input, delimiters);
}

std::map<std::string, int> isoform_count(std::vector<std::string> isoform,
                                         const std::string sep) {
  std::map<std::string, int> counts;
  for (const auto& iso : isoform) {
    std::string mid = get_sub_string(iso, sep, sep);
    counts[mid]++;
  }
  return counts;
}

std::unordered_map<std::string, int> splice_site_count_cpp(
    std::vector<std::string> isoform, const std::string split,
    const std::string sep) {
  auto ic = isoform_count(isoform, sep);
  std::unordered_map<std::string, int> sc;
  for (const auto& p : ic) {
    auto sites = split_by_chars(p.first, split + sep);
    for (const auto& s : sites) sc[s] += p.second;
  }
  return sc;
}

std::string getSubString(const std::string& strValue,
                         const std::string& startChar,
                         const std::string& endChar) {
  return get_sub_string(strValue, startChar, endChar);
}

bool isin(const std::string element, const std::vector<std::string>& vec) {
  for (const auto& v : vec)
    if (v == element) return true;
  return false;
}

// --- exon_corres.h ---

std::pair<int, int> exonstr2bin(std::string exon,
                                const std::string& delimiters) {
  auto sites = split_by_chars(exon, delimiters);
  return {std::stoi(sites[0]), std::stoi(sites[1])};
}

std::vector<std::pair<int, int>> isostr2bins(std::string isoform,
                                             const std::string& delimiters) {
  auto sites = split_by_chars(isoform, delimiters);
  std::vector<std::pair<int, int>> out;
  for (size_t i = 0; i + 1 < sites.size(); i += 2)
    out.emplace_back(std::stoi(sites[i]), std::stoi(sites[i + 1]));
  return out;
}

std::string paste(std::vector<std::string> s, std::string sep) {
  std::string out;
  for (size_t i = 0; i < s.size(); i++) {
    out += s[i];
    if (i != s.size() - 1) out += sep;
  }
  return out;
}

std::vector<std::string> exon_status(std::vector<std::string> exons,
                                     std::string split) {
  if (exons.size() == 1) exons[0] = exons[0] + split + "2";
  for (size_t i = 0; i < exons.size(); i++) {
    if (i == 0) exons[i] = exons[i] + split + "-1";
    else if (i == exons.size() - 1) exons[i] = exons[i] + split + "1";
    else exons[i] = exons[i] + split + "0";
  }
  return exons;
}

// --- umi_dist.h ---

int iso_len(std::vector<int> sites) {
  int len = 0;
  for (size_t i = 0; i + 1 < sites.size(); i += 2)
    len += sites[i + 1] - sites[i] + 1;
  return len;
}

int bin2_intersect(int a_start, int a_end, int b_start, int b_end) {
  if (a_start <= b_end && a_end >= b_start)
    return std::min(a_end, b_end) - std::max(a_start, b_start) + 1;
  return 0;
}

std::vector<int> sites_chop(std::vector<int> sites, int start, int end) {
  std::vector<int> chop;
  for (size_t i = 0; i + 1 < sites.size(); i += 2) {
    int s = sites[i], e = sites[i + 1];
    if (e <= start || s >= end) continue;
    chop.push_back(std::max(s, start));
    chop.push_back(std::min(e, end));
  }
  return chop;
}

int iso2_mid_dist(std::string a, std::string b, const std::string split,
                  const std::string sep) {
  if (a == b) return 0;
  auto sa = longcellsrc::isoform2sites_core(a, split, sep);
  auto sb = longcellsrc::isoform2sites_core(b, split, sep);
  int start = std::max(sa[0], sb[0]);
  int end = std::min(sa[sa.size() - 1], sb[sb.size() - 1]);
  if (start >= end) return -1;
  auto ca = sites_chop(sa, start, end);
  auto cb = sites_chop(sb, start, end);
  int intersect = 0, next_start = 0;
  for (int i = 0; i < static_cast<int>(ca.size()) / 2; i++) {
    int j = next_start;
    bool flag = false;
    while (j < static_cast<int>(cb.size()) / 2) {
      int inter =
          bin2_intersect(ca[2 * i], ca[2 * i + 1], cb[2 * j], cb[2 * j + 1]);
      if (inter > 0) {
        if (!flag) { next_start = j; flag = true; }
        intersect += inter;
      } else if (flag) {
        break;
      }
      j++;
    }
  }
  return iso_len(ca) + iso_len(cb) - 2 * intersect;
}

int minEditDis(std::string seq1, std::string seq2, int k) {
  auto k1 = kmer(seq1, k, 1);
  auto k2 = kmer(seq2, k, 1);
  int best = static_cast<int>(seq1.size());
  for (const auto& s1 : k1) {
    for (const auto& s2 : k2) {
      int d = editDist(s1, s2);
      if (d < best) best = d;
      if (d == 0) break;
    }
  }
  return best;
}

std::vector<std::string> flatten(
    std::vector<std::vector<std::string>> const& vec) {
  std::vector<std::string> out;
  for (const auto& v : vec) out.insert(out.end(), v.begin(), v.end());
  return out;
}

std::vector<std::string> str_split(std::string s, std::string split) {
  std::vector<std::string> out;
  size_t pos;
  while ((pos = s.find(split)) != std::string::npos) {
    out.push_back(s.substr(0, pos));
    s.erase(0, pos + split.length());
  }
  if (!s.empty()) out.push_back(s);
  return out;
}

// --- umi_cluster.h ---

int selectSum(std::map<std::string, int> count, std::vector<std::string> id) {
  int sum = 0;
  for (const auto& i : id) sum += count[i];
  return sum;
}

// --- bc.h ---

double cos_sim(std::vector<int> a, std::vector<int> b) {
  if (a.empty() || b.empty()) return 0;
  std::vector<int> inter;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::back_inserter(inter));
  return inter.size() / (std::sqrt(a.size()) * std::sqrt(b.size()));
}

std::vector<int> kmer_include(std::string seq, std::set<std::string> dic) {
  std::vector<int> sc;
  int id = 0;
  for (const auto& k : dic) {
    if (seq.find(k) != std::string::npos) sc.push_back(id);
    id++;
  }
  return sc;
}

std::vector<std::vector<int>> barcodes_cos_vec(
    std::vector<std::string> barcodes, std::set<std::string> dic) {
  std::vector<std::vector<int>> out;
  for (const auto& b : barcodes) out.push_back(kmer_include(b, dic));
  return out;
}

std::vector<std::string> barcode_cand_cos(
    std::string seq, std::vector<std::string> barcodes,
    std::set<std::string> dic, std::vector<std::vector<int>> index, int top,
    double thresh) {
  auto sc = kmer_include(seq, dic);
  std::vector<std::pair<int, double>> bar_cand;
  for (size_t j = 0; j < barcodes.size(); j++) {
    double cs = cos_sim(sc, index[j]);
    if (cs >= thresh) bar_cand.emplace_back(static_cast<int>(j), cs);
  }
  std::sort(bar_cand.begin(), bar_cand.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });
  std::vector<std::string> out;
  int limit = std::min(top, static_cast<int>(bar_cand.size()));
  for (int i = 0; i < limit; i++) out.push_back(barcodes[bar_cand[i].first]);
  return out;
}

std::vector<std::vector<std::string>> NeighborExtract(
    std::vector<std::string> reads, std::vector<int> id, std::vector<int> pos,
    const int UMI_len, const int flank, const int bar_len) {
  std::vector<std::string> umis, adapters;
  for (size_t i = 0; i < id.size(); i++) {
    int p = pos[i];
    std::string umi, adapter;
    if (p >= UMI_len + flank)
      umi = reads[id[i]].substr(p - UMI_len - flank, UMI_len + 2 * flank);
    if (p + bar_len + UMI_len + 2 * flank <=
        static_cast<int>(reads[id[i]].size()))
      adapter = reads[id[i]].substr(p + bar_len, UMI_len + 2 * flank);
    umis.push_back(umi);
    adapters.push_back(adapter);
  }
  return {umis, adapters};
}

// ============================================================
// longcellsrc namespace — core functions
// ============================================================
namespace longcellsrc {

// qnorm_core, cigar_process_core, ... etc.

// ============================================================
// Normal quantile
// ============================================================

double qnorm_core(double p, double mu, double sigma) {
  if (p <= 0.0) return mu - 10.0 * sigma;
  if (p >= 1.0) return mu + 10.0 * sigma;

  const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                       -2.759285104469687e+02, 1.383577518672690e+02,
                       -3.066479806614716e+01, 2.506628277459239e+00};
  const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                       -1.556989798598866e+02, 6.680131188771972e+01,
                       -1.328068155288572e+01};
  const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                       -2.400758277161838e+00, -2.549732539343734e+00,
                       4.374664141464968e+00,  2.938163982698783e+00};
  const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                       2.445134137142996e+00, 3.754408661907416e+00};

  double q = p - 0.5;
  double r, x;
  if (std::fabs(q) <= 0.425) {
    r = 0.180625 - q * q;
    x = q * (((((a[5] * r + a[4]) * r + a[3]) * r + a[2]) * r + a[1]) * r +
             a[0]) /
        (((((b[4] * r + b[3]) * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
  } else {
    r = (q < 0.0) ? p : 1.0 - p;
    r = std::sqrt(-std::log(r));
    if (r <= 5.0) {
      r -= 1.6;
      x = (((((c[5] * r + c[4]) * r + c[3]) * r + c[2]) * r + c[1]) * r +
           c[0]) /
          ((((d[3] * r + d[2]) * r + d[1]) * r + d[0]) * r + 1.0);
    } else {
      r -= 5.0;
      x = (((((c[5] * r + c[4]) * r + c[3]) * r + c[2]) * r + c[1]) * r +
           c[0]) /
          ((((d[3] * r + d[2]) * r + d[1]) * r + d[0]) * r + 1.0);
    }
    if (q < 0.0) x = -x;
  }
  return mu + sigma * x;
}

// ============================================================
// reads_extraction
// ============================================================

CigarResult cigar_process_core(const std::string& cigar) {
  std::regex pattern(R"(([0-9]+(?:\.[0-9]+)?)\s*([a-zA-Z]+))");
  std::sregex_iterator iter(cigar.begin(), cigar.end(), pattern);
  std::sregex_iterator end;
  CigarResult out;
  for (; iter != end; ++iter) {
    std::smatch match = *iter;
    out.marks.push_back(match[2].str());
    out.counts.push_back(std::stoi(match[1].str()));
  }
  return out;
}

bool cigar_check_core(const std::vector<std::string>& marks,
                      const std::vector<int>& counts) {
  for (size_t i = 0; i < counts.size(); i++)
    if ((marks[i] == "I" || marks[i] == "D") && counts[i] > 6) return false;
  return true;
}

int seq_end_core(int start, const std::vector<std::string>& marks,
                 const std::vector<int>& counts) {
  int end = start;
  for (size_t i = 0; i < marks.size(); i++)
    if (marks[i] == "M" || marks[i] == "D" || marks[i] == "N")
      end += counts[i];
  return end;
}

IsoformMatrix extract_isoform_core(int read_start,
                                   const std::vector<std::string>& marks,
                                   const std::vector<int>& counts,
                                   int refer_start, int refer_end, int flank) {
  IsoformMatrix blocks;
  int start = read_start, end = start;
  bool flag = false;
  for (size_t i = 0; i < marks.size(); i++) {
    if (marks[i] == "M" || marks[i] == "D") {
      end += counts[i];
    } else if (marks[i] == "N") {
      if (start >= refer_start - flank && end <= refer_end + flank)
        blocks.emplace_back(start, end - 1);
      else
        flag = true;
      start = end + counts[i];
      end = start;
    }
  }
  if (start >= refer_start - flank && end <= refer_end + flank)
    blocks.emplace_back(start, end - 1);
  else
    flag = true;
  if (flag && blocks.empty()) blocks.emplace_back(-1, -1);
  return blocks;
}

IsoformMatrix chop_out_bound_core(const IsoformMatrix& isoform,
                                  const IsoformMatrix& annotation,
                                  const std::string& strand, bool& chopped) {
  if (isoform.empty() || isoform[0].first == -1) return isoform;
  int s = 0, e = static_cast<int>(isoform.size()) - 1;
  while (s <= e && isoform[s].second < annotation[0].first) s++;
  while (e >= s && isoform[e].first > annotation.back().second) e--;
  if (strand == "+" && e < static_cast<int>(isoform.size()) - 1) chopped = true;
  if (strand == "-" && s > 0) chopped = true;
  if (s <= e)
    return IsoformMatrix(isoform.begin() + s, isoform.begin() + e + 1);
  return {{-1, -1}};
}

IsoformMatrix splicesite_correct_core(const IsoformMatrix& isoform,
                                      const IsoformMatrix& annotation,
                                      const std::string& strand, int bin) {
  bool chopped = false;
  auto c = chop_out_bound_core(isoform, annotation, strand, chopped);
  if (c.empty() || c[0].first == -1) return c;
  int i_len = static_cast<int>(c.size());
  int r_len = static_cast<int>(annotation.size());
  if (i_len == 1) {
    if (c[0].first < annotation[0].first) c[0].first = annotation[0].first;
    if (c[0].second > annotation[r_len - 1].second)
      c[0].second = annotation[r_len - 1].second;
    return c;
  }
  if (c[0].first < annotation[0].first) c[0].first = annotation[0].first;
  if (c[i_len - 1].second > annotation[r_len - 1].second)
    c[i_len - 1].second = annotation[r_len - 1].second;
  for (int i = 1; i < i_len; i++) {
    int j = 0;
    while (j < r_len && c[i].first > annotation[j].first + bin) j++;
    int best = c[i].first, mindis = bin * 10;
    while (j < r_len &&
           c[i].first >= annotation[j].first - bin &&
           c[i].first <= annotation[j].first + bin) {
      int d = std::abs(c[i].first - annotation[j].first);
      if (d < mindis) { mindis = d; best = annotation[j].first; }
      j++;
    }
    c[i].first = best;
  }
  for (int i = 0; i < i_len - 1; i++) {
    int j = 0;
    while (j < r_len && c[i].second > annotation[j].second + bin) j++;
    int best = c[i].second, mindis = bin * 10;
    while (j < r_len &&
           c[i].second >= annotation[j].second - bin &&
           c[i].second <= annotation[j].second + bin) {
      int d = std::abs(c[i].second - annotation[j].second);
      if (d < mindis) { mindis = d; best = annotation[j].second; }
      j++;
    }
    c[i].second = best;
  }
  return c;
}

int isoform_end_core(const IsoformMatrix& isoform, const std::string& strand) {
  if (strand == "+") return isoform.back().second;
  return isoform[0].first;
}

std::string isoform_to_string_core(const IsoformMatrix& isoform,
                                   const std::string& sep) {
  std::string s;
  for (size_t i = 0; i < isoform.size(); i++) {
    s += std::to_string(isoform[i].first) + "," +
         std::to_string(isoform[i].second);
    if (i != isoform.size() - 1) s += sep;
  }
  return s;
}

ReadsResult extract_reads_core(const std::vector<std::string>& seqs,
                               const std::vector<std::string>& cigars,
                               const std::vector<int>& positions,
                               const IsoformMatrix& annotation,
                               const std::string& strand, int toolkit,
                               int end_flank, int splice_site_bin) {
  ReadsResult out;
  int rs = annotation[0].first;
  int re = annotation.back().second;
  for (size_t i = 0; i < seqs.size(); i++) {
    auto cr = cigar_process_core(cigars[i]);
    if (!cigar_check_core(cr.marks, cr.counts)) continue;
    auto iso = extract_isoform_core(positions[i], cr.marks, cr.counts, rs, re,
                                    end_flank);
    iso = splicesite_correct_core(iso, annotation, strand, splice_site_bin);
    if (iso.empty() || iso[0].first == -1) continue;
    int ie = isoform_end_core(iso, strand);
    std::string is = isoform_to_string_core(iso);
    if (!is.empty()) {
      out.ids.push_back(static_cast<int>(i) + 1);
      out.isoforms.push_back(is);
      out.isoends.push_back(ie);
      out.polyAs.push_back("1");
    }
  }
  return out;
}

// ============================================================
// reads_filter
// ============================================================

std::vector<double> size_filter_core(const std::vector<double>& size,
                                     double ratio) {
  int n = static_cast<int>(size.size());
  double total = 0;
  for (auto s : size) total += s;
  double left = 0, right = total;
  int id = 0;
  double diff = 100;
  for (int i = 0; i < n; i++) {
    left += size[i];
    right -= size[i];
    double t = left - right * ratio / (1.0 - ratio);
    if (std::abs(t) < diff) { diff = std::abs(t); id = i; }
    if (t >= 0) break;
  }
  int eq = 0, le = 0;
  for (int i = 0; i < n; i++) {
    if (size[i] == size[id]) eq++;
    if (size[i] <= size[id]) le++;
  }
  std::vector<double> w(n, 0);
  for (int i = 0; i < n; i++) {
    if (size[i] > size[id]) w[i] = 1;
    else if (size[i] == size[id] && eq > 0)
      w[i] = static_cast<double>(le - id - 1) / eq;
  }
  return w;
}

// ============================================================
// splice_site_correct
// ============================================================

SpliceSiteTableResult splice_site_table_core(
    const std::vector<std::string>& isoforms, const std::string& split,
    const std::string& sep, int thresh) {
  SpliceSiteTableResult out;
  auto ssc = splice_site_count_cpp(isoforms, split, sep);
  std::vector<std::string> ss;
  for (const auto& p : ssc)
    if (p.second >= thresh) ss.push_back(p.first);

  int n_iso = static_cast<int>(isoforms.size());
  int n_ss = static_cast<int>(ss.size());

  if (n_ss > 0) {
    for (int i = 0; i < n_iso; i++) {
      auto sites = split_by_chars(isoforms[i], split + sep);
      int n_sites = static_cast<int>(sites.size());
      std::string start = sites[0], end = sites[n_sites - 1];
      bool bad = false;
      for (int k = 1; k < n_sites - 1; k++) {
        if (std::find(ss.begin(), ss.end(), sites[k]) == ss.end()) {
          bad = true; break;
        }
      }
      if (!bad) {
        out.ids.push_back(i + 1);
        out.starts.push_back(start);
        out.ends.push_back(end);
        std::vector<int> row(n_ss, 0);
        for (int k = 1; k < n_sites - 1; k++) {
          auto it = std::find(ss.begin(), ss.end(), sites[k]);
          if (it != ss.end()) row[it - ss.begin()] = 1;
        }
        out.mid.push_back(row);
      }
    }
  } else {
    for (int i = 0; i < n_iso; i++) {
      auto sites = split_by_chars(isoforms[i], split + sep);
      if (sites.size() == 2) {
        out.ids.push_back(i + 1);
        out.starts.push_back(sites[0]);
        out.ends.push_back(sites[1]);
      }
    }
  }
  out.splice_site_names = ss;
  return out;
}

std::vector<std::vector<bool>> matrix_xor_core(
    const std::vector<std::vector<int>>& mat) {
  int n = static_cast<int>(mat.size());
  if (n == 0) return {};
  int m = static_cast<int>(mat[0].size());
  std::vector<std::vector<bool>> r(n, std::vector<bool>(n, false));
  for (int i = 0; i < n; i++) {
    r[i][i] = true;
    for (int j = i + 1; j < n; j++) {
      bool conflict = false;
      for (int k = 0; k < m; k++) {
        if (mat[i][k] != -1 && mat[j][k] != -1 && mat[i][k] != mat[j][k]) {
          conflict = true; break;
        }
      }
      r[i][j] = r[j][i] = !conflict;
    }
  }
  return r;
}

// ============================================================
// exon_corres
// ============================================================

int bin_sum_core(const std::vector<std::pair<int, int>>& bins) {
  int sum = 0;
  for (const auto& b : bins) {
    if (b.first > b.second)
      throw std::runtime_error(
          "The end position of each bin should be larger than its start "
          "position!");
    sum += b.second - b.first + 1;
  }
  return sum;
}

std::string bin2exonid_core(const std::string& bin_str, int status,
                            const std::vector<int>& start,
                            const std::vector<int>& end,
                            const std::vector<std::string>& exon_id,
                            int mid_bias, int end_bias, int end_overlap,
                            const std::string& nonsense_label,
                            const std::string& split, const std::string& sep) {
  if (status < -1 || status > 2)
    throw std::runtime_error(
        "status must be -1 (start), 0 (middle), 1 (end), or 2 (only one)");
  auto bin = exonstr2bin(bin_str, sep);
  int n_exon = static_cast<int>(exon_id.size());
  std::vector<std::string> eseq;
  std::vector<std::pair<int, int>> left;
  bool flag = false;
  for (int i = 0; i < n_exon; i++) {
    if (start[i] <= bin.second && end[i] >= bin.first) {
      flag = true;
      if (status == 0) {
        eseq.push_back(exon_id[i]);
      } else {
        if (std::min(bin.second, end[i]) - std::max(bin.first, start[i]) + 1 >=
            end_overlap)
          eseq.push_back(exon_id[i]);
      }
      if (start[i] != bin.first) {
        if (status == 0 || status == 1)
          left.emplace_back(std::min(start[i], bin.first),
                            std::max(start[i], bin.first) - 1);
        else if ((status == -1 || status == 2) &&
                 bin.first < start[i] - end_bias)
          left.emplace_back(bin.first, start[i] - 1);
      }
      bin.first = end[i] + 1;
      if (end[i] > bin.second) {
        if (status == 0 || status == -1)
          left.emplace_back(bin.second + 1, end[i]);
        break;
      }
    } else if (flag) {
      break;
    }
  }
  int ll = bin_sum_core(left);
  int el = std::max(bin.second - bin.first + 1, 0);
  if (status == 1 || status == 2) {
    if (el > end_bias) ll += el;
  } else {
    ll += el;
  }
  if (ll > mid_bias) return nonsense_label;
  std::string r;
  for (size_t i = 0; i < eseq.size(); i++) {
    r += eseq[i];
    if (i != eseq.size() - 1) r += split;
  }
  return r;
}

std::map<std::string, std::string> bins2exonids_core(
    const std::vector<std::string>& bins, const std::vector<int>& status,
    const std::vector<int>& start, const std::vector<int>& end,
    const std::vector<std::string>& exon_id, int mid_bias, int end_bias,
    int end_overlap, const std::string& nonsense_label,
    const std::string& split, const std::string& sep) {
  if (bins.size() != status.size())
    throw std::runtime_error(
        "Each bin should have its corresponding status indication!");
  std::map<std::string, std::string> out;
  for (size_t i = 0; i < bins.size(); i++) {
    out[bins[i] + std::to_string(status[i])] = bin2exonid_core(
        bins[i], status[i], start, end, exon_id, mid_bias, end_bias,
        end_overlap, nonsense_label, split, sep);
  }
  return out;
}

static std::vector<std::string> exon_status_helper(
    const std::vector<std::string>& exons, const std::string& split) {
  auto r = exons;
  if (r.size() == 1) r[0] = r[0] + split + "2";
  for (size_t i = 0; i < r.size(); i++) {
    if (i == 0) r[i] = r[i] + split + "-1";
    else if (i == r.size() - 1) r[i] = r[i] + split + "1";
    else r[i] = r[i] + split + "0";
  }
  return r;
}

std::map<std::string, std::string> isos2exonids_index_core(
    const std::vector<std::string>& isoforms, const std::vector<int>& start,
    const std::vector<int>& end, const std::vector<std::string>& exon_id,
    int mid_bias, int end_bias, int end_overlap,
    const std::string& nonsense_label, const std::string& split,
    const std::string& sep) {
  std::set<std::string> eset;
  std::map<std::string, std::vector<std::string>> iso_exonstr;
  for (const auto& iso : isoforms) {
    auto exons = split_by_chars(iso, split);
    auto estatus = exon_status_helper(exons, split);
    iso_exonstr[iso] = estatus;
    for (const auto& es : estatus) eset.insert(es);
  }
  std::vector<std::string> evec(eset.begin(), eset.end());
  std::map<std::string, std::string> e2id;
  for (const auto& es : evec) {
    auto parts = split_by_chars(es, split);
    e2id[es] = bin2exonid_core(parts[0], std::stoi(parts[1]), start, end,
                               exon_id, mid_bias, end_bias, end_overlap,
                               nonsense_label, split, sep);
  }
  std::map<std::string, std::string> iso_exonid;
  for (const auto& iso : isoforms) {
    std::string exonid;
    for (const auto& es : iso_exonstr[iso]) exonid += e2id[es] + split;
    exonid = exonid.substr(0, exonid.size() - 1);
    iso_exonid[iso] = exonid;
  }
  return iso_exonid;
}

// ============================================================
// umi_dist
// ============================================================

std::pair<double, double> iso2_mid_diff_core(const std::string& a,
                                             const std::string& b,
                                             int end_bias,
                                             const std::string& split,
                                             const std::string& sep) {
  if (a == b) return {0.0, 1.0};
  auto sa = isoform2sites_core(a, split, sep);
  auto sb = isoform2sites_core(b, split, sep);
  int a_size = static_cast<int>(sa.size());
  int start = sa.front(), end = sa.back();
  auto cb = sites_chop(sb, start, end);
  int b_sz = static_cast<int>(cb.size());
  if (b_sz == 0) return {-1.0, -1.0};
  if (start <= cb[0] && start >= cb[0] - end_bias && cb[0] <= sa[1])
    start = cb[0];
  if (end >= cb[b_sz - 1] && end <= cb[b_sz - 1] + end_bias &&
      cb[b_sz - 1] >= sa[a_size - 2])
    end = cb[b_sz - 1];
  if (start >= end) return {-1.0, -1.0};
  auto ca = sites_chop(sa, start, end);
  int ca_sz = static_cast<int>(ca.size());
  int cb_sz = static_cast<int>(cb.size());
  int intersect = 0, next_start = 0;
  for (int i = 0; i < ca_sz / 2; i++) {
    int j = next_start;
    bool flag = false;
    while (j < cb_sz / 2) {
      int inter =
          bin2_intersect(ca[2 * i], ca[2 * i + 1], cb[2 * j], cb[2 * j + 1]);
      if (inter > 0) {
        if (!flag) { next_start = j; flag = true; }
        intersect += inter;
      } else if (flag) { break; }
      j++;
    }
  }
  double dis = static_cast<double>(iso_len(ca) + iso_len(cb) - 2 * intersect);
  double ratio = static_cast<double>(intersect) / iso_len(sb);
  return {dis, ratio};
}

int needle_core(const std::string& A, const std::string& B, int match_score,
                int mismatch_score, int gap_score) {
  int n = static_cast<int>(A.size()), m = static_cast<int>(B.size());
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
  for (int i = 0; i <= n; i++) dp[i][0] = i * gap_score;
  for (int j = 0; j <= m; j++) dp[0][j] = j * gap_score;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= m; j++) {
      int S = (A[i - 1] == B[j - 1]) ? match_score : mismatch_score;
      dp[i][j] = std::max({dp[i - 1][j - 1] + S, dp[i - 1][j] + gap_score,
                           dp[i][j - 1] + gap_score});
    }
  return dp[n][m];
}

std::vector<std::vector<int>> index_core(
    const std::vector<std::string>& data,
    const std::vector<std::string>& uniq) {
  std::vector<std::vector<int>> out;
  for (const auto& u : uniq) {
    std::vector<int> idx;
    for (size_t j = 0; j < data.size(); j++)
      if (data[j] == u) idx.push_back(static_cast<int>(j));
    out.push_back(idx);
  }
  return out;
}

std::vector<std::vector<int>> umi_graph_table_core(
    const std::vector<std::string>& umi,
    const std::vector<std::string>& isoform, const std::vector<int>& count,
    int sim_thresh, int iso_thresh, const std::string& split,
    const std::string& sep) {
  int n = static_cast<int>(umi.size());
  if (static_cast<int>(isoform.size()) != n ||
      static_cast<int>(count.size()) != n)
    throw std::runtime_error("The size of isoforms and umi don't match!");

  std::set<std::string> iset(isoform.begin(), isoform.end());
  std::vector<std::string> iuniq(iset.begin(), iset.end());
  std::vector<std::pair<int, int>> ipairs;
  for (size_t i = 0; i < iuniq.size(); i++) {
    for (size_t j = i; j < iuniq.size(); j++) {
      if (i == j) {
        ipairs.emplace_back(static_cast<int>(i), static_cast<int>(j));
      } else {
        int d = iso2_mid_dist(iuniq[i], iuniq[j], split, sep);
        if (d >= 0 && d <= iso_thresh)
          ipairs.emplace_back(static_cast<int>(i), static_cast<int>(j));
      }
    }
  }
  auto iidx = index_core(isoform, iuniq);
  std::set<std::string> uset(umi.begin(), umi.end());
  std::vector<std::string> uuniq(uset.begin(), uset.end());
  int uu_sz = static_cast<int>(uuniq.size());
  std::vector<int> corr(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < uu_sz; j++)
      if (umi[i] == uuniq[j]) { corr[i] = j; break; }

  const int NONE = -100;
  std::vector<std::vector<int>> cache(uu_sz, std::vector<int>(uu_sz, NONE));
  for (int i = 0; i < uu_sz; i++)
    cache[i][i] = static_cast<int>(umi[0].size());

  std::vector<std::vector<int>> result;
  for (const auto& p : ipairs) {
    for (int i : iidx[p.first]) {
      for (int j : iidx[p.second]) {
        int ns = cache[corr[i]][corr[j]];
        if (ns == NONE) {
          ns = needle_core(umi[i], umi[j]);
          cache[corr[i]][corr[j]] = cache[corr[j]][corr[i]] = ns;
        }
        if (ns >= sim_thresh) {
          int score = (i == j) ? count[i] * (count[i] - 1) / 2
                               : count[i] * count[j];
          result.push_back({i, j, ns, score});
        }
      }
    }
  }
  return result;
}

// ============================================================
// umi_cluster
// ============================================================

std::vector<ShareNeighborRow> share_neighbor_core(
    const std::vector<std::string>& index,
    const std::vector<std::vector<std::string>>& neighbor,
    const std::vector<int>& count) {
  int n = static_cast<int>(index.size());
  std::map<std::string, int> mapping;
  for (int i = 0; i < n; i++) mapping[index[i]] = count[i];
  std::map<std::string, std::vector<std::string>> nb;
  for (int i = 0; i < n; i++) {
    auto node = neighbor[i];
    std::sort(node.begin(), node.end());
    nb[index[i]] = node;
  }
  std::vector<ShareNeighborRow> out;
  for (int i = 0; i < n; i++) {
    const auto& n1 = nb[index[i]];
    for (const auto& j : n1) {
      const auto& n2 = nb[j];
      std::vector<std::string> inter;
      std::set_intersection(n1.begin(), n1.end(), n2.begin(), n2.end(),
                            std::back_inserter(inter));
      int sum = 0;
      for (const auto& s : inter) sum += mapping[s];
      out.push_back({index[i], j, sum});
    }
  }
  return out;
}

// ============================================================
// bc
// ============================================================

std::vector<int> pos_filter_core(const std::vector<double>& start,
                                 const std::vector<int>& edit,
                                 std::string& warning_msg) {
  int num = static_cast<int>(start.size());
  std::vector<double> correct;
  std::vector<int> correct_id;
  int max_edit = 0;
  double edit_sum = 0;
  int edit_count[16] = {};
  for (int i = 0; i < num; i++) {
    if (edit[i] == 0) { correct.push_back(start[i]); correct_id.push_back(i); }
    if (edit[i] > max_edit) max_edit = edit[i];
    if (edit[i] < 16) edit_count[edit[i]]++;
    edit_sum += edit[i];
  }
  if (correct.size() < 20) {
    warning_msg = "Too few reads identified with confident cell barcode";
    return correct_id;
  }
  double mu = bc_mean(correct);
  int max_offset = 0;
  int offset_count[100] = {};
  std::vector<int> offsets(num);
  for (int i = 0; i < num; i++) {
    offsets[i] = static_cast<int>(std::round(std::abs(start[i] - mu)));
    if (offsets[i] > max_offset) max_offset = offsets[i];
    if (offsets[i] < 100) offset_count[offsets[i]]++;
  }
  double em = edit_sum / num;
  double erat[16] = {}, orat[100] = {};
  for (int i = max_edit; i >= 0; i--) {
    double c = 0;
    for (int j = i; j <= max_edit; j++) c += edit_count[j];
    erat[i] = c / num;
  }
  for (int i = max_offset; i >= 0; i--) {
    double c = 0;
    for (int j = i; j <= max_offset; j++) c += offset_count[j];
    orat[i] = c / num;
  }
  std::vector<std::pair<int, int>> preserve;
  for (int i = 0; i <= max_edit; i++)
    for (int j = 0; j <= max_offset; j++)
      if (std::pow(erat[i], em) * std::pow(orat[j], 0.5) >= 0.05 || i == 0)
        preserve.emplace_back(i, j);
  std::vector<int> pid;
  for (const auto& p : preserve)
    for (int j = 0; j < num; j++)
      if (edit[j] == p.first && offsets[j] == p.second) pid.push_back(j);
  std::sort(pid.begin(), pid.end());
  return pid;
}

BarcodeMatchResult barcode_match_core(
    const std::vector<std::string>& seq,
    const std::vector<std::string>& barcodes, double mu, double sigma,
    double sigma_start, int k, int batch, int top, double cos_thresh,
    double alpha, int edit_thresh, int UMI_len, int flank) {
  auto dic = kmer(barcodes, k, 1);
  auto idx = barcodes_cos_vec(barcodes, dic);
  int n_seq = static_cast<int>(seq.size());
  int bar_len = static_cast<int>(barcodes[0].length());
  int times = n_seq / batch;

  std::vector<int> rid;
  std::vector<std::string> rbar;
  std::vector<double> rpos;
  std::vector<int> redit;

  for (int t = 0; t <= times; t++) {
    int ss = batch * t, se = batch * (t + 1);
    if (t == times) {
      if (n_seq % batch) se = n_seq;
      else break;
    }
    double iv_s = std::max(0.0, qnorm_core(alpha / 2.0, mu, sigma));
    double iv_e = qnorm_core(1.0 - alpha / 2.0, mu, sigma) + bar_len;
    iv_s = std::round(iv_s);
    iv_e = std::round(iv_e);

    for (int j = ss; j < se; j++) {
      if (static_cast<int>(seq[j].size()) < iv_s + bar_len - 1) continue;
      std::string sub = seq[j].substr(static_cast<size_t>(iv_s),
                                      static_cast<size_t>(iv_e - iv_s + 1));
      if (static_cast<int>(sub.size()) < bar_len) continue;
      auto cand = barcode_cand_cos(sub, barcodes, dic, idx, top, cos_thresh);
      if (cand.empty()) continue;

      int best_edit = bar_len;
      std::pair<int, int> best(best_edit, -100);
      int best_id = 0;
      bool found = false;
      for (size_t p = 0; p < cand.size(); p++) {
        auto temp = minEditDist(sub, cand[p]);
        if (temp.first < best_edit) {
          found = true;
          best_id = static_cast<int>(p);
          best = temp;
          best_edit = temp.first;
        } else if (temp.first == best_edit &&
                   std::abs(temp.second - mu) < std::abs(best.second - mu)) {
          best = temp;
          best_id = static_cast<int>(p);
        }
      }
      if (!found || best.first > edit_thresh) continue;
      rid.push_back(j);
      rbar.push_back(cand[best_id]);
      rpos.push_back(iv_s + best.second);
      redit.push_back(best.first);
    }

    if (rpos.size() > static_cast<size_t>(batch) && rpos.size() < 10000) {
      double s = 0;
      for (auto v : rpos) s += v;
      mu = s / rpos.size();
      double ssq = 0;
      for (auto v : rpos) ssq += (v - mu) * (v - mu);
      int nn = static_cast<int>(rpos.size());
      if (nn > 1)
        sigma = std::sqrt(ssq * (nn + 1) / ((nn + 0.5) * (nn - 1))) +
                sigma_start / nn;
    }
  }

  BarcodeMatchResult out;
  if (!rpos.empty()) {
    std::string wmsg;
    auto pid = pos_filter_core(rpos, redit, wmsg);
    for (int p : pid) {
      out.ids.push_back(rid[p]);
      out.barcodes.push_back(rbar[p]);
      out.positions.push_back(static_cast<int>(rpos[p]));
      out.edits.push_back(redit[p]);
      int pos = out.positions.back();
      int ridx = out.ids.back();
      std::string umi, adapter;
      if (pos >= UMI_len + flank)
        umi = seq[ridx].substr(pos - UMI_len - flank, UMI_len + 2 * flank);
      if (pos + bar_len + UMI_len + 2 * flank <=
          static_cast<int>(seq[ridx].size()))
        adapter = seq[ridx].substr(pos + bar_len, UMI_len + 2 * flank);
      out.umis.push_back(umi);
      out.adapters.push_back(adapter);
    }
  }
  return out;
}

}  // namespace longcellsrc
