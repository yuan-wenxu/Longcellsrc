#include "isoform_core.h"

#include <algorithm>
#include <stdexcept>

namespace longcellsrc {
namespace {

std::vector<std::string> split_exact(const std::string& value, const std::string& delim) {
  if (delim.empty()) {
    throw std::invalid_argument("Delimiter must not be empty.");
  }

  std::vector<std::string> out;
  size_t start = 0;
  while (true) {
    size_t pos = value.find(delim, start);
    if (pos == std::string::npos) {
      if (start < value.size()) {
        out.push_back(value.substr(start));
      }
      break;
    }
    out.push_back(value.substr(start, pos - start));
    start = pos + delim.size();
  }
  return out;
}

int iso_len_core(const std::vector<int>& sites) {
  if (sites.size() % 2 != 0) {
    throw std::invalid_argument("Isoform site vectors must contain an even number of coordinates.");
  }

  int len = 0;
  for (size_t i = 0; i < sites.size(); i += 2) {
    len += sites[i + 1] - sites[i] + 1;
  }
  return len;
}

int bin2_intersect_core(int a_start, int a_end, int b_start, int b_end) {
  if (a_start <= b_end && a_end >= b_start) {
    return std::min(a_end, b_end) - std::max(a_start, b_start) + 1;
  }
  return 0;
}

std::vector<int> sites_chop_core(const std::vector<int>& sites, int start, int end) {
  if (sites.size() % 2 != 0) {
    throw std::invalid_argument("Isoform site vectors must contain an even number of coordinates.");
  }

  std::vector<int> chop;
  for (size_t i = 0; i < sites.size(); i += 2) {
    const int sites_s = sites[i];
    const int sites_e = sites[i + 1];
    if (sites_e <= start || sites_s >= end) {
      continue;
    }
    chop.push_back(std::max(sites_s, start));
    chop.push_back(std::min(sites_e, end));
  }
  return chop;
}

int iso2_dis_core_impl(const std::string& a,
                       const std::string& b,
                       const std::string& split,
                       const std::string& sep) {
  if (a == b) {
    return 0;
  }

  const std::vector<int> sites_a = isoform2sites_core(a, split, sep);
  const std::vector<int> sites_b = isoform2sites_core(b, split, sep);

  const int a_size = static_cast<int>(sites_a.size());
  const int b_size = static_cast<int>(sites_b.size());

  int i = 0;
  int j = 0;
  int intersect = 0;
  while (i < a_size / 2) {
    if (sites_a[2 * i] > sites_b[2 * j + 1]) {
      if (j == b_size / 2 - 1) {
        break;
      }
      j++;
    } else if (sites_a[2 * i + 1] < sites_b[2 * j]) {
      i++;
    } else {
      intersect += bin2_intersect_core(
          sites_a[2 * i], sites_a[2 * i + 1], sites_b[2 * j], sites_b[2 * j + 1]);
      i++;
    }
  }

  return iso_len_core(sites_a) + iso_len_core(sites_b) - 2 * intersect;
}

std::pair<double, double> iso2_mid_diff_core(const std::string& a,
                                             const std::string& b,
                                             int end_bias,
                                             const std::string& split,
                                             const std::string& sep) {
  if (a == b) {
    return {0.0, 1.0};
  }

  const std::vector<int> sites_a = isoform2sites_core(a, split, sep);
  const std::vector<int> sites_b = isoform2sites_core(b, split, sep);

  const int a_size = static_cast<int>(sites_a.size());
  int start = sites_a.front();
  int end = sites_a.back();

  const std::vector<int> chop_b = sites_chop_core(sites_b, start, end);
  const int b_size = static_cast<int>(chop_b.size());
  if (b_size == 0) {
    return {-1.0, -1.0};
  }

  if (start <= chop_b[0] && start >= chop_b[0] - end_bias && chop_b[0] <= sites_a[1]) {
    start = chop_b[0];
  }
  if (end >= chop_b[b_size - 1] && end <= chop_b[b_size - 1] + end_bias &&
      chop_b[b_size - 1] >= sites_a[a_size - 2]) {
    end = chop_b[b_size - 1];
  }

  if (start >= end) {
    return {-1.0, -1.0};
  }

  const std::vector<int> chop_a = sites_chop_core(sites_a, start, end);
  const int chop_a_size = static_cast<int>(chop_a.size());
  const int chop_b_size = static_cast<int>(chop_b.size());

  int intersect = 0;
  int next_start = 0;
  for (int i = 0; i < chop_a_size / 2; i++) {
    int j = next_start;
    bool flag = false;

    while (j < chop_b_size / 2) {
      const int a_start = chop_a[2 * i];
      const int a_end = chop_a[2 * i + 1];
      const int b_start = chop_b[2 * j];
      const int b_end = chop_b[2 * j + 1];

      const int inter = bin2_intersect_core(a_start, a_end, b_start, b_end);
      if (inter > 0) {
        if (!flag) {
          next_start = j;
          flag = true;
        }
        intersect += inter;
      } else if (flag) {
        break;
      }
      j++;
    }
  }

  const double dis = static_cast<double>(iso_len_core(chop_a) + iso_len_core(chop_b) - 2 * intersect);
  const double ratio = static_cast<double>(intersect) / static_cast<double>(iso_len_core(sites_b));
  return {dis, ratio};
}

void parse_sites_any_delimiter(const std::string& iso, std::vector<double>& values) {
  values.clear();
  double cur = 0.0;
  bool in_num = false;
  bool seen_dot = false;
  double frac = 0.1;

  for (char c : iso) {
    if (c == ' ') {
      continue;
    }

    if (c >= '0' && c <= '9') {
      if (!in_num) {
        in_num = true;
        cur = 0.0;
        seen_dot = false;
        frac = 0.1;
      }
      if (!seen_dot) {
        cur = cur * 10.0 + (c - '0');
      } else {
        cur += (c - '0') * frac;
        frac *= 0.1;
      }
    } else if (c == '.') {
      if (!in_num) {
        in_num = true;
        cur = 0.0;
      }
      if (!seen_dot) {
        seen_dot = true;
        frac = 0.1;
      }
    } else if (in_num) {
      values.push_back(cur);
      in_num = false;
    }
  }

  if (in_num) {
    values.push_back(cur);
  }
}

}  // namespace

std::vector<int> isoform2sites_core(const std::string& iso,
                                    const std::string& split,
                                    const std::string& sep) {
  std::vector<std::string> exons = split_exact(iso, split);
  std::vector<int> sites_int;
  for (const auto& exon : exons) {
    std::vector<std::string> bounds = split_exact(exon, sep);
    for (const auto& bound : bounds) {
      sites_int.push_back(static_cast<int>(std::stod(bound)));
    }
  }
  return sites_int;
}

double iso_length_core(const std::string& iso) {
  std::vector<double> buf;
  parse_sites_any_delimiter(iso, buf);

  double len = 0.0;
  for (size_t i = 0; i + 1 < buf.size(); i += 2) {
    len += buf[i + 1] - buf[i] + 1.0;
  }
  return len;
}

int iso2_dis_core(const std::string& a,
                  const std::string& b,
                  const std::string& split,
                  const std::string& sep) {
  return iso2_dis_core_impl(a, b, split, sep);
}

std::vector<IsoDistanceRow> isos_dis_core(const std::vector<std::string>& isoforms,
                                          int thresh,
                                          const std::string& split,
                                          const std::string& sep) {
  std::vector<IsoDistanceRow> out;
  const int n = static_cast<int>(isoforms.size());
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      const int dis = iso2_dis_core_impl(isoforms[i], isoforms[j], split, sep);
      if (dis <= thresh) {
        out.push_back({i, j, dis});
      }
    }
  }
  return out;
}

std::vector<IsoSetMidDiffRow> isoset_mid_diff_core(
    const std::vector<std::string>& iso_set1,
    const std::vector<std::string>& iso_set2,
    int thresh,
    double overlap_thresh,
    int end_bias,
    const std::string& split,
    const std::string& sep) {
  std::vector<IsoSetMidDiffRow> out;
  const int iso1_size = static_cast<int>(iso_set1.size());
  const int iso2_size = static_cast<int>(iso_set2.size());
  for (int i = 0; i < iso1_size; i++) {
    for (int j = 0; j < iso2_size; j++) {
      const auto [dis, overlap] = iso2_mid_diff_core(iso_set1[i], iso_set2[j], end_bias, split, sep);
      if (dis >= 0 && dis <= thresh && overlap >= overlap_thresh) {
        out.push_back({i, j, dis, overlap});
      }
    }
  }
  return out;
}

}  // namespace longcellsrc
