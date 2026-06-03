#pragma once

#include <string>
#include <vector>

namespace longcellsrc {

struct IsoDistanceRow {
  int node1;
  int node2;
  int dis;
};

struct IsoSetMidDiffRow {
  int index1;
  int index2;
  double dis;
  double overlap;
};

std::vector<int> isoform2sites_core(const std::string& iso,
                                    const std::string& split = "|",
                                    const std::string& sep = ",");

double iso_length_core(const std::string& iso);

int iso2_dis_core(const std::string& a,
                  const std::string& b,
                  const std::string& split = "|",
                  const std::string& sep = ",");

std::vector<IsoDistanceRow> isos_dis_core(const std::vector<std::string>& isoforms,
                                          int thresh,
                                          const std::string& split = "|",
                                          const std::string& sep = ",");

std::vector<IsoSetMidDiffRow> isoset_mid_diff_core(
    const std::vector<std::string>& iso_set1,
    const std::vector<std::string>& iso_set2,
    int thresh,
    double overlap_thresh,
    int end_bias,
    const std::string& split = "|",
    const std::string& sep = ",");

}  // namespace longcellsrc
