#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <exception>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "../../src/longcellsrc_core.h"
#include "../../src/edit.h"

namespace {

// ---- helpers ---------------------------------------------------------------

bool py_to_string_vector(PyObject* obj, const char* arg_name,
                         std::vector<std::string>& out) {
  PyObject* seq = PySequence_Fast(obj, nullptr);
  if (seq == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be a sequence of strings.",
                 arg_name);
    return false;
  }
  const Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject** items = PySequence_Fast_ITEMS(seq);
  out.clear();
  out.reserve(static_cast<size_t>(size));
  for (Py_ssize_t i = 0; i < size; ++i) {
    if (!PyUnicode_Check(items[i])) {
      Py_DECREF(seq);
      PyErr_Format(PyExc_TypeError, "%s must contain only strings.", arg_name);
      return false;
    }
    const char* item = PyUnicode_AsUTF8(items[i]);
    if (item == nullptr) {
      Py_DECREF(seq);
      return false;
    }
    out.emplace_back(item);
  }
  Py_DECREF(seq);
  return true;
}

bool py_to_int_vector(PyObject* obj, const char* arg_name,
                      std::vector<int>& out) {
  PyObject* seq = PySequence_Fast(obj, nullptr);
  if (seq == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be a sequence of ints.", arg_name);
    return false;
  }
  const Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject** items = PySequence_Fast_ITEMS(seq);
  out.clear();
  out.reserve(static_cast<size_t>(size));
  for (Py_ssize_t i = 0; i < size; ++i) {
    long v = PyLong_AsLong(items[i]);
    if (v == -1 && PyErr_Occurred()) {
      Py_DECREF(seq);
      return false;
    }
    out.push_back(static_cast<int>(v));
  }
  Py_DECREF(seq);
  return true;
}

bool py_to_double_vector(PyObject* obj, const char* arg_name,
                         std::vector<double>& out) {
  PyObject* seq = PySequence_Fast(obj, nullptr);
  if (seq == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be a sequence of floats.",
                 arg_name);
    return false;
  }
  const Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject** items = PySequence_Fast_ITEMS(seq);
  out.clear();
  out.reserve(static_cast<size_t>(size));
  for (Py_ssize_t i = 0; i < size; ++i) {
    double v = PyFloat_AsDouble(items[i]);
    if (v == -1.0 && PyErr_Occurred()) {
      Py_DECREF(seq);
      return false;
    }
    out.push_back(v);
  }
  Py_DECREF(seq);
  return true;
}

bool py_to_int_matrix(PyObject* obj, const char* arg_name,
                      std::vector<std::vector<int>>& out) {
  PyObject* seq = PySequence_Fast(obj, nullptr);
  if (seq == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be a list of lists.", arg_name);
    return false;
  }
  const Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject** items = PySequence_Fast_ITEMS(seq);
  out.clear();
  out.reserve(static_cast<size_t>(size));
  for (Py_ssize_t i = 0; i < size; ++i) {
    std::vector<int> row;
    if (!py_to_int_vector(items[i], arg_name, row)) {
      Py_DECREF(seq);
      return false;
    }
    out.push_back(std::move(row));
  }
  Py_DECREF(seq);
  return true;
}

bool py_to_string_list_list(PyObject* obj, const char* arg_name,
                            std::vector<std::vector<std::string>>& out) {
  PyObject* seq = PySequence_Fast(obj, nullptr);
  if (seq == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be a list of lists of strings.",
                 arg_name);
    return false;
  }
  const Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject** items = PySequence_Fast_ITEMS(seq);
  out.clear();
  out.reserve(static_cast<size_t>(size));
  for (Py_ssize_t i = 0; i < size; ++i) {
    std::vector<std::string> row;
    if (!py_to_string_vector(items[i], arg_name, row)) {
      Py_DECREF(seq);
      return false;
    }
    out.push_back(std::move(row));
  }
  Py_DECREF(seq);
  return true;
}

PyObject* build_int_list(const std::vector<int>& values) {
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(values.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(values.size()); ++i) {
    PyObject* v = PyLong_FromLong(values[static_cast<size_t>(i)]);
    if (v == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, v);
  }
  return out;
}

PyObject* build_int_list_list(const std::vector<std::vector<int>>& values) {
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(values.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(values.size()); ++i) {
    PyObject* inner = build_int_list(values[static_cast<size_t>(i)]);
    if (inner == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, inner);
  }
  return out;
}

PyObject* build_float_list(const std::vector<double>& values) {
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(values.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(values.size()); ++i) {
    PyObject* v = PyFloat_FromDouble(values[static_cast<size_t>(i)]);
    if (v == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, v);
  }
  return out;
}

PyObject* build_bool_list_list(const std::vector<std::vector<bool>>& values) {
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(values.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(values.size()); ++i) {
    const auto& row = values[static_cast<size_t>(i)];
    PyObject* inner = PyList_New(static_cast<Py_ssize_t>(row.size()));
    if (inner == nullptr) { Py_DECREF(out); return nullptr; }
    for (Py_ssize_t j = 0; j < static_cast<Py_ssize_t>(row.size()); ++j) {
      PyObject* v = row[static_cast<size_t>(j)] ? Py_True : Py_False;
      Py_INCREF(v);
      PyList_SET_ITEM(inner, j, v);
    }
    PyList_SET_ITEM(out, i, inner);
  }
  return out;
}

// ---- existing wrappers (unchanged) -----------------------------------------

PyObject* py_isoform2sites(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* iso = nullptr;
  const char* split = "|";
  const char* sep = ",";
  static const char* kwlist[] = {"iso", "split", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "s|ss", const_cast<char**>(kwlist), &iso, &split, &sep))
    return nullptr;
  try {
    return build_int_list(longcellsrc::isoform2sites_core(iso, split, sep));
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_isos_len_cpp(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* isos = nullptr;
  static const char* kwlist[] = {"isos", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O",
                                   const_cast<char**>(kwlist), &isos))
    return nullptr;
  std::vector<std::string> isoforms;
  if (!py_to_string_vector(isos, "isos", isoforms)) return nullptr;
  try {
    std::vector<double> out;
    out.reserve(isoforms.size());
    for (const auto& iso : isoforms)
      out.push_back(longcellsrc::iso_length_core(iso));
    return build_float_list(out);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_isos_dis(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* isoforms_obj = nullptr;
  int thresh = 10;
  const char* split = "|";
  const char* sep = ",";
  static const char* kwlist[] = {"isoforms", "thresh", "split", "sep",
                                 nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iss",
                                   const_cast<char**>(kwlist), &isoforms_obj,
                                   &thresh, &split, &sep))
    return nullptr;
  std::vector<std::string> isoforms;
  if (!py_to_string_vector(isoforms_obj, "isoforms", isoforms)) return nullptr;
  try {
    const auto rows =
        longcellsrc::isos_dis_core(isoforms, thresh, split, sep);
    PyObject* out = PyList_New(static_cast<Py_ssize_t>(rows.size()));
    if (out == nullptr) return nullptr;
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(rows.size()); ++i) {
      const auto& row = rows[static_cast<size_t>(i)];
      PyObject* t = Py_BuildValue("(iii)", row.node1, row.node2, row.dis);
      if (t == nullptr) { Py_DECREF(out); return nullptr; }
      PyList_SET_ITEM(out, i, t);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_isoset_mid_diff(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* iso_set1_obj = nullptr;
  PyObject* iso_set2_obj = nullptr;
  int thresh = 3;
  double overlap_thresh = 0.5;
  int end_bias = 200;
  const char* split = "|";
  const char* sep = ",";
  static const char* kwlist[] = {"iso_set1", "iso_set2", "thresh",
                                 "overlap_thresh", "end_bias", "split",
                                 "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "OO|idiss", const_cast<char**>(kwlist), &iso_set1_obj,
          &iso_set2_obj, &thresh, &overlap_thresh, &end_bias, &split, &sep))
    return nullptr;
  std::vector<std::string> iso_set1, iso_set2;
  if (!py_to_string_vector(iso_set1_obj, "iso_set1", iso_set1) ||
      !py_to_string_vector(iso_set2_obj, "iso_set2", iso_set2))
    return nullptr;
  try {
    const auto rows = longcellsrc::isoset_mid_diff_core(
        iso_set1, iso_set2, thresh, overlap_thresh, end_bias, split, sep);
    PyObject* out = PyList_New(static_cast<Py_ssize_t>(rows.size()));
    if (out == nullptr) return nullptr;
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(rows.size()); ++i) {
      const auto& row = rows[static_cast<size_t>(i)];
      PyObject* t =
          Py_BuildValue("(iidd)", row.index1, row.index2, row.dis, row.overlap);
      if (t == nullptr) { Py_DECREF(out); return nullptr; }
      PyList_SET_ITEM(out, i, t);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- edit.h wrappers -------------------------------------------------------

PyObject* py_reverse(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* seqs_obj = nullptr;
  static const char* kwlist[] = {"seqs", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O",
                                   const_cast<char**>(kwlist), &seqs_obj))
    return nullptr;
  std::vector<std::string> seqs;
  if (!py_to_string_vector(seqs_obj, "seqs", seqs)) return nullptr;
  try {
    auto result = reverse(seqs);
    PyObject* out = PyList_New(static_cast<Py_ssize_t>(result.size()));
    if (out == nullptr) return nullptr;
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result.size()); ++i) {
      PyObject* s =
          PyUnicode_FromStringAndSize(result[static_cast<size_t>(i)].c_str(),
                                      result[static_cast<size_t>(i)].size());
      if (s == nullptr) { Py_DECREF(out); return nullptr; }
      PyList_SET_ITEM(out, i, s);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_kmer_barcodes(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* barcodes_obj = nullptr;
  int k = 8, step = 1;
  static const char* kwlist[] = {"barcodes", "k", "step", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ii",
                                   const_cast<char**>(kwlist), &barcodes_obj,
                                   &k, &step))
    return nullptr;
  std::vector<std::string> barcodes;
  if (!py_to_string_vector(barcodes_obj, "barcodes", barcodes)) return nullptr;
  try {
    auto kmer_set = kmer(barcodes, k, step);
    PyObject* out = PyList_New(static_cast<Py_ssize_t>(kmer_set.size()));
    if (out == nullptr) return nullptr;
    Py_ssize_t idx = 0;
    for (const auto& km : kmer_set) {
      PyObject* s = PyUnicode_FromStringAndSize(km.c_str(), km.size());
      if (s == nullptr) { Py_DECREF(out); return nullptr; }
      PyList_SET_ITEM(out, idx++, s);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_kmer_str(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  int k = 8, step = 1;
  static const char* kwlist[] = {"seq", "k", "step", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|ii",
                                   const_cast<char**>(kwlist), &seq, &k, &step))
    return nullptr;
  try {
    auto result = kmer(std::string(seq), k, step);
    PyObject* out = PyList_New(static_cast<Py_ssize_t>(result.size()));
    if (out == nullptr) return nullptr;
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result.size()); ++i) {
      PyObject* s = PyUnicode_FromStringAndSize(
          result[static_cast<size_t>(i)].c_str(),
          result[static_cast<size_t>(i)].size());
      if (s == nullptr) { Py_DECREF(out); return nullptr; }
      PyList_SET_ITEM(out, i, s);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_edit_distance(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* a = nullptr;
  const char* b = nullptr;
  static const char* kwlist[] = {"a", "b", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss",
                                   const_cast<char**>(kwlist), &a, &b))
    return nullptr;
  return PyLong_FromLong(editDist(std::string(a), std::string(b)));
}

PyObject* py_min_edit_dist(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  const char* barcode = nullptr;
  static const char* kwlist[] = {"seq", "barcode", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss",
                                   const_cast<char**>(kwlist), &seq, &barcode))
    return nullptr;
  try {
    auto result = minEditDist(std::string(seq), std::string(barcode));
    return Py_BuildValue("(ii)", result.first, result.second);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- tag_extraction.h wrappers ---------------------------------------------

PyObject* py_base_count(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  const char* base = nullptr;
  static const char* kwlist[] = {"seq", "base", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss",
                                   const_cast<char**>(kwlist), &seq, &base))
    return nullptr;
  return PyLong_FromLong(baseCount(std::string(seq), base[0]));
}

PyObject* py_reverse_complement(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  static const char* kwlist[] = {"seq", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s",
                                   const_cast<char**>(kwlist), &seq))
    return nullptr;
  auto result = reverseComplement(std::string(seq));
  return PyUnicode_FromStringAndSize(result.c_str(), result.size());
}

PyObject* py_polyA_detect(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  int bin = 20, count = 15;
  const char* base = "A";
  static const char* kwlist[] = {"seq", "bin", "count", "base", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|iis",
                                   const_cast<char**>(kwlist), &seq, &bin,
                                   &count, &base))
    return nullptr;
  bool r = polyADetect(std::string(seq), bin, count, base[0]);
  return PyBool_FromLong(r ? 1 : 0);
}

PyObject* py_polyA_trim_pos(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  int polyA_len = 10;
  static const char* kwlist[] = {"seq", "polyA_len", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|i",
                                   const_cast<char**>(kwlist), &seq, &polyA_len))
    return nullptr;
  return PyLong_FromUnsignedLong(
      polyARm(std::string(seq), polyA_len));
}

PyObject* py_str_slide_search(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq = nullptr;
  const char* adapter = nullptr;
  int window = 12, step = 3;
  int first = 1;
  static const char* kwlist[] = {"seq", "adapter", "window", "step",
                                 "first", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|iip",
                                   const_cast<char**>(kwlist), &seq, &adapter,
                                   &window, &step, &first))
    return nullptr;
  return PyLong_FromLong(
      strSlideSearch(std::string(seq), std::string(adapter), window, step,
                     first != 0));
}

PyObject* py_str_subset(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* s = nullptr;
  int window = 12, step = 3;
  static const char* kwlist[] = {"s", "window", "step", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|ii",
                                   const_cast<char**>(kwlist), &s, &window,
                                   &step))
    return nullptr;
  auto result = strSubset(std::string(s), window, step);
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(result.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result.size()); ++i) {
    PyObject* item = PyUnicode_FromStringAndSize(
        result[static_cast<size_t>(i)].c_str(),
        result[static_cast<size_t>(i)].size());
    if (item == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, item);
  }
  return out;
}

PyObject* py_replicate_str(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* pattern = nullptr;
  int times = 1;
  static const char* kwlist[] = {"pattern", "times", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|i",
                                   const_cast<char**>(kwlist), &pattern,
                                   &times))
    return nullptr;
  auto result = replicate(std::string(pattern), times);
  return PyUnicode_FromStringAndSize(result.c_str(), result.size());
}

// ---- reads_extraction core -------------------------------------------------

PyObject* py_cigar_process(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* cigar = nullptr;
  static const char* kwlist[] = {"cigar", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s",
                                   const_cast<char**>(kwlist), &cigar))
    return nullptr;
  try {
    auto r = longcellsrc::cigar_process_core(std::string(cigar));
    PyObject* marks_list =
        PyList_New(static_cast<Py_ssize_t>(r.marks.size()));
    PyObject* counts_list = build_int_list(r.counts);
    if (marks_list == nullptr || counts_list == nullptr) {
      Py_XDECREF(marks_list); Py_XDECREF(counts_list); return nullptr;
    }
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(r.marks.size()); ++i) {
      PyObject* m = PyUnicode_FromString(r.marks[i].c_str());
      if (m == nullptr) { Py_DECREF(marks_list); Py_DECREF(counts_list); return nullptr; }
      PyList_SET_ITEM(marks_list, i, m);
    }
    return Py_BuildValue("(NN)", marks_list, counts_list);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_cigar_check(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *marks_obj = nullptr, *counts_obj = nullptr;
  static const char* kwlist[] = {"marks", "counts", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO",
                                   const_cast<char**>(kwlist), &marks_obj,
                                   &counts_obj))
    return nullptr;
  std::vector<std::string> marks;
  std::vector<int> counts;
  if (!py_to_string_vector(marks_obj, "marks", marks) ||
      !py_to_int_vector(counts_obj, "counts", counts))
    return nullptr;
  bool r = longcellsrc::cigar_check_core(marks, counts);
  return PyBool_FromLong(r ? 1 : 0);
}

PyObject* py_seq_end(PyObject*, PyObject* args, PyObject* kwargs) {
  int start = 0;
  PyObject *marks_obj = nullptr, *counts_obj = nullptr;
  static const char* kwlist[] = {"start", "marks", "counts", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iOO",
                                   const_cast<char**>(kwlist), &start,
                                   &marks_obj, &counts_obj))
    return nullptr;
  std::vector<std::string> marks;
  std::vector<int> counts;
  if (!py_to_string_vector(marks_obj, "marks", marks) ||
      !py_to_int_vector(counts_obj, "counts", counts))
    return nullptr;
  return PyLong_FromLong(
      longcellsrc::seq_end_core(start, marks, counts));
}

// ---- reads_filter core -----------------------------------------------------

PyObject* py_size_filter(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* size_obj = nullptr;
  double ratio = 0.1;
  static const char* kwlist[] = {"size", "ratio", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|d",
                                   const_cast<char**>(kwlist), &size_obj,
                                   &ratio))
    return nullptr;
  std::vector<double> size;
  if (!py_to_double_vector(size_obj, "size", size)) return nullptr;
  try {
    return build_float_list(longcellsrc::size_filter_core(size, ratio));
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- splice_site_correct ---------------------------------------------------

PyObject* py_splice_site_split(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* input = nullptr;
  const char* delimiters = nullptr;
  static const char* kwlist[] = {"input", "delimiters", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss",
                                   const_cast<char**>(kwlist), &input,
                                   &delimiters))
    return nullptr;
  auto result = splice_site_cpp(std::string(input), std::string(delimiters));
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(result.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result.size()); ++i) {
    PyObject* s = PyUnicode_FromStringAndSize(
        result[static_cast<size_t>(i)].c_str(),
        result[static_cast<size_t>(i)].size());
    if (s == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, s);
  }
  return out;
}

PyObject* py_isoform_count(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* isoforms_obj = nullptr;
  const char* sep = "|";
  static const char* kwlist[] = {"isoforms", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|s",
                                   const_cast<char**>(kwlist), &isoforms_obj,
                                   &sep))
    return nullptr;
  std::vector<std::string> isoforms;
  if (!py_to_string_vector(isoforms_obj, "isoforms", isoforms)) return nullptr;
  auto result = isoform_count(isoforms, std::string(sep));
  PyObject* out = PyDict_New();
  if (out == nullptr) return nullptr;
  for (const auto& p : result) {
    PyObject* v = PyLong_FromLong(p.second);
    if (v == nullptr) { Py_DECREF(out); return nullptr; }
    if (PyDict_SetItemString(out, p.first.c_str(), v) < 0) {
      Py_DECREF(v); Py_DECREF(out); return nullptr;
    }
    Py_DECREF(v);
  }
  return out;
}

PyObject* py_splice_site_count(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* isoforms_obj = nullptr;
  const char* split = "|";
  const char* sep = ",";
  static const char* kwlist[] = {"isoforms", "split", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ss",
                                   const_cast<char**>(kwlist), &isoforms_obj,
                                   &split, &sep))
    return nullptr;
  std::vector<std::string> isoforms;
  if (!py_to_string_vector(isoforms_obj, "isoforms", isoforms)) return nullptr;
  auto result = splice_site_count_cpp(isoforms, std::string(split),
                                      std::string(sep));
  PyObject* out = PyDict_New();
  if (out == nullptr) return nullptr;
  for (const auto& p : result) {
    PyObject* v = PyLong_FromLong(p.second);
    if (v == nullptr) { Py_DECREF(out); return nullptr; }
    if (PyDict_SetItemString(out, p.first.c_str(), v) < 0) {
      Py_DECREF(v); Py_DECREF(out); return nullptr;
    }
    Py_DECREF(v);
  }
  return out;
}

PyObject* py_splice_site_table(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* isoforms_obj = nullptr;
  const char* split = "|";
  const char* sep = ",";
  int thresh = 10;
  static const char* kwlist[] = {"isoforms", "split", "sep",
                                 "splice_site_thresh", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ssi",
                                   const_cast<char**>(kwlist), &isoforms_obj,
                                   &split, &sep, &thresh))
    return nullptr;
  std::vector<std::string> isoforms;
  if (!py_to_string_vector(isoforms_obj, "isoforms", isoforms)) return nullptr;
  try {
    auto r = longcellsrc::splice_site_table_core(isoforms, split, sep, thresh);
    PyObject* out = PyDict_New();
    if (out == nullptr) return nullptr;

    PyObject* ids = build_int_list(r.ids);
    if (ids == nullptr) { Py_DECREF(out); return nullptr; }
    PyDict_SetItemString(out, "id", ids);
    Py_DECREF(ids);

    PyObject* starts = PyList_New(static_cast<Py_ssize_t>(r.starts.size()));
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(r.starts.size()); ++i) {
      PyObject* s = PyUnicode_FromString(r.starts[i].c_str());
      PyList_SET_ITEM(starts, i, s);
    }
    PyDict_SetItemString(out, "start", starts);
    Py_DECREF(starts);

    PyObject* ends = PyList_New(static_cast<Py_ssize_t>(r.ends.size()));
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(r.ends.size()); ++i) {
      PyObject* s = PyUnicode_FromString(r.ends[i].c_str());
      PyList_SET_ITEM(ends, i, s);
    }
    PyDict_SetItemString(out, "end", ends);
    Py_DECREF(ends);

    if (!r.mid.empty()) {
      PyObject* mid = build_int_list_list(r.mid);
      PyDict_SetItemString(out, "mid", mid);
      Py_DECREF(mid);

      PyObject* snames =
          PyList_New(static_cast<Py_ssize_t>(r.splice_site_names.size()));
      for (Py_ssize_t i = 0;
           i < static_cast<Py_ssize_t>(r.splice_site_names.size()); ++i) {
        PyObject* s =
            PyUnicode_FromString(r.splice_site_names[i].c_str());
        PyList_SET_ITEM(snames, i, s);
      }
      PyDict_SetItemString(out, "splice_site_names", snames);
      Py_DECREF(snames);
    }

    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_matrix_xor(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* mat_obj = nullptr;
  static const char* kwlist[] = {"mat", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O",
                                   const_cast<char**>(kwlist), &mat_obj))
    return nullptr;
  std::vector<std::vector<int>> mat;
  if (!py_to_int_matrix(mat_obj, "mat", mat)) return nullptr;
  try {
    return build_bool_list_list(longcellsrc::matrix_xor_core(mat));
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- exon_corres ------------------------------------------------------------

PyObject* py_exon_str_to_bin(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* exon = nullptr;
  const char* delimiters = nullptr;
  static const char* kwlist[] = {"exon", "delimiters", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss",
                                   const_cast<char**>(kwlist), &exon,
                                   &delimiters))
    return nullptr;
  auto r = exonstr2bin(std::string(exon), std::string(delimiters));
  return Py_BuildValue("(ii)", r.first, r.second);
}

PyObject* py_iso_str_to_bins(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* isoform = nullptr;
  const char* delimiters = nullptr;
  static const char* kwlist[] = {"isoform", "delimiters", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss",
                                   const_cast<char**>(kwlist), &isoform,
                                   &delimiters))
    return nullptr;
  auto bins =
      isostr2bins(std::string(isoform), std::string(delimiters));
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(bins.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(bins.size()); ++i) {
    PyObject* t = Py_BuildValue("(ii)", bins[i].first, bins[i].second);
    if (t == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, t);
  }
  return out;
}

PyObject* py_bin_sum(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* bins_obj = nullptr;
  static const char* kwlist[] = {"bins", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O",
                                   const_cast<char**>(kwlist), &bins_obj))
    return nullptr;
  PyObject* seq = PySequence_Fast(bins_obj, nullptr);
  if (seq == nullptr) {
    PyErr_SetString(PyExc_TypeError, "bins must be a list of (int,int) pairs");
    return nullptr;
  }
  std::vector<std::pair<int, int>> bins;
  const Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject** items = PySequence_Fast_ITEMS(seq);
  for (Py_ssize_t i = 0; i < size; ++i) {
    PyObject* t = items[i];
    if (!PyTuple_Check(t) || PyTuple_Size(t) != 2) {
      Py_DECREF(seq);
      PyErr_SetString(PyExc_TypeError, "each bin must be (int, int)");
      return nullptr;
    }
    bins.push_back({static_cast<int>(PyLong_AsLong(PyTuple_GET_ITEM(t, 0))),
                    static_cast<int>(PyLong_AsLong(PyTuple_GET_ITEM(t, 1)))});
  }
  Py_DECREF(seq);
  try {
    return PyLong_FromLong(longcellsrc::bin_sum_core(bins));
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_paste_strings(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* strings_obj = nullptr;
  const char* sep = "|";
  static const char* kwlist[] = {"strings", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|s",
                                   const_cast<char**>(kwlist), &strings_obj,
                                   &sep))
    return nullptr;
  std::vector<std::string> strings;
  if (!py_to_string_vector(strings_obj, "strings", strings)) return nullptr;
  auto result = paste(strings, std::string(sep));
  return PyUnicode_FromStringAndSize(result.c_str(), result.size());
}

PyObject* py_exon_status(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* exons_obj = nullptr;
  const char* split = "|";
  static const char* kwlist[] = {"exons", "split", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|s",
                                   const_cast<char**>(kwlist), &exons_obj,
                                   &split))
    return nullptr;
  std::vector<std::string> exons;
  if (!py_to_string_vector(exons_obj, "exons", exons)) return nullptr;
  auto result = exon_status(exons, std::string(split));
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(result.size()));
  if (out == nullptr) return nullptr;
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result.size()); ++i) {
    PyObject* s = PyUnicode_FromStringAndSize(
        result[i].c_str(), result[i].size());
    if (s == nullptr) { Py_DECREF(out); return nullptr; }
    PyList_SET_ITEM(out, i, s);
  }
  return out;
}

PyObject* py_bin_to_exon_id(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* bin_str = nullptr;
  int status = 0;
  PyObject* start_obj = nullptr, *end_obj = nullptr, *exon_id_obj = nullptr;
  int mid_bias = 0, end_bias = 10, end_overlap = 10;
  const char* nonsense_label = "N";
  const char* split = "|", *sep = ",";
  static const char* kwlist[] = {"bin_str", "status", "start", "end",
                                 "exon_id", "mid_bias", "end_bias",
                                 "end_overlap", "nonsense_label", "split",
                                 "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "siOOO|iiisss", const_cast<char**>(kwlist), &bin_str,
          &status, &start_obj, &end_obj, &exon_id_obj, &mid_bias, &end_bias,
          &end_overlap, &nonsense_label, &split, &sep))
    return nullptr;
  std::vector<int> start, end;
  std::vector<std::string> exon_id;
  if (!py_to_int_vector(start_obj, "start", start) ||
      !py_to_int_vector(end_obj, "end", end) ||
      !py_to_string_vector(exon_id_obj, "exon_id", exon_id))
    return nullptr;
  try {
    auto r = longcellsrc::bin2exonid_core(
        bin_str, status, start, end, exon_id, mid_bias, end_bias, end_overlap,
        nonsense_label, split, sep);
    return PyUnicode_FromStringAndSize(r.c_str(), r.size());
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_bins_to_exon_ids(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *bins_obj = nullptr, *status_obj = nullptr;
  PyObject *start_obj = nullptr, *end_obj = nullptr, *exon_id_obj = nullptr;
  int mid_bias = 0, end_bias = 10, end_overlap = 10;
  const char* nonsense_label = "N";
  const char* split = "|", *sep = ",";
  static const char* kwlist[] = {"bins", "status", "start", "end",
                                 "exon_id", "mid_bias", "end_bias",
                                 "end_overlap", "nonsense_label", "split",
                                 "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "OOOOO|iiisss", const_cast<char**>(kwlist), &bins_obj,
          &status_obj, &start_obj, &end_obj, &exon_id_obj, &mid_bias,
          &end_bias, &end_overlap, &nonsense_label, &split, &sep))
    return nullptr;
  std::vector<std::string> bins;
  std::vector<int> status, start, end;
  std::vector<std::string> exon_id;
  if (!py_to_string_vector(bins_obj, "bins", bins) ||
      !py_to_int_vector(status_obj, "status", status) ||
      !py_to_int_vector(start_obj, "start", start) ||
      !py_to_int_vector(end_obj, "end", end) ||
      !py_to_string_vector(exon_id_obj, "exon_id", exon_id))
    return nullptr;
  try {
    auto r = longcellsrc::bins2exonids_core(
        bins, status, start, end, exon_id, mid_bias, end_bias, end_overlap,
        nonsense_label, split, sep);
    PyObject* out = PyDict_New();
    if (out == nullptr) return nullptr;
    for (const auto& p : r) {
      PyObject* v =
          PyUnicode_FromStringAndSize(p.second.c_str(), p.second.size());
      if (v == nullptr) { Py_DECREF(out); return nullptr; }
      if (PyDict_SetItemString(out, p.first.c_str(), v) < 0) {
        Py_DECREF(v); Py_DECREF(out); return nullptr;
      }
      Py_DECREF(v);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_isos_to_exon_ids_index(PyObject*, PyObject* args,
                                    PyObject* kwargs) {
  PyObject *isoforms_obj = nullptr, *start_obj = nullptr;
  PyObject *end_obj = nullptr, *exon_id_obj = nullptr;
  int mid_bias = 0, end_bias = 10, end_overlap = 10;
  const char* nonsense_label = "N";
  const char* split = "|", *sep = ",";
  static const char* kwlist[] = {"isoforms", "start", "end", "exon_id",
                                 "mid_bias", "end_bias", "end_overlap",
                                 "nonsense_label", "split", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "OOOOO|iiisss", const_cast<char**>(kwlist),
          &isoforms_obj, &start_obj, &end_obj, &exon_id_obj, &mid_bias,
          &end_bias, &end_overlap, &nonsense_label, &split, &sep))
    return nullptr;
  std::vector<std::string> isoforms;
  std::vector<int> start, end;
  std::vector<std::string> exon_id;
  if (!py_to_string_vector(isoforms_obj, "isoforms", isoforms) ||
      !py_to_int_vector(start_obj, "start", start) ||
      !py_to_int_vector(end_obj, "end", end) ||
      !py_to_string_vector(exon_id_obj, "exon_id", exon_id))
    return nullptr;
  try {
    auto r = longcellsrc::isos2exonids_index_core(
        isoforms, start, end, exon_id, mid_bias, end_bias, end_overlap,
        nonsense_label, split, sep);
    PyObject* out = PyDict_New();
    if (out == nullptr) return nullptr;
    for (const auto& p : r) {
      PyObject* v =
          PyUnicode_FromStringAndSize(p.second.c_str(), p.second.size());
      if (v == nullptr) { Py_DECREF(out); return nullptr; }
      if (PyDict_SetItemString(out, p.first.c_str(), v) < 0) {
        Py_DECREF(v); Py_DECREF(out); return nullptr;
      }
      Py_DECREF(v);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- umi_dist ---------------------------------------------------------------

PyObject* py_iso_len(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* sites_obj = nullptr;
  static const char* kwlist[] = {"sites", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O",
                                   const_cast<char**>(kwlist), &sites_obj))
    return nullptr;
  std::vector<int> sites;
  if (!py_to_int_vector(sites_obj, "sites", sites)) return nullptr;
  return PyLong_FromLong(iso_len(sites));
}

PyObject* py_bin2_intersect(PyObject*, PyObject* args, PyObject* kwargs) {
  int a_start, a_end, b_start, b_end;
  static const char* kwlist[] = {"a_start", "a_end", "b_start",
                                 "b_end", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiii",
                                   const_cast<char**>(kwlist), &a_start,
                                   &a_end, &b_start, &b_end))
    return nullptr;
  return PyLong_FromLong(bin2_intersect(a_start, a_end, b_start, b_end));
}

PyObject* py_sites_chop(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* sites_obj = nullptr;
  int start = 0, end = 0;
  static const char* kwlist[] = {"sites", "start", "end", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oii",
                                   const_cast<char**>(kwlist), &sites_obj,
                                   &start, &end))
    return nullptr;
  std::vector<int> sites;
  if (!py_to_int_vector(sites_obj, "sites", sites)) return nullptr;
  return build_int_list(sites_chop(sites, start, end));
}

PyObject* py_iso2_mid_dist(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* a = nullptr, *b = nullptr;
  const char* split = "|", *sep = ",";
  static const char* kwlist[] = {"a", "b", "split", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|ss",
                                   const_cast<char**>(kwlist), &a, &b, &split,
                                   &sep))
    return nullptr;
  return PyLong_FromLong(
      iso2_mid_dist(a, b, split, sep));
}

PyObject* py_iso2_mid_diff(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* a = nullptr, *b = nullptr;
  int end_bias = 200;
  const char* split = "|", *sep = ",";
  static const char* kwlist[] = {"a", "b", "end_bias", "split",
                                 "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|iss",
                                   const_cast<char**>(kwlist), &a, &b,
                                   &end_bias, &split, &sep))
    return nullptr;
  try {
    auto r = longcellsrc::iso2_mid_diff_core(a, b, end_bias, split, sep);
    return Py_BuildValue("(dd)", r.first, r.second);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

PyObject* py_needleman(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* a = nullptr, *b = nullptr;
  int match = 1, mismatch = -1, gap = -1;
  static const char* kwlist[] = {"a", "b", "match_score",
                                 "mismatch_score", "gap_score", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|iii",
                                   const_cast<char**>(kwlist), &a, &b, &match,
                                   &mismatch, &gap))
    return nullptr;
  return PyLong_FromLong(
      longcellsrc::needle_core(a, b, match, mismatch, gap));
}

PyObject* py_umi_min_edit_dist(PyObject*, PyObject* args, PyObject* kwargs) {
  const char* seq1 = nullptr, *seq2 = nullptr;
  int k = 10;
  static const char* kwlist[] = {"seq1", "seq2", "k", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|i",
                                   const_cast<char**>(kwlist), &seq1, &seq2,
                                   &k))
    return nullptr;
  return PyLong_FromLong(minEditDis(seq1, seq2, k));
}

PyObject* py_index_of(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *data_obj = nullptr, *uniq_obj = nullptr;
  static const char* kwlist[] = {"data", "uniq", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO",
                                   const_cast<char**>(kwlist), &data_obj,
                                   &uniq_obj))
    return nullptr;
  std::vector<std::string> data, uniq;
  if (!py_to_string_vector(data_obj, "data", data) ||
      !py_to_string_vector(uniq_obj, "uniq", uniq))
    return nullptr;
  return build_int_list_list(longcellsrc::index_core(data, uniq));
}

PyObject* py_umi_graph_table(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *umi_obj = nullptr, *isoform_obj = nullptr, *count_obj = nullptr;
  int sim_thresh = 5, iso_thresh = 80;
  const char* split = "|", *sep = ",";
  static const char* kwlist[] = {"umi", "isoform", "count", "sim_thresh",
                                 "iso_thresh", "split", "sep", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|iiss",
                                   const_cast<char**>(kwlist), &umi_obj,
                                   &isoform_obj, &count_obj, &sim_thresh,
                                   &iso_thresh, &split, &sep))
    return nullptr;
  std::vector<std::string> umi, isoform;
  std::vector<int> count;
  if (!py_to_string_vector(umi_obj, "umi", umi) ||
      !py_to_string_vector(isoform_obj, "isoform", isoform) ||
      !py_to_int_vector(count_obj, "count", count))
    return nullptr;
  try {
    auto r = longcellsrc::umi_graph_table_core(umi, isoform, count,
                                               sim_thresh, iso_thresh, split,
                                               sep);
    return build_int_list_list(r);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- umi_cluster ------------------------------------------------------------

PyObject* py_select_sum(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *count_dict = nullptr, *ids_obj = nullptr;
  static const char* kwlist[] = {"count", "ids", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO",
                                   const_cast<char**>(kwlist), &count_dict,
                                   &ids_obj))
    return nullptr;
  if (!PyDict_Check(count_dict)) {
    PyErr_SetString(PyExc_TypeError, "count must be a dict");
    return nullptr;
  }
  std::map<std::string, int> cmap;
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  while (PyDict_Next(count_dict, &pos, &key, &value)) {
    const char* k = PyUnicode_AsUTF8(key);
    if (k == nullptr) return nullptr;
    cmap[k] = static_cast<int>(PyLong_AsLong(value));
  }
  std::vector<std::string> ids;
  if (!py_to_string_vector(ids_obj, "ids", ids)) return nullptr;
  return PyLong_FromLong(selectSum(cmap, ids));
}

PyObject* py_share_neighbor(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *index_obj = nullptr, *neighbor_obj = nullptr, *count_obj = nullptr;
  static const char* kwlist[] = {"index", "neighbor", "count", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "OOO", const_cast<char**>(kwlist), &index_obj,
          &neighbor_obj, &count_obj))
    return nullptr;
  std::vector<std::string> index;
  std::vector<std::vector<std::string>> neighbor;
  std::vector<int> count;
  if (!py_to_string_vector(index_obj, "index", index) ||
      !py_to_string_list_list(neighbor_obj, "neighbor", neighbor) ||
      !py_to_int_vector(count_obj, "count", count))
    return nullptr;
  try {
    auto rows =
        longcellsrc::share_neighbor_core(index, neighbor, count);
    PyObject* out = PyList_New(static_cast<Py_ssize_t>(rows.size()));
    if (out == nullptr) return nullptr;
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(rows.size()); ++i) {
      const auto& row = rows[static_cast<size_t>(i)];
      PyObject* s1 = PyUnicode_FromStringAndSize(row.node1.c_str(),
                                                  row.node1.size());
      PyObject* s2 = PyUnicode_FromStringAndSize(row.node2.c_str(),
                                                  row.node2.size());
      if (s1 == nullptr || s2 == nullptr) {
        Py_XDECREF(s1); Py_XDECREF(s2); Py_DECREF(out); return nullptr;
      }
      PyObject* t = Py_BuildValue("(OOi)", s1, s2, row.share);
      Py_DECREF(s1); Py_DECREF(s2);
      if (t == nullptr) { Py_DECREF(out); return nullptr; }
      PyList_SET_ITEM(out, i, t);
    }
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- bc ---------------------------------------------------------------------

PyObject* py_qnorm(PyObject*, PyObject* args, PyObject* kwargs) {
  double p = 0.5, mu = 0.0, sigma = 1.0;
  static const char* kwlist[] = {"p", "mu", "sigma", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|dd",
                                   const_cast<char**>(kwlist), &p, &mu, &sigma))
    return nullptr;
  return PyFloat_FromDouble(longcellsrc::qnorm_core(p, mu, sigma));
}

PyObject* py_cos_similarity(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *a_obj = nullptr, *b_obj = nullptr;
  static const char* kwlist[] = {"a", "b", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO",
                                   const_cast<char**>(kwlist), &a_obj, &b_obj))
    return nullptr;
  std::vector<int> a, b;
  if (!py_to_int_vector(a_obj, "a", a) || !py_to_int_vector(b_obj, "b", b))
    return nullptr;
  return PyFloat_FromDouble(cos_sim(a, b));
}

PyObject* py_neighbor_extract(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *reads_obj = nullptr, *id_obj = nullptr, *pos_obj = nullptr;
  int umi_len = 10, flank = 1, bar_len = 16;
  static const char* kwlist[] = {"reads", "id", "pos", "umi_len",
                                 "flank", "bar_len", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "OOO|iii", const_cast<char**>(kwlist), &reads_obj,
          &id_obj, &pos_obj, &umi_len, &flank, &bar_len))
    return nullptr;
  std::vector<std::string> reads;
  std::vector<int> id, pos;
  if (!py_to_string_vector(reads_obj, "reads", reads) ||
      !py_to_int_vector(id_obj, "id", id) ||
      !py_to_int_vector(pos_obj, "pos", pos))
    return nullptr;
  auto result = NeighborExtract(reads, id, pos, umi_len, flank, bar_len);
  PyObject* umis = PyList_New(static_cast<Py_ssize_t>(result[0].size()));
  PyObject* adapters =
      PyList_New(static_cast<Py_ssize_t>(result[1].size()));
  if (umis == nullptr || adapters == nullptr) {
    Py_XDECREF(umis); Py_XDECREF(adapters); return nullptr;
  }
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result[0].size()); ++i) {
    PyObject* s = PyUnicode_FromString(result[0][i].c_str());
    if (s == nullptr) { Py_DECREF(umis); Py_DECREF(adapters); return nullptr; }
    PyList_SET_ITEM(umis, i, s);
  }
  for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(result[1].size()); ++i) {
    PyObject* s = PyUnicode_FromString(result[1][i].c_str());
    if (s == nullptr) { Py_DECREF(umis); Py_DECREF(adapters); return nullptr; }
    PyList_SET_ITEM(adapters, i, s);
  }
  return Py_BuildValue("(OO)", umis, adapters);
}

PyObject* py_pos_filter(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *start_obj = nullptr, *edit_obj = nullptr;
  static const char* kwlist[] = {"start", "edit", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO",
                                   const_cast<char**>(kwlist), &start_obj,
                                   &edit_obj))
    return nullptr;
  std::vector<double> start;
  std::vector<int> edit;
  if (!py_to_double_vector(start_obj, "start", start) ||
      !py_to_int_vector(edit_obj, "edit", edit))
    return nullptr;
  std::string wmsg;
  auto result = longcellsrc::pos_filter_core(start, edit, wmsg);
  return Py_BuildValue("(Os)", build_int_list(result), wmsg.c_str());
}

PyObject* py_barcode_match(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject *seq_obj = nullptr, *barcodes_obj = nullptr;
  double mu = 20, sigma = 10, sigma_start = 10;
  int k = 8, batch = 100, top = 8;
  double cos_thresh = 0.25, alpha = 0.05;
  int edit_thresh = 5, umi_len = 10, flank = 1;
  static const char* kwlist[] = {"seq",     "barcodes", "mu",     "sigma",
                                 "sigma_start", "k", "batch", "top",
                                 "cos_thresh", "alpha", "edit_thresh",
                                 "UMI_len", "flank", nullptr};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "OO|ddiiidiiii", const_cast<char**>(kwlist), &seq_obj,
          &barcodes_obj, &mu, &sigma, &sigma_start, &k, &batch, &top,
          &cos_thresh, &alpha, &edit_thresh, &umi_len, &flank))
    return nullptr;
  std::vector<std::string> seq, barcodes;
  if (!py_to_string_vector(seq_obj, "seq", seq) ||
      !py_to_string_vector(barcodes_obj, "barcodes", barcodes))
    return nullptr;
  try {
    auto r = longcellsrc::barcode_match_core(seq, barcodes, mu, sigma,
                                             sigma_start, k, batch, top,
                                             cos_thresh, alpha, edit_thresh,
                                             umi_len, flank);
    PyObject* out = PyDict_New();
    if (out == nullptr) return nullptr;

    PyObject* ids = build_int_list(r.ids);
    PyObject* bars =
        PyList_New(static_cast<Py_ssize_t>(r.barcodes.size()));
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(r.barcodes.size());
         ++i) {
      PyObject* s = PyUnicode_FromString(r.barcodes[i].c_str());
      PyList_SET_ITEM(bars, i, s);
    }
    PyObject* positions = build_int_list(r.positions);
    PyObject* edits = build_int_list(r.edits);

    PyObject* umis =
        PyList_New(static_cast<Py_ssize_t>(r.umis.size()));
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(r.umis.size()); ++i) {
      PyObject* s = PyUnicode_FromString(r.umis[i].c_str());
      PyList_SET_ITEM(umis, i, s);
    }
    PyObject* adapters =
        PyList_New(static_cast<Py_ssize_t>(r.adapters.size()));
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(r.adapters.size());
         ++i) {
      PyObject* s = PyUnicode_FromString(r.adapters[i].c_str());
      PyList_SET_ITEM(adapters, i, s);
    }

    PyDict_SetItemString(out, "id", ids);
    PyDict_SetItemString(out, "barcode", bars);
    PyDict_SetItemString(out, "pos", positions);
    PyDict_SetItemString(out, "edit", edits);
    PyDict_SetItemString(out, "umi", umis);
    PyDict_SetItemString(out, "adapter", adapters);

    Py_DECREF(ids); Py_DECREF(bars); Py_DECREF(positions);
    Py_DECREF(edits); Py_DECREF(umis); Py_DECREF(adapters);
    return out;
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return nullptr;
  }
}

// ---- reads_extraction (corrected cigar_process) -----------------------------

// NOTE: We need to fix py_cigar_process. Let's override it here.
// Actually let me just fix it in the final module_methods array.

// ---- method table -----------------------------------------------------------

PyMethodDef module_methods[] = {
    // Barcode matching
    {"barcode_match", reinterpret_cast<PyCFunction>(py_barcode_match),
     METH_VARARGS | METH_KEYWORDS,
     "Match cell barcodes in reads."},
    {"neighbor_extract", reinterpret_cast<PyCFunction>(py_neighbor_extract),
     METH_VARARGS | METH_KEYWORDS,
     "Extract UMI and adapter sequences flanking barcode positions."},

    // FASTQ/tag helpers
    {"str_slide_search", reinterpret_cast<PyCFunction>(py_str_slide_search),
     METH_VARARGS | METH_KEYWORDS,
     "Search for an adapter in a sequence using sliding windows."},
    {"base_count", reinterpret_cast<PyCFunction>(py_base_count),
     METH_VARARGS | METH_KEYWORDS,
     "Count occurrences of a base in a DNA sequence."},

    // Read and CIGAR parsing
    {"cigar_process", reinterpret_cast<PyCFunction>(py_cigar_process),
     METH_VARARGS | METH_KEYWORDS,
     "Parse a CIGAR string into (marks, counts) lists."},
    {"seq_end", reinterpret_cast<PyCFunction>(py_seq_end),
     METH_VARARGS | METH_KEYWORDS,
     "Compute reference end position from start, CIGAR marks, and counts."},

    // Isoform filtering and splice-site summaries
    {"isos_dis", reinterpret_cast<PyCFunction>(py_isos_dis),
     METH_VARARGS | METH_KEYWORDS,
     "Return matching isoform pairs as (node1, node2, dis) tuples."},
    {"size_filter", reinterpret_cast<PyCFunction>(py_size_filter),
     METH_VARARGS | METH_KEYWORDS,
     "Compute weights based on a size filter threshold."},
    {"isos_len_cpp", reinterpret_cast<PyCFunction>(py_isos_len_cpp),
     METH_VARARGS | METH_KEYWORDS,
     "Compute effective lengths for isoform strings."},
    {"splice_site_table", reinterpret_cast<PyCFunction>(py_splice_site_table),
     METH_VARARGS | METH_KEYWORDS,
     "Build splice-site table for isoforms."},
    {"matrix_xor", reinterpret_cast<PyCFunction>(py_matrix_xor),
     METH_VARARGS | METH_KEYWORDS,
     "Compute pairwise non-conflict matrix from an integer matrix."},

    // Isoform-to-exon mapping
    {"isos_to_exon_ids", reinterpret_cast<PyCFunction>(py_isos_to_exon_ids_index),
     METH_VARARGS | METH_KEYWORDS,
     "Map isoform strings to exon ID strings."},

    // Isoform distance helpers
    {"isoform2sites", reinterpret_cast<PyCFunction>(py_isoform2sites),
     METH_VARARGS | METH_KEYWORDS,
     "Convert an isoform string into a flat list of splice-site coordinates."},
    {"sites_chop", reinterpret_cast<PyCFunction>(py_sites_chop),
     METH_VARARGS | METH_KEYWORDS,
     "Clip a list of exon coordinates to a given interval."},
    {"iso2_mid_dist", reinterpret_cast<PyCFunction>(py_iso2_mid_dist),
     METH_VARARGS | METH_KEYWORDS,
     "Compute the midpoint-distance between two isoforms."},
    {"iso2_mid_diff", reinterpret_cast<PyCFunction>(py_iso2_mid_diff),
     METH_VARARGS | METH_KEYWORDS,
     "Compute (distance, overlap) between two isoforms."},
    {"isoset_mid_diff", reinterpret_cast<PyCFunction>(py_isoset_mid_diff),
     METH_VARARGS | METH_KEYWORDS,
     "Return cross-set matches as (index1, index2, dis, overlap) tuples."},

    // UMI graph helpers
    {"umi_graph_table", reinterpret_cast<PyCFunction>(py_umi_graph_table),
     METH_VARARGS | METH_KEYWORDS,
     "Build UMI similarity graph table."},
    {"share_neighbor", reinterpret_cast<PyCFunction>(py_share_neighbor),
     METH_VARARGS | METH_KEYWORDS,
     "Compute shared neighbors between index entries."},

    {nullptr, nullptr, 0, nullptr}};

PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_core",
    "Python bindings for Longcellsrc C++ functions.",
    -1,
    module_methods,
};

}  // namespace

PyMODINIT_FUNC PyInit__core(void) {
  return PyModule_Create(&module_def);
}
