/**
 * @file build.hxx
 * @author Annie Robison (amrobison@ucdavis.edu)
 * @brief
 * @date 2020-12-04
 *
 * @copyright Copyright (c) 2020
 *
 */

#include <gunrock/algorithms/algorithms.hxx>
#include <gunrock/formats/formats.hxx>

namespace gunrock {
namespace graph {
namespace detail {

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t,
          typename new_weight_t>
auto builder(graph::graph_t<
                 space,
                 vertex_t,
                 edge_t,
                 weight_t,
                 graph::graph_csr_t<space, vertex_t, edge_t, weight_t>>& old_G,
             new_weight_t* new_value) {
  // Enable CSR.
  using new_csr_v_t = graph::graph_csr_t<space, vertex_t, edge_t, new_weight_t>;
  using old_csr_v_t = graph::graph_csr_t<space, vertex_t, edge_t, weight_t>;
  using new_graph_type =
      graph::graph_t<space, vertex_t, edge_t, new_weight_t, new_csr_v_t>;

  auto number_of_edges = old_G.template get_number_of_edges<old_csr_v_t>();
  auto number_of_vertices =
      old_G.template get_number_of_vertices<old_csr_v_t>();
  auto v1 = old_G.template getV1<old_csr_v_t>();
  auto v2 = old_G.template getV2<old_csr_v_t>();

  new_graph_type G(old_G.properties);
  G.template set<new_csr_v_t>(number_of_vertices, number_of_edges, v1, v2,
                              new_value);

  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::csr_t<space, vertex_t, edge_t, weight_t>& csr) {
  // Enable CSR.
  using csr_v_t = graph::graph_csr_t<space, vertex_t, edge_t, weight_t>;
  using csr_f_t = format::csr_t<space, vertex_t, edge_t, weight_t>;
  using graph_type = graph::graph_t<space, vertex_t, edge_t, weight_t, csr_v_t>;

  graph_type G(properties);
  G.template set<csr_v_t, csr_f_t>(csr);

  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::coo_t<space, vertex_t, edge_t, weight_t>& coo) {
  // Enable COO.
  using coo_v_t = graph::graph_coo_t<space, vertex_t, edge_t, weight_t>;
  using coo_f_t = format::coo_t<space, vertex_t, edge_t, weight_t>;

  using graph_type = graph::graph_t<space, vertex_t, edge_t, weight_t, coo_v_t>;
  graph_type G(properties);

  G.template set<coo_v_t, coo_f_t>(coo);

  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::csc_t<space, vertex_t, edge_t, weight_t>& csc) {
  // Enable CSC.
  using csc_v_t = graph::graph_csc_t<space, vertex_t, edge_t, weight_t>;
  using csc_f_t = format::csc_t<space, vertex_t, edge_t, weight_t>;

  using graph_type = graph::graph_t<space, vertex_t, edge_t, weight_t, csc_v_t>;
  graph_type G(properties);

  G.template set<csc_v_t, csc_f_t>(csc);
  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::coo_t<space, vertex_t, edge_t, weight_t>& coo,
             format::csr_t<space, vertex_t, edge_t, weight_t>& csr) {
  // Enable CSR.
  using csr_v_t = graph::graph_csr_t<space, vertex_t, edge_t, weight_t>;
  using csr_f_t = format::csr_t<space, vertex_t, edge_t, weight_t>;

  // Enable COO.
  using coo_v_t = graph::graph_coo_t<space, vertex_t, edge_t, weight_t>;
  using coo_f_t = format::coo_t<space, vertex_t, edge_t, weight_t>;

  using graph_type =
      graph::graph_t<space, vertex_t, edge_t, weight_t, csr_v_t, coo_v_t>;

  graph_type G(properties);

  G.template set<csr_v_t, csr_f_t>(csr);
  G.template set<coo_v_t, coo_f_t>(coo);

  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::csc_t<space, vertex_t, edge_t, weight_t>& csc,
             format::csr_t<space, vertex_t, edge_t, weight_t>& csr) {
  // Enable CSR.
  using csr_v_t = graph::graph_csr_t<space, vertex_t, edge_t, weight_t>;
  using csr_f_t = format::csr_t<space, vertex_t, edge_t, weight_t>;

  // Enable CSC.
  using csc_v_t = graph::graph_csc_t<space, vertex_t, edge_t, weight_t>;
  using csc_f_t = format::csc_t<space, vertex_t, edge_t, weight_t>;

  using graph_type =
      graph::graph_t<space, vertex_t, edge_t, weight_t, csr_v_t, csc_v_t>;

  graph_type G(properties);

  G.template set<csr_v_t, csr_f_t>(csr);
  G.template set<csc_v_t, csc_f_t>(csc);

  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::coo_t<space, vertex_t, edge_t, weight_t>& coo,
             format::csc_t<space, vertex_t, edge_t, weight_t>& csc) {
  // Enable COO.
  using coo_v_t = graph::graph_coo_t<space, vertex_t, edge_t, weight_t>;
  using coo_f_t = format::coo_t<space, vertex_t, edge_t, weight_t>;

  // Enable CSC.
  using csc_v_t = graph::graph_csc_t<space, vertex_t, edge_t, weight_t>;
  using csc_f_t = format::csc_t<space, vertex_t, edge_t, weight_t>;

  using graph_type =
      graph::graph_t<space, vertex_t, edge_t, weight_t, csc_v_t, coo_v_t>;

  graph_type G(properties);

  G.template set<coo_v_t, coo_f_t>(coo);
  G.template set<csc_v_t, csc_f_t>(csc);

  return G;
}

template <memory_space_t space,
          typename edge_t,
          typename vertex_t,
          typename weight_t>
auto builder(graph::graph_properties_t properties,
             format::coo_t<space, vertex_t, edge_t, weight_t>& coo,
             format::csc_t<space, vertex_t, edge_t, weight_t>& csc,
             format::csr_t<space, vertex_t, edge_t, weight_t>& csr) {
  // Enable CSR.
  using csr_v_t = graph::graph_csr_t<space, vertex_t, edge_t, weight_t>;
  using csr_f_t = format::csr_t<space, vertex_t, edge_t, weight_t>;

  // Enable COO.
  using coo_v_t = graph::graph_coo_t<space, vertex_t, edge_t, weight_t>;
  using coo_f_t = format::coo_t<space, vertex_t, edge_t, weight_t>;

  // Enable CSC.
  using csc_v_t = graph::graph_csc_t<space, vertex_t, edge_t, weight_t>;
  using csc_f_t = format::csc_t<space, vertex_t, edge_t, weight_t>;

  using graph_type = graph::graph_t<space, vertex_t, edge_t, weight_t, csr_v_t,
                                    csc_v_t, coo_v_t>;

  graph_type G(properties);

  G.template set<csr_v_t, csr_f_t>(csr);
  G.template set<coo_v_t, coo_f_t>(coo);
  G.template set<csc_v_t, csc_f_t>(csc);

  return G;
}
}  // namespace detail
}  // namespace graph
}  // namespace gunrock