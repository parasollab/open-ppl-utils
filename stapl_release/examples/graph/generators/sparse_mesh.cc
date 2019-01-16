#include <stapl/containers/graph/multidigraph.hpp>
#include <stapl/containers/graph/views/graph_view.hpp>
#include <stapl/containers/graph/generators/sparse_mesh.hpp>

stapl::exit_code stapl_main(int argc, char* argv[])
{
  //! [Example]
  using view_type = stapl::graph_view<stapl::multidigraph<int>>;

  auto v = stapl::generators::make_sparse_mesh<view_type>(16, 16, 0.5, true);
  //! [Example]

  std::cout << v.size() << std::endl;

  return 0;
}
