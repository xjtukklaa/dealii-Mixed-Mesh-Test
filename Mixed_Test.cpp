#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
 
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
 
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
 
#include <fstream>
#include <iostream>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
 
using namespace dealii;
const int dim = 3;
const int degree = 1;
class Mixed_Mesh
{
public:
  Mixed_Mesh();
 
  void run();
 
private:
  void make_grid();
  void set_up_system();
  void assemble_system();
  void solve();
  void output_results();

  Triangulation<dim> triangulation;
  
  const hp::MappingCollection<dim> mapping;
  const hp::FECollection<dim>      fe;
  const hp::QCollection<dim>       quadrature_formula;
  
  DoFHandler<dim> dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;

  double E;
  double possion;
  double rho;
};
 
Mixed_Mesh::Mixed_Mesh()
  // correct oder 
  // : mapping(MappingFE<dim>(FE_PyramidP<dim>(degree)), MappingFE<dim>(FE_SimplexP<dim>(degree)), MappingFE<dim>(FE_Q<dim>(degree)))
  // , fe(FESystem<dim>(FE_PyramidP<dim>(degree)^3), FESystem<dim>(FE_SimplexP<dim>(degree)^3), FESystem<dim>(FE_Q<dim>(degree)^3))
  // , quadrature_formula(QGaussPyramid<dim>(degree + 1), QGaussSimplex<dim>(degree + 1), QGauss<dim>(degree + 1))
  // error oder
  : mapping(MappingFE<dim>(FE_SimplexP<dim>(degree)), MappingFE<dim>(FE_PyramidP<dim>(degree)), MappingFE<dim>(FE_Q<dim>(degree)))
  , fe(FESystem<dim>(FE_SimplexP<dim>(degree)^3), FESystem<dim>(FE_PyramidP<dim>(degree)^3), FESystem<dim>(FE_Q<dim>(degree)^3))
  , quadrature_formula(QGaussSimplex<dim>(degree + 1), QGaussPyramid<dim>(degree + 1), QGauss<dim>(degree + 1))
  // 
  , dof_handler(triangulation)
{
  E = 2e11;
  possion = 0.3;
  rho = 2701;
}
 
void Mixed_Mesh::make_grid()
{
  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  const std::string  input_file = "modle.msh";
  gridin.read_msh(input_file);
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}
 
void Mixed_Mesh::set_up_system()
{
  // correct oder
  // for (const auto &cell : dof_handler.active_cell_iterators())
  // {
  //   if (cell->reference_cell() == ReferenceCells::Pyramid)
  //     cell->set_active_fe_index(0);
  //   else if (cell->reference_cell() == ReferenceCells::Tetrahedron)
  //     cell->set_active_fe_index(1);
  //   else if (cell->reference_cell() == ReferenceCells::Hexahedron)
  //     cell->set_active_fe_index(2);
  //   else
  //     Assert(false, ExcNotImplemented());
  // }
  // error oder
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->reference_cell() == ReferenceCells::Tetrahedron)
      cell->set_active_fe_index(0);
    else if (cell->reference_cell() == ReferenceCells::Pyramid)
      cell->set_active_fe_index(1);
    else if (cell->reference_cell() == ReferenceCells::Hexahedron)
      cell->set_active_fe_index(2);
    else
      Assert(false, ExcNotImplemented());
  }
  // 
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl; 
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
}

void Mixed_Mesh::assemble_system()
{
  hp::FEValues<dim> hp_fe_values(mapping,
                               fe,
                               quadrature_formula,
                               update_values | update_gradients |
                               update_JxW_values);
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  unsigned int dofs_per_cell;
  double lambda = E * possion / ((1 + possion) * (1 - 2 * possion));
  double mu = E / (2 * (1 + possion));
  Tensor<1, dim> rhs_values;
  rhs_values[0] = 1e6;
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    hp_fe_values.reinit(cell);
    const auto &fe_values = hp_fe_values.get_present_fe_values();

    dofs_per_cell = cell->get_fe().n_dofs_per_cell();
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);

    cell_matrix = 0;
    cell_rhs    = 0;
    for (const unsigned int i : fe_values.dof_indices())
    {
      const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
      for (const unsigned int j : fe_values.dof_indices())
      {
        const unsigned int component_j = cell->get_fe().system_to_component_index(j).first;
        for (const unsigned int q_point : fe_values.quadrature_point_indices())
        {
          cell_matrix(i, j) += 
              ((fe_values.shape_grad(i, q_point)[component_i] *
                fe_values.shape_grad(j, q_point)[component_j] *
                lambda) +
                (fe_values.shape_grad(i, q_point)[component_j] *
                fe_values.shape_grad(j, q_point)[component_i] *
                mu) +
                ((component_i == component_j) ? (fe_values.shape_grad(i, q_point) *
                                                fe_values.shape_grad(j, q_point) *
                                                mu)
                                              : 0)) *
              fe_values.JxW(q_point);
        }
      }
    }
    for (const unsigned int i : fe_values.dof_indices())
    {
      const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
      for (const unsigned int q_point : fe_values.quadrature_point_indices())
      {
        cell_rhs(i) += fe_values.shape_value(i, q_point) *
                        rhs_values[component_i] *
                        fe_values.JxW(q_point);
      }
    }
    cell->get_dof_indices(local_dof_indices);
    for (const unsigned int i : fe_values.dof_indices())
      for (const unsigned int j : fe_values.dof_indices())
        system_matrix.add(local_dof_indices[i],
                          local_dof_indices[j],
                          cell_matrix(i, j));

    for (const unsigned int i : fe_values.dof_indices())
      system_rhs(local_dof_indices[i]) += cell_rhs(i);
  }
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 1, Functions::ZeroFunction<dim>(dim), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                    system_matrix,
                                    solution,
                                    system_rhs);  
}

void Mixed_Mesh::solve()
{
  SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

void Mixed_Mesh::output_results()
{
  DataOut<dim> data_out;
 
  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(mapping, 2);
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
}

void Mixed_Mesh::run()
{
  make_grid();
  set_up_system();
  assemble_system();
  solve();
  output_results();
}

int main()
{
  deallog.depth_console(2);
  Mixed_Mesh Mixed_Mesh_Problem;
  Mixed_Mesh_Problem.run();
  return 0;
}