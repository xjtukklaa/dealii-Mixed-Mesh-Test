# dealii-Mixed-Mesh-Test
In order to read mixed mesh, the DEALII_WITH_GMSH_API should be on 
```c++
// correct order 
// : mapping(MappingFE<dim>(FE_PyramidP<dim>(degree)), MappingFE<dim>(FE_SimplexP<dim>(degree)), MappingFE<dim>(FE_Q<dim>(degree)))
// , fe(FESystem<dim>(FE_PyramidP<dim>(degree)^3), FESystem<dim>(FE_SimplexP<dim>(degree)^3), FESystem<dim>(FE_Q<dim>(degree)^3))
// , quadrature_formula(QGaussPyramid<dim>(degree + 1), QGaussSimplex<dim>(degree + 1), QGauss<dim>(degree + 1))
// error order
: mapping(MappingFE<dim>(FE_SimplexP<dim>(degree)), MappingFE<dim>(FE_PyramidP<dim>(degree)), MappingFE<dim>(FE_Q<dim>(degree)))
, fe(FESystem<dim>(FE_SimplexP<dim>(degree)^3), FESystem<dim>(FE_PyramidP<dim>(degree)^3), FESystem<dim>(FE_Q<dim>(degree)^3))
, quadrature_formula(QGaussSimplex<dim>(degree + 1), QGaussPyramid<dim>(degree + 1), QGauss<dim>(degree + 1))
```
```c++
// correct order
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
// error order
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
```
comment and uncomment this two parts.
the output of correct order (release mode) as follows:
```
Number of active cells: 11728
Number of degrees of freedom: 15207
DEAL:cg::Starting value 3.24343e+06
DEAL:cg::Convergence step 287 value 0.0317760
```
the output of error order (release mode) as follows:
```
Number of active cells: 11728
Number of degrees of freedom: 20430
DEAL:cg::Starting value 3.21878e+06
DEAL:cg::Failure step 1000 value 4.64660e+13
terminate called after throwing an instance of 'dealii::SolverControl::NoConvergence'
```
the correct dofs should be 15207.
