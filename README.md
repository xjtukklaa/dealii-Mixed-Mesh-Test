# dealii-Mixed-Mesh-Test
dealii Mixed Mesh Test
```c++
// correct oder 
// : mapping(MappingFE<dim>(FE_PyramidP<dim>(degree)), MappingFE<dim>(FE_SimplexP<dim>(degree)), MappingFE<dim>(FE_Q<dim>(degree)))
// , fe(FESystem<dim>(FE_PyramidP<dim>(degree)^3), FESystem<dim>(FE_SimplexP<dim>(degree)^3), FESystem<dim>(FE_Q<dim>(degree)^3))
// , quadrature_formula(QGaussPyramid<dim>(degree + 1), QGaussSimplex<dim>(degree + 1), QGauss<dim>(degree + 1))
// error oder
: mapping(MappingFE<dim>(FE_SimplexP<dim>(degree)), MappingFE<dim>(FE_PyramidP<dim>(degree)), MappingFE<dim>(FE_Q<dim>(degree)))
, fe(FESystem<dim>(FE_SimplexP<dim>(degree)^3), FESystem<dim>(FE_PyramidP<dim>(degree)^3), FESystem<dim>(FE_Q<dim>(degree)^3))
, quadrature_formula(QGaussSimplex<dim>(degree + 1), QGaussPyramid<dim>(degree + 1), QGauss<dim>(degree + 1))
```
```c++
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
```
comment this two parts.
