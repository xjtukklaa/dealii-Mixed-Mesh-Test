# dealii-Mixed-Mesh-Test
In order to read mixed mesh, the DEAL_II_GMSH_WITH_API should be on 
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

Futher Test In Debug mode(errot order)
```c++
Number of active cells: 11728

--------------------------------------------------------
An error occurred in line <657> of file </home/zxl/workspace/dealii/tmp/unpack/deal.II-v9.5.2/source/fe/fe_simplex_p.cc> in function
    dealii::FiniteElementDomination::Domination dealii::FE_SimplexP<dim, spacedim>::compare_for_domination(const dealii::FiniteElement<dim, spacedim>&, unsigned int) const [with int dim = 3; int spacedim = 3]
The violated condition was: 
    false
Additional information: 
    You are trying to use functionality in deal.II that is currently not
    implemented. In many cases, this indicates that there simply didn't
    appear much of a need for it, or that the author of the original code
    did not have the time to implement a particular case. If you hit this
    exception, it is therefore worth the time to look into the code to
    find out whether you may be able to implement the missing
    functionality. If you do, please consider providing a patch to the
    deal.II development sources (see the deal.II website on how to
    contribute).

Stacktrace:
-----------
#0  /home/zxl/workspace/dealii/deal.II-v9.5.2/lib/libdeal_II.g.so.9.5.2: dealii::FE_SimplexP<3, 3>::compare_for_domination(dealii::FiniteElement<3, 3> const&, unsigned int) const
#1  /home/zxl/workspace/dealii/deal.II-v9.5.2/lib/libdeal_II.g.so.9.5.2: dealii::FESystem<3, 3>::compare_for_domination(dealii::FiniteElement<3, 3> const&, unsigned int) const
#2  /home/zxl/workspace/dealii/deal.II-v9.5.2/lib/libdeal_II.g.so.9.5.2: dealii::hp::FECollection<3, 3>::find_dominating_fe(std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> > const&, unsigned int) const
...
```
correct order
```c++
Number of active cells: 11728
Number of degrees of freedom: 15207
```
source/fe/fe_pyramid_p.cc<196>
```c++
template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_PyramidP<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // (if fe_other is derived from FE_SimplexDGP)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(&fe_other) !=
        nullptr)
      // there are no requirements between continuous and discontinuous
      // elements
      return FiniteElementDomination::no_requirements;

  // vertex/line/face domination
  // (if fe_other is not derived from FE_SimplexDGP)
  // & cell domination
  // ----------------------------------------
  if (const FE_PyramidP<dim, spacedim> *fe_pp_other =
        dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_pp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_pp_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_SimplexP<dim, spacedim> *fe_p_other =
             dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_p_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_p_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q<dim, spacedim> *fe_q_other =
             dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim, spacedim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}
```

source/fe/fe_simplex_p.cc<804>
```c++
template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_SimplexP<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // (if fe_other is derived from FE_SimplexDGP)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(&fe_other) !=
        nullptr)
      // there are no requirements between continuous and discontinuous
      // elements
      return FiniteElementDomination::no_requirements;

  // vertex/line/face domination
  // (if fe_other is not derived from FE_SimplexDGP)
  // & cell domination
  // ----------------------------------------
  if (const FE_SimplexP<dim, spacedim> *fe_p_other =
        dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_p_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_p_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q<dim, spacedim> *fe_q_other =
             dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim, spacedim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}
```
In fe_simplex_P.cc, the funciton--compare_for_domination won't check the pyramid cell type. this maybe the one reason for this problem. but i not very sure.
