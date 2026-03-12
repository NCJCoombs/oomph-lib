// Driver code for benchmark of axisym non-linear solid equations
// in cylindrical coordinates: Thick spherical shell with external pressure
// Compare against analytical solution for an incompressible material

// Generic routines
#include "generic.h"

// Axisymmetric non-linear solid equations in cylindrical coords
#include "axisym_solid.h"

// Constitutive laws
#include "neo_hookean_constitutive.h"

// The mesh: annular
#include "meshes/annular_mesh.h"

using namespace std;
using namespace oomph;

//=======start_namespace==========================================
/// Global physical variables
//================================================================
namespace GlobalPhysicalVariables
{
  // Poisson ratio, should be 0.5 since using incompressible formulation
  double nu = 0.49;

  // The external pressure
  double ext_pressure = 0.05;

  // Traction function for appyling a pressure on the outer surface
  void external_pressure_fct(const double& time,
                             const Vector<double>& x,
                             const Vector<double>& n,
                             Vector<double>& traction)
  {
    traction[0] = -1.0 * ext_pressure * n[0];
    traction[1] = -1.0 * ext_pressure * n[1];
  }
} // namespace GlobalPhysicalVariables

//=======start_namespace==========================================
/// Global simulation settings
//================================================================
namespace GlobalSimSettings
{
  // Folder for output
  string result_folder = "RESLT";

  /// Flag for using a pressure formulation
  unsigned use_pressure_formulation = 0;

  /// If using a pressure formulation, is the solid incompressible?
  unsigned incompressible = 0;

  // Number of elements in radial direction
  unsigned n_radial = 10;

  // Number of elements in azimuthal direction
  unsigned n_azimuthal = 10;

  // Inner radius
  double a = 1.0;

  // Thickness
  double h = 0.5;
} // namespace GlobalSimSettings

//======================start_mesh================================
/// Upgrade the mesh to a solidmesh
//================================================================
template<class ELEMENT>
class ElasticTwoDAnnularMesh : public virtual TwoDAnnularMesh<ELEMENT>,
                               public virtual SolidMesh
{
public:
  /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
  /// ones so that the initial configuration is the stress-free one.
  ElasticTwoDAnnularMesh<ELEMENT>(
    const bool& periodic,
    const double& azimuthal_fraction,
    const unsigned& ntheta,
    const unsigned& nr,
    const double& a,
    const double& h,
    const double& phi,
    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    : RectangularQuadMesh<ELEMENT>(
        ntheta, nr, 1.0, 1.0, periodic, time_stepper_pt),
      TwoDAnnularMesh<ELEMENT>(
        periodic, azimuthal_fraction, ntheta, nr, a, h, phi, time_stepper_pt)
  {
    // Make the current configuration the undeformed one by
    // setting the nodal Lagrangian coordinates to their current
    // Eulerian ones
    set_lagrangian_nodal_coordinates();
  }
};

//==============start_problem=========================================
/// Unstructured solid problem
//====================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
class ShellProblem : public Problem
{
public:
  /// Constructor:
  ShellProblem();

  /// Destructor
  ~ShellProblem() {}

  /// Doc the solution
  void doc_solution();

  /// Compute pressure integral
  void compute_pressure_integral();

private:
  /// Create the boundary elements
  void create_boundary_elements();

  /// Delete boundary elements and wipe the associated mesh
  void delete_boundary_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Kill surface element
      delete Surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Surface_mesh_pt->flush_element_and_node_storage();
  } // end of delete_boundary_elements

  /// Complete the problem setup (attach function/physical parameter pointers)
  void complete_problem_setup();

  /// Bulk mesh
  ElasticTwoDAnnularMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;

  /// Pointer to auxillary surface mesh, which exists solely for post-processing
  Mesh* Aux_surface_mesh_pt;

  /// Pointer to constitutive law
  ConstitutiveLaw* Constitutive_law_pt;

  // DocInfo object
  DocInfo doc_info;
};

//===============start_constructor========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
ShellProblem<ELEMENT, INTERFACE_ELEMENT>::ShellProblem()
{
#ifdef OOMPH_HAS_MUMPS
  // Use mumps if available
  linear_solver_pt() = new MumpsSolver;
#endif

  // Set output directory
  doc_info.set_directory(GlobalSimSettings::result_folder);

  // Initialise counter for solutions
  doc_info.number() = 0;

  // Create the mesh: Use constuctor which allows rotation
  bool periodic = false;
  double azimuthal_fraction = 0.25;
  double phi = MathematicalConstants::Pi;
  Bulk_mesh_pt =
    new ElasticTwoDAnnularMesh<ELEMENT>(periodic,
                                        azimuthal_fraction,
                                        GlobalSimSettings::n_azimuthal,
                                        GlobalSimSettings::n_radial,
                                        GlobalSimSettings::a,
                                        GlobalSimSettings::h,
                                        phi);

  // Create the "surface mesh" that will contain only the interface
  // elements. The constructor just creates the mesh without giving
  // it any elements, nodes, etc.
  Surface_mesh_pt = new Mesh;
  Aux_surface_mesh_pt = new Mesh;

  // Create the boundary elements
  create_boundary_elements();

  // Add the two submeshes
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine all sub-meshes into a single mesh
  build_global_mesh();

  // Define a constitutive law from a strain energy density function.
  // Version depends on the incompressibility flag
  StrainEnergyFunction* Strain_energy_function_pt;
  if (GlobalSimSettings::incompressible)
  {
    Strain_energy_function_pt =
      new IncompressibleNeoHookean(&GlobalPhysicalVariables::nu);
  }
  else
  {
    Strain_energy_function_pt = new NeoHookean(&GlobalPhysicalVariables::nu);
  }
  Constitutive_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(Strain_energy_function_pt);

  // Complete the problem set up
  complete_problem_setup();

  // Assign eqn numbers
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // end constructor

//========start_of_complete_problem_setup==================================
/// Assign function/parameter pointers, set boundary conditions etc.
//=========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void ShellProblem<ELEMENT, INTERFACE_ELEMENT>::complete_problem_setup()
{
  // Determine number of bulk elements in mesh
  const unsigned n_element_bulk = Bulk_mesh_pt->nelement();

  // Loop over the bulk elements
  for (unsigned e = 0; e < n_element_bulk; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the constitutive law
    el_pt->constitutive_law_pt() = Constitutive_law_pt;

    // Cast to element with pressure dofs
    QAxisymmetricCylindricalPVDWithPressureElement* el_with_pres_pt =
      dynamic_cast<QAxisymmetricCylindricalPVDWithPressureElement*>(el_pt);

    // Can set compressible/incompressible condition here
    if (el_with_pres_pt != 0)
    {
      if (GlobalSimSettings::incompressible)
      {
        el_with_pres_pt->set_incompressible();
      }
      else
      {
        el_with_pres_pt->set_compressible();
      }
    }
  } // End of loop over bulk elements

  // Determine number of 1D interface elements in mesh
  const unsigned n_interface_element = Surface_mesh_pt->nelement();

  // Loop over the interface elements
  for (unsigned e = 0; e < n_interface_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    INTERFACE_ELEMENT* el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(e));

    // Set the traction function for applying pressure
    el_pt->traction_fct_pt() = &GlobalPhysicalVariables::external_pressure_fct;

  } // End of loop over interface elements

  // Set the Dirichlet boundary conditions for this problem
  // Find number of nodes on boundary 3
  unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(3);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position on symmetry axis only
    Bulk_mesh_pt->boundary_node_pt(3, n)->pin_position(0);
  }
  n_boundary_node = Bulk_mesh_pt->nboundary_node(1);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position on symmetry axis only
    Bulk_mesh_pt->boundary_node_pt(1, n)->pin_position(1);
  }
}

//============start_of_create_boundary_elements===============
/// Create boundary elements
//=======================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void ShellProblem<ELEMENT, INTERFACE_ELEMENT>::create_boundary_elements()
{
  // Loop over elements on boundary 2
  unsigned n_element = Bulk_mesh_pt->nboundary_element(2);
  for (unsigned e = 0; e < n_element; e++)
  {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT* bulk_element_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(2, e));

    // Find the index of the face of element e along boundary b
    int face_index = Bulk_mesh_pt->face_index_at_boundary(2, e);

    // Create the interface element
    INTERFACE_ELEMENT* interface_element_pt =
      new INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(2);
  }

  // Loop over elements on boundary 1
  n_element = Bulk_mesh_pt->nboundary_element(1);
  for (unsigned e = 0; e < n_element; e++)
  {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT* bulk_element_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(1, e));

    // Find the index of the face of element e along boundary b
    int face_index = Bulk_mesh_pt->face_index_at_boundary(1, e);

    // Create the interface element
    INTERFACE_ELEMENT* interface_element_pt =
      new INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Aux_surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(1);
  }
} // end of create_boundary_elements

//========================================================================
/// Compute the integral which should equal the pressure
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void ShellProblem<ELEMENT, INTERFACE_ELEMENT>::compute_pressure_integral()
{
  // Effective shear modulus
  double mu = 0.5 / (1.0 + GlobalPhysicalVariables::nu);

  // Initialize integral
  double intgrl = 0.0;

  // Loop over elements in aux surface mesh
  unsigned n_el = Aux_surface_mesh_pt->nelement();
  for (unsigned e = 0; e < n_el; e++)
  {
    // Cast to element
    INTERFACE_ELEMENT* el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(Aux_surface_mesh_pt->element_pt(e));

    // Loop over integration points
    unsigned n_intpt = el_pt->integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = el_pt->integral_pt()->weight(ipt);

      // Shape fns at integration point
      Shape psi(el_pt->nnode());
      el_pt->shape_at_knot(ipt, psi);

      // Calculate the global position and lagrangian coordinate
      Vector<double> interpolated_x(2, 0.0), interpolated_xi(2, 0.0);
      for (unsigned l = 0; l < el_pt->nnode(); l++)
      {
        // Loop over the number of lagrangian coordinates (2)
        for (unsigned i = 0; i < 2; i++)
        {
          // Calculate the global position
          interpolated_x[i] += el_pt->nodal_position(l, i) * psi(l);
          interpolated_xi[i] += el_pt->lagrangian_position(l, i) * psi(l);
        }
      }

      // Calculate Jacobian of mapping and pre-multiply with integral weight
      double J = el_pt->J_eulerian_at_knot(ipt);
      double W = w * J;

      // The principal stretch
      double lambda = interpolated_x[0] / interpolated_xi[0];

      // Add integral contributions
      intgrl += (lambda - pow(lambda, -5.0)) / interpolated_xi[0] * W;
    }
  }

  // Pressure
  double p = -2.0 * mu * intgrl;

  std::cout << "Computed pressure: " << p << std::endl;
}

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
void ShellProblem<ELEMENT, INTERFACE_ELEMENT>::doc_solution()
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 2;

  sprintf(filename,
          "%s/bulk_soln%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file.close();

  // Increment the counter
  doc_info.number()++;
}

//===========start_main===================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main(int argc, char* argv[])
{
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
    "--use_solid_pressure",
    &GlobalSimSettings::use_pressure_formulation,
    "Use solid pressure");

  CommandLineArgs::specify_command_line_flag(
    "--incompressible", &GlobalSimSettings::incompressible, "Incompressible");

  CommandLineArgs::specify_command_line_flag(
    "--poisson", &GlobalPhysicalVariables::nu, "Poisson");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  // Create the problem
  if (GlobalSimSettings::use_pressure_formulation)
  {
    // Build the problem with an additional degree of freedom for pressure
    ShellProblem<QAxisymmetricCylindricalPVDWithPressureElement,
                 AxisymmetricCylindricalSolidTractionElement<
                   QAxisymmetricCylindricalPVDWithPressureElement>>
      problem;

    problem.newton_solve();
    problem.doc_solution();
    problem.compute_pressure_integral();
  }
  else
  {
    // Standard build
    ShellProblem<QAxisymmetricCylindricalPVDElement<3>,
                 AxisymmetricCylindricalSolidTractionElement<
                   QAxisymmetricCylindricalPVDElement<3>>>
      problem;

    problem.newton_solve();
    problem.doc_solution();
    problem.compute_pressure_integral();
  }
} // end main
