// Neo-Hookean constitutive law

#include "constitutive.h"

namespace oomph
{
  //=====================================================================
  /// Neo Hookean constitutive law in terms of poisson ratio
  //====================================================================
  class NeoHookean : public StrainEnergyFunction
  {
  public:
    /// Constructor takes the pointers to the constitutive parameter:
    /// Poisson's ratio. Young's modulus is set
    /// to 1, implying that it has been used to scale the stresses
    NeoHookean(double* nu_pt) : StrainEnergyFunction(), Nu_pt(nu_pt) {}

    /// Virtual destructor (empty)
    virtual ~NeoHookean() {}

    /// Return the strain energy in terms of strain tensor
    double W(const DenseMatrix<double>& gamma)
    {
      return StrainEnergyFunction::W(gamma);
    }

    /// Return the strain energy in terms of the strain invariants
    double W(const Vector<double>& I)
    {
      double nu = *Nu_pt;
      return 0.5 *
             (0.5 * (I[0] - 3.0) - 0.5 * log(I[2]) +
              (nu / (1.0 - 2.0 * nu)) * 0.25 * log(I[2]) * log(I[2])) /
             (1.0 + nu);
    }

    /// Return the derivatives of the strain energy function with
    /// respect to the strain invariants
    void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      double nu = *Nu_pt;
      dWdI[0] = 0.25 / (1.0 + nu);
      dWdI[1] = 0.0;
      dWdI[2] = 0.5 * (-0.5 + (nu / (1.0 - 2.0 * nu)) * 0.5 * log(I[2])) /
                ((1.0 + nu) * I[2]);
    }

    /// Pure virtual function in which the user must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode. False.
    bool requires_incompressibility_constraint()
    {
      return false;
    }

  private:
    /// Poisson's ratio
    double* Nu_pt;
  };

  //=====================================================================
  /// Neo Hookean constitutive law (incompressible version)
  //====================================================================
  class IncompressibleNeoHookean : public StrainEnergyFunction
  {
  public:
    /// Constructor takes the pointers to the constitutive parameter:
    /// Poisson's ratio. Young's modulus is set
    /// to 1, implying that it has been used to scale the stresses
    IncompressibleNeoHookean(double* nu_pt)
      : StrainEnergyFunction(), Nu_pt(nu_pt)
    {
    }

    /// Virtual destructor (empty)
    virtual ~IncompressibleNeoHookean() {}

    /// Return the strain energy in terms of strain tensor
    double W(const DenseMatrix<double>& gamma)
    {
      return StrainEnergyFunction::W(gamma);
    }

    /// Return the strain energy in terms of the strain invariants
    double W(const Vector<double>& I)
    {
      double nu = *Nu_pt;
      return 0.25 * (I[0] - 3.0) / (1.0 + nu);
    }

    /// Return the derivatives of the strain energy function with
    /// respect to the strain invariants
    void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      double nu = *Nu_pt;
      dWdI[0] = 0.25 / (1.0 + nu);
      dWdI[1] = 0.0;
      dWdI[2] = 0.0;
    }

    /// Pure virtual function in which the user must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode. True.
    bool requires_incompressibility_constraint()
    {
      return true;
    }

  private:
    /// Poisson's ratio
    double* Nu_pt;
  };
} // namespace oomph
