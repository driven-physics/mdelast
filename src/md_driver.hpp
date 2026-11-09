#ifndef MD_DRIVER_HPP_
#define MD_DRIVER_HPP_

#include <vector>
#include <array>
#include <cassert>
#include <stdio.h>

#include "potential.hpp"
#include "eam.hpp"
#include "vector3d.hpp"
#include "thermostat.hpp"
#include "barostat.hpp"
#include "barostat_flex.hpp"

#include <Eigen/Dense>

struct MDDriver {

    enum Ensamble { NVE, NVT, NPT, NPT_FLEX };
    enum ThermostatCoupling { THERMOSTAT_KINETIC, THERMOSTAT_CONFIGURATIONAL, THERMOSTAT_AVERAGE };

    static constexpr std::array<unsigned int, 6> voigt_first_index = {0, 1, 2, 1, 0, 0};
    static constexpr std::array<unsigned int, 6> voigt_second_index = {0, 1, 2, 2, 2, 1};
    using Matrix6d = Eigen::Matrix<double, 6, 6>;
    using Vector6d = Eigen::Matrix<double, 6, 1>;

    struct Box {
        Box(Eigen::Matrix3d const &h_) : h(h_) {}
        Eigen::Matrix3d h;
        double volume() const {
            return h.determinant();
        }
        double length() const {
            return h.trace() / 3.;
        }
    };

    MDDriver(std::vector<Vector3d> const &positions_, std::vector<Vector3d> const &velocities_,
             std::vector<potential::PotentialPtr> const &pot_,
             std::vector<potential::EAMPtr> const &eam_, Ensamble const ensamble_,
             ThermostatCoupling const thermostat_coupling_, Thermostat const &thermostat_,
             Barostat const &barostat_, Box const &box_, std::vector<double> const &rg_,
             double const Nf_, std::vector<double> const &mass_,
             std::vector<unsigned int> const &types_, double const dt_global_,
             bool const neigbors_list_, double const r_skin_, double const fatal_conservation_);

    static double nint(double const x) {
        return static_cast<double>(x < 0. ? static_cast<int>(x - 0.5) : static_cast<int>(x + 0.5));
    }

    Vector3d apply_bc(Vector3d const &c) {
        return c - Vector3d(nint(c.x), nint(c.y), nint(c.z));
    }

    Vector3d get_scaled_coordinates(Vector3d const &p);
    Vector3d get_position(Vector3d const &c);

    void set_kinetic_energies();
    void set_momentum_tot();
    void calculate_configurational_temperature();
    void calculate_thermodynamics();

    void calculate_forces();
    void calculate_eam();
    void calculate_pair_correlation();
    void calculate_pair_correlations();
    void calculate_partial_pair_correlation();
    void integrate();

    void calculate_neighbors();
    std::pair<unsigned int, unsigned int>
    get_neighbors_list_indices(int const particle_index) const;

    void U_NHC_thermo_propagator(double const dt);
    void U4_NHC_thermo_propagator(double const dt, bool const outward);
    void U3_NHC_thermo_propagator(double const dt);

    void U_NHC_baro_propagator(double const dt);
    void U4_NHC_baro_propagator(double const dt, bool const outward);
    void U3_NHC_baro_propagator(double const dt);

    void U2_propagator(double const dt);
    void U2p_propagator(double const dt);
    void U1_propagator(double const dt);

    bool const neighbors_list;
    double const r_skin;
    std::vector<int> point;
    std::vector<int> neighbors;

    unsigned int const Nparticles;
    unsigned int const Ntypes;
    double const Nf;                   // number of degrees of freedom
    double const alpha;                // 1+3/Nf, used in NPT calculations
    std::vector<Vector3d> positions;   // real positions
    std::vector<Vector3d> coordinates; // scaled on [-0.5, 0.5]
    std::vector<Vector3d> velocities;
    std::vector<Vector3d> forces;
    std::vector<double> rho_eam;      // EAM densities
    std::vector<double> F_eam;        // EAM energies F(rho_eam)
    std::vector<double> dFdrho_eam;   // EAM F'
    std::vector<double> d2Fdrho2_eam; // EAM F''
    std::vector<Vector6d> g_eam;
    std::vector<Vector3d> h_eam;
    std::vector<double> hsq_eam;
    std::vector<double> energies_pot;
    std::vector<double> energies_kin;
    Vector3d momentum_tot; //-- total linear momentum

    std::vector<double> const rg; // cell cenetered, must be evenly spaced
    std::vector<double> gr;
    std::vector<std::vector<double>> gr_partial;
    double const drg;
    double const rmaxg;

    std::vector<potential::PotentialPtr> pot;
    double const r_cut_sq;
    double const r_cutr_skin_sq;
    std::vector<potential::EAMPtr> eam;
    bool const has_eam;

    Ensamble const ensamble;
    ThermostatCoupling const thermostat_coupling;
    Thermostat thermostat;
    Barostat barostat;

    Box const box_initial;
    Box box;

    std::vector<double> const mass;
    std::vector<unsigned int> const types;
    std::vector<unsigned int> Nparticles_per_type;

    double epot_eam_atom; // EAM energy per atom
    double epot_atom;
    double ekin_atom;
    double etot_atom;

    double ekin_tot;
    double etot_cons_initial;
    double etot_cons;
    double etot_cons_err;

    double virial;
    Vector6d born_stress_tensor;     // Born stress tensor
    Matrix6d born_elasticity_tensor; // Born term of the elastic constant

    double total_laplacian;
    double total_force2;
    double temperature_kin;
    double temperature_conf_inverse;
    double pressure;
    Vector6d stress_tensor;

    double const dt_global;
    double time;
    unsigned long long int cycle;

    double const
        fatal_conservation; // fatal current or initial conservation for which we exit if occured

    // When false, the next integrate() computes only what the dynamics needs (forces, energies,
    // virial and stress). The Born elasticity tensor, the configurational-temperature Laplacian
    // and |F|^2, and the EAM sums g and h that feed them are skipped and keep their last values.
    // The trajectory does not depend on this flag (unless the thermostat is coupled to T_conf).
    bool sample_observables = true;
};

#endif /* MD_DRIVER_HPP_ */
