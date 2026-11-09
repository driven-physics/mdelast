#include "md_driver.hpp"

#include <limits>
#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <boost/math/special_functions/pow.hpp>

#include "units.hpp"
#include "numerics.hpp"

using boost::math::pow;

MDDriver::MDDriver(std::vector<Vector3d> const &positions_,
                   std::vector<Vector3d> const &velocities_,
                   std::vector<potential::PotentialPtr> const &pot_,
                   std::vector<potential::EAMPtr> const &eam_, Ensamble const ensamble_,
                   ThermostatCoupling const thermostat_coupling_, Thermostat const &thermostat_,
                   Barostat const &barostat_, Box const &box_, std::vector<double> const &rg_,
                   double const Nf_, std::vector<double> const &mass_,
                   std::vector<unsigned int> const &types_, double const dt_global_,
                   bool const neighbors_list_, double const r_skin_,
                   double const fatal_conservation_)
    : neighbors_list(neighbors_list_), r_skin(r_skin_),
      point(neighbors_list ? positions_.size() + 1 : 0, 0),
      neighbors(neighbors_list ? positions_.size() * 40 : 0, 0),

      Nparticles(positions_.size()), Ntypes(mass_.size()), Nf(Nf_), alpha(1. + 3. / Nf),
      positions(positions_), coordinates(positions_.size(), Vector3d(0., 0., 0.)),
      velocities(velocities_), forces(positions_.size(), Vector3d(0., 0., 0.)),
      rho_eam(positions_.size(), 0.), F_eam(positions_.size(), 0.),
      dFdrho_eam(positions_.size(), 0.), d2Fdrho2_eam(positions_.size(), 0.),
      g_eam(positions_.size()), h_eam(positions_.size(), Vector3d(0., 0., 0.)),
      hsq_eam(positions.size(), 0.), energies_pot(positions_.size(), 0.),
      energies_kin(positions_.size(), 0.), momentum_tot(0., 0., 0.),

      rg(rg_), gr(rg_.size(), 0.),
      gr_partial(Ntypes * (Ntypes + 1) / 2, std::vector<double>(rg_.size(), 0.)),
      drg(rg_[1] - rg_[0]), rmaxg(rg[rg_.size() - 1] + 0.5 * (rg_[1] - rg_[0])), pot(pot_),
      r_cut_sq(pow<2>(pot[0]->rcut)), r_cutr_skin_sq(pow<2>(pot[0]->rcut * r_skin)), eam(eam_),
      has_eam(eam.size() > 0), ensamble(ensamble_), thermostat_coupling(thermostat_coupling_),
      thermostat(thermostat_), barostat(barostat_), box_initial(box_), box(box_),

      mass(mass_), types(types_), Nparticles_per_type(Ntypes, 0),

      epot_eam_atom(0.), epot_atom(0.), ekin_atom(0.), etot_atom(0.), ekin_tot(0.),
      etot_cons_initial(0.), etot_cons(0.), etot_cons_err(0.), virial(0.), born_stress_tensor(),
      born_elasticity_tensor(), total_laplacian(0.), total_force2(0.), temperature_kin(0.),
      temperature_conf_inverse(0.), pressure(0.), stress_tensor(), dt_global(dt_global_), time(0.),
      cycle(0), fatal_conservation(fatal_conservation_) {

    if (has_eam)
        printf("EAM object was given to MDDriver constructor\n");
    else
        printf("EAM object was not given to MDDriver constructor\n");

    assert(dt_global > 0.);
    assert(box_.h.isDiagonal());
    assert(positions.size() == Nparticles);
    assert(coordinates.size() == Nparticles);
    assert(velocities.size() == Nparticles);
    assert(forces.size() == Nparticles);
    assert(energies_pot.size() == Nparticles);
    assert(energies_kin.size() == Nparticles);
    assert(types.size() == Nparticles);
    if (has_eam)
        assert(eam.size() == Ntypes);
    assert(pot.size() == Ntypes * (Ntypes + 1) / 2);

    for (const auto &p : pot) {
        if (!p) {
        } // harmless use, prevents unused warning
        assert(p != nullptr);
    }
    for (const auto &p : eam) {
        if (!p) {
        } // harmless use, prevents unused warning
        assert(p != nullptr);
    }

    assert(r_skin >= 1);

    assert(rg.size() > 1);
    assert(rg[0] > 0.); // rg is cell centered
    assert(drg > 0.);

    Eigen::Vector3d edge = box.h.diagonal();

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto const p = positions[i];
        if (std::abs(p.x) >= 0.5 * edge.x() or std::abs(p.y) >= 0.5 * edge.y() or
            std::abs(p.z) >= 0.5 * edge.z()) {
            printf("bad initial position\n");
            exit(1);
        }

        coordinates[i] = get_scaled_coordinates(p);

        if (types[i] >= Ntypes) {
            printf("bad type\n");
            exit(1);
        }
        Nparticles_per_type[types[i]] += 1;
    }

    born_stress_tensor.fill(0.);
    born_elasticity_tensor.fill(0.);
    stress_tensor.fill(0.);

    if (neighbors_list)
        calculate_neighbors();
    calculate_forces();
    set_kinetic_energies();
    calculate_thermodynamics();

    etot_cons_initial = etot_cons;

    printf("###############\n"
           "MDDriver constructor:\n"
           "N=%d\n"
           "Nf=%g\n"
           "Ntypes=%d\n"
           "dt_global=%g\n"
           "drg=%g\n"
           "r_skin=%f\n"
           "point_size=%ld\n"
           "mass:   ",
           Nparticles, Nf, Ntypes, dt_global, drg, r_skin, point.size());
    for (unsigned int i = 0; i < Ntypes; ++i) {
        printf("%g   ", mass[i]);
    }
    printf("\n");
    if (ensamble == NVE)
        printf("Esamble: NVE\n");
    else if (ensamble == NVT)
        printf("Esamble: NVT\n");
    else if (ensamble == NPT)
        printf("Esamble: NPT\n");
    else {
        printf("bad ensamble %d\n", ensamble);
        exit(1);
    }
    printf("###############\n");

    for (unsigned int i = 0; i < thermostat.nchain; ++i) {
        printf("THERMOSTAT chain %d/%d q=%g p_eta=%g eta=%g g=%g T=%g\n", i + 1, thermostat.nchain,
               thermostat.q[i], thermostat.p_eta[i], thermostat.eta[i], thermostat.Nf,
               thermostat.temperature);
    }
    auto const &baro_thermostat = barostat.thermostat;
    for (unsigned int i = 0; i < baro_thermostat.nchain; ++i) {
        printf("BAROSTAT THERMOSTAT chain %d/%d q=%g p_eta=%g eta=%g g=%g T=%g[K]\n", i + 1,
               baro_thermostat.nchain, baro_thermostat.q[i], baro_thermostat.p_eta[i],
               baro_thermostat.eta[i], baro_thermostat.Nf, baro_thermostat.temperature);
    }
    printf("BAROSTAT eps=%g peps=%g w=%g P=%g[GPa]\n", barostat.epsilon, barostat.p_epsilon,
           barostat.w, barostat.pressure / Units::GPa);
}

Vector3d MDDriver::get_scaled_coordinates(Vector3d const &p) {
    Vector3d ptag = p;
    ptag.x /= box.h(0, 0);
    ptag.y /= box.h(1, 1);
    ptag.z /= box.h(2, 2);

    Vector3d const c = apply_bc(ptag);
    assert(std::abs(c.x) < 0.5);
    assert(std::abs(c.y) < 0.5);
    assert(std::abs(c.z) < 0.5);
    return c;
}

Vector3d MDDriver::get_position(Vector3d const &c) {
    Vector3d p = apply_bc(c);
    p.x *= box.h(0, 0);
    p.y *= box.h(1, 1);
    p.z *= box.h(2, 2);

    assert(std::abs(p.x) <= 0.5 * box.h(0, 0));
    assert(std::abs(p.y) <= 0.5 * box.h(1, 1));
    assert(std::abs(p.z) <= 0.5 * box.h(2, 2));
    return p;
}

std::pair<unsigned int, unsigned int>
MDDriver::get_neighbors_list_indices(int const particle_index) const {
    if (neighbors_list) {
        return {point[particle_index], point[particle_index + 1]};
    }
    return {particle_index + 1, Nparticles};
}

void MDDriver::calculate_eam() {
    if (not has_eam)
        return;

    // g and h (and F'') feed only the Born tensor and the T_conf Laplacian.
    bool const observables = sample_observables;

    std::fill(rho_eam.begin(), rho_eam.end(), 0.);
    if (observables) {
        std::fill(h_eam.begin(), h_eam.end(), Vector3d(0., 0., 0.));
        for (auto &g : g_eam)
            g.fill(0.);
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto const &ci = coordinates[i];
        double &rho_i = rho_eam[i];
        Vector6d &g_i = g_eam[i];
        Vector3d &h_i = h_eam[i];
        unsigned int const type_i = types[i];

        auto const [start, end] = get_neighbors_list_indices(i);

        for (unsigned int k = start; k < end; ++k) {
            int j = neighbors_list ? neighbors[k] : k;
            unsigned int const type_j = types[j];

            auto const dp = get_position(ci - coordinates[j]);
            // Cache the squared distance to avoid redundant square roots.
            double const r2 = dp.length2();
            if (r2 > r_cut_sq)
                continue;
            double const dr = std::sqrt(r2);

            // Cache EAM spline values used by multiple terms.
            double const rho_ij = eam[type_j]->rho(dr); // contribution to rho_i from atom j
            double const rho_ji = eam[type_i]->rho(dr); // contribution to rho_j from atom i

            if (rho_ij > 1e25 || rho_ji > 1e25) {
                printf("dr = %f         i = %d         j = %d \n", dr, i, j);
                printf("eam[type_j]->rho(dr) = %f         eam[type_i]->rho(dr) = %f  \n",
                       rho_ij * 1e-25, rho_ji * 1e-25);
                printf("r_ij = (%f,    %f,    %f) \n", dp.x, dp.y, dp.z);
                printf("r_i = (%f,    %f,    %f) \n", ci.x, ci.y, ci.z);
                printf("r_j = (%f,    %f,    %f) \n", coordinates[j].x, coordinates[j].y,
                       coordinates[j].z);
            }
            rho_i += rho_ij;
            rho_eam[j] += rho_ji;

            if (not observables)
                continue;

            Vector6d &g_j = g_eam[j];
            // For rho_i we add rho_{type_j}(r_ij); for rho_j we add rho_{type_i}(r_ij).
            // Therefore derivatives must follow the same type mapping.
            double const drhodr_over_r_ij = eam[type_j]->drhodr(dr) / dr;
            double const drhodr_over_r_ji = eam[type_i]->drhodr(dr) / dr;
            for (unsigned int iv = 0; iv < 6; ++iv) {
                unsigned int const ii = voigt_first_index[iv];
                unsigned int const jj = voigt_second_index[iv];

                g_i(iv) += drhodr_over_r_ij * dp[ii] * dp[jj];
                g_j(iv) += drhodr_over_r_ji * dp[ii] * dp[jj];
            }

            h_i += drhodr_over_r_ij * dp;
            // dp = r_i - r_j. For atom j we need (r_j - r_i) = -dp.
            h_eam[j] -= drhodr_over_r_ji * dp;
        }
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        double const rho = rho_eam[i];
        unsigned int const type_i = types[i];
        F_eam[i] = eam[type_i]->F(rho);
        dFdrho_eam[i] = eam[type_i]->dFdrho(rho);
        if (observables) {
            d2Fdrho2_eam[i] = eam[type_i]->d2Fdrho2(rho);
            hsq_eam[i] = h_eam[i].length2();
        }
    }
}

void MDDriver::calculate_neighbors() {
    long unsigned int k = 0;
    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto const &ci = coordinates[i];
        point[i] = k;
        for (unsigned int j = i + 1; j < Nparticles; ++j) {
            auto const dp = get_position(ci - coordinates[j]);
            if (dp.length2() > r_cutr_skin_sq)
                continue;
            if (k < neighbors.size()) {
                neighbors[k] = j;
            } else {
                neighbors.push_back(j);
            }
            k += 1;
        }
    }
    // end marker: total number of stored neighbor entries
    point[Nparticles] = k;
}

void MDDriver::calculate_forces() {

    calculate_eam();

    // Forces, energies, the virial and the Born stress are needed every step. The Born
    // elasticity tensor and the Laplacian are computed only when observables are sampled.
    bool const observables = sample_observables;

    std::fill(forces.begin(), forces.end(), Vector3d(0., 0., 0.));
    std::fill(energies_pot.begin(), energies_pot.end(), 0.);
    virial = 0.;
    born_stress_tensor.fill(0.);
    if (observables) {
        born_elasticity_tensor.fill(0.);
        total_laplacian = 0.;
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto &fi = forces[i];
        auto const &ci = coordinates[i];
        double &epot_i = energies_pot[i];
        double const dFdrho_i = dFdrho_eam[i];
        double const d2Fdrho2_i = d2Fdrho2_eam[i];
        unsigned int const type_i = types[i];
        //-- add embedding energy
        epot_i += F_eam[i];

        auto const [start, end] = get_neighbors_list_indices(i);

        for (unsigned int k = start; k < end; ++k) {
            int j = neighbors_list ? neighbors[k] : k;
            unsigned int const type_j = types[j];
            // periodic B.C. - calculate forces from mirror particles
            auto const dp = get_position(ci - coordinates[j]);

            // Cache the squared distance to avoid redundant square roots.
            double const r2 = dp.length2();
            if (r2 > r_cut_sq)
                continue;

            // Packed triangular indexing expects (min(type_i,type_j), max(type_i,type_j))
            unsigned int const tmax = std::max(type_i, type_j);
            unsigned int const tmin = std::min(type_i, type_j);
            unsigned int const pot_index = tmax * (tmax + 1) / 2 + tmin;

            double const dr = std::sqrt(r2);
            auto const result = pot[pot_index]->calculate(dr);

            double dedr = result.dVdr;
            double Xij = 0.;
            if (observables) {
                Xij = result.d2Vdr2 - result.dVdr / dr;
                total_laplacian += 2. * (result.d2Vdr2 + 2. * result.dVdr / dr);
            }

            if (has_eam) {
                double const drhodr_i = eam[type_i]->drhodr(dr);
                double const drhodr_j = eam[type_j]->drhodr(dr);
                double const dFdrho_j = dFdrho_eam[j];

                dedr += dFdrho_i * drhodr_j + dFdrho_j * drhodr_i;

                if (observables) {
                    double const d2rhodr2_i = eam[type_i]->d2rhodr2(dr);
                    double const d2rhodr2_j = eam[type_j]->d2rhodr2(dr);

                    total_laplacian += 2. * (dFdrho_j * (d2rhodr2_i + 2. * drhodr_i / dr) +
                                             dFdrho_i * (d2rhodr2_j + 2. * drhodr_j / dr));
                    total_laplacian +=
                        d2Fdrho2_i * pow<2>(drhodr_j) + d2Fdrho2_eam[j] * pow<2>(drhodr_i);

                    // EAM contribution to the Born term (matches the single-component form):
                    // X_ij += (dF_i * d²ρ_j/dr² + dF_j * d²ρ_i/dr²) - (dF_i * dρ_j/dr + dF_j *
                    // dρ_i/dr)/r
                    Xij += dFdrho_i * (d2rhodr2_j - drhodr_j / dr) +
                           dFdrho_j * (d2rhodr2_i - drhodr_i / dr);
                }
            }

            Xij /= pow<2>(dr);

            Vector3d const df = -dedr / dr * dp;
            fi += df;
            forces[j] -= df;

            //-- potential energy
            epot_i += 0.5 * result.V;
            energies_pot[j] += 0.5 * result.V;

            //-- virial
            virial += dr * dedr;

            // virial tensor
            for (unsigned int iv = 0; iv < 6; ++iv) {
                unsigned int const ii = voigt_first_index[iv];
                unsigned int const jj = voigt_second_index[iv];
                double const num = dp[ii] * dp[jj];
                born_stress_tensor(iv) += dedr * num / dr;

                if (not observables)
                    continue;
                for (unsigned int ivp = iv; ivp < 6; ++ivp) {
                    unsigned int const iip = voigt_first_index[ivp];
                    unsigned int const jjp = voigt_second_index[ivp];
                    double const num2 = num * dp[iip] * dp[jjp];
                    born_elasticity_tensor(iv, ivp) += num2 * Xij;
                }
            }
        }

        if (has_eam and observables) {
            Vector6d const &g_i = g_eam[i];
            for (unsigned int iv = 0; iv < 6; ++iv) {
                for (unsigned int ivp = iv; ivp < 6; ++ivp) {
                    born_elasticity_tensor(iv, ivp) += d2Fdrho2_i * g_i(iv) * g_i(ivp);
                }
            }

            total_laplacian += d2Fdrho2_i * hsq_eam[i]; /// 4.;
        }
    }

    virial /= -3.;

    double const vir_trace =
        -(born_stress_tensor(0) + born_stress_tensor(1) + born_stress_tensor(2));
    double const v_err = std::abs(virial / (vir_trace / 3.) - 1.);
    if (v_err > 1e-7) {
        printf("FATAL VIRIALS err=%g %1.16E %1.16E\n", v_err, virial, vir_trace / 3.);
    }
}

//---- U2: velocity Verlet kick step propagator
void MDDriver::U2_propagator(double const dt) {

    double c = 1., e = 1.;

    if (ensamble == NPT) {

        double const x = alpha * dt * barostat.p_epsilon / barostat.w;

        auto const res = numerics::expx_operator(x);
        e = res.expx;
        c = res.op;
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto &v = velocities[i];
        v = v * e + c * dt * forces[i] / mass[types[i]];
    }
    set_kinetic_energies(); // the kinetic energy is used in the NHC U4 propagator
    calculate_configurational_temperature();
}

void MDDriver::U2p_propagator(double const dt) {
    assert(ensamble == NPT);
    //----- Eq. 31 in: T.-Q. Yu et al. / Chemical Physics 370 (2010) 294305
    double const pv_baro = barostat.pressure * box.volume();
    barostat.p_epsilon += (2. * alpha * ekin_tot + 3. * virial - 3. * pv_baro) * dt;
}

//---- U1: velocity Verlet drift step propagator
void MDDriver::U1_propagator(double const dt) {

    double cf = 1.;

    if (ensamble == NPT) {
        double const x = dt * barostat.p_epsilon / barostat.w;
        barostat.epsilon += x;
        cf = numerics::expx_operator(x).op;
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        Vector3d dc = cf * dt * velocities[i];
        dc.x /= box.h(0, 0);
        dc.y /= box.h(1, 1);
        dc.z /= box.h(2, 2);

        // timestep must be small enough so that this is true:
        assert(std::abs(dc.x) < 0.01);
        assert(std::abs(dc.y) < 0.01);
        assert(std::abs(dc.z) < 0.01);

        auto &c = coordinates[i];
        c += dc;
        positions[i] = get_position(c);
    }

    if (ensamble == NPT) {
        double const expe = std::exp(barostat.epsilon);
        double const vol_old = box.volume();
        box.h = box_initial.h * expe;

        if (std::abs(box.volume() / vol_old - 1.) > 0.05) {
            printf("FATAL: volume change too large... expe=%g V0=%g Vold=%g Vnew=%g\n", expe,
                   box_initial.volume(), vol_old, box.volume());
            exit(1);
        }
    }
}

void MDDriver::U3_NHC_thermo_propagator(double const dt) {
    assert(ensamble == NVT or ensamble == NPT);
    thermostat.U3_propagator(dt, velocities.data(), velocities.size());
    set_kinetic_energies(); // the kinetic energy is used in the NHC U4 propagator
}

void MDDriver::U3_NHC_baro_propagator(double const dt) {
    assert(ensamble == NPT);
    barostat.thermostat.U3_propagator(dt, &barostat.p_epsilon, 1);
}

void MDDriver::U4_NHC_thermo_propagator(double const dt, bool const outward) {
    assert(ensamble == NVT or ensamble == NPT);

    double e_couple = ekin_tot;

    if (thermostat_coupling != THERMOSTAT_KINETIC) {
        if (temperature_conf_inverse == 0.) {
            printf("WARNING: temperature_conf_inverse is zero; falling back to kinetic coupling\n");
        } else {
            double const econf = (Units::k_boltz / temperature_conf_inverse) * Nf / 2.;
            if (thermostat_coupling == THERMOSTAT_CONFIGURATIONAL) {
                e_couple = econf;
            } else { // THERMOSTAT_AVERAGE
                e_couple = 0.5 * (ekin_tot + econf);
            }
        }
    }

    thermostat.U4_propagator(dt, e_couple, outward);
}

void MDDriver::U4_NHC_baro_propagator(double const dt, bool const outward) {
    assert(ensamble == NPT);
    double const e_epsilon = 0.5 * pow<2>(barostat.p_epsilon) / barostat.w;
    barostat.thermostat.U4_propagator(dt, e_epsilon, outward);
}

void MDDriver::U_NHC_thermo_propagator(double const dt) {
    assert(ensamble == NVT or ensamble == NPT);

    // Suzuki-Yoshida method
    for (double const w : numerics::wsy7) {
        double const dt_w = w * dt;
        U4_NHC_thermo_propagator(0.5 * dt_w, false);
        U3_NHC_thermo_propagator(dt_w);
        U4_NHC_thermo_propagator(0.5 * dt_w, true);
    }
}

void MDDriver::U_NHC_baro_propagator(double const dt) {
    assert(ensamble == NPT);

    // Suzuki-Yoshida method
    for (double const w : numerics::wsy7) {
        double const dt_w = w * dt;
        U4_NHC_baro_propagator(0.5 * dt_w, false);
        U3_NHC_baro_propagator(dt_w);
        U4_NHC_baro_propagator(0.5 * dt_w, true);
    }
}

void MDDriver::integrate() {
    // scheme is based on:
    // https://github.com/Allen-Tildesley/examples/blob/master/md_nvt_lj.f90

    if (ensamble == NPT) {
        U_NHC_baro_propagator(0.5 * dt_global);
    }

    if (ensamble == NVT or ensamble == NPT) {
        U_NHC_thermo_propagator(0.5 * dt_global);
    }

    if (ensamble == NPT) {
        U2p_propagator(0.5 * dt_global);
    }

    // Velocity Verlet algorithm
    U2_propagator(0.5 * dt_global);
    U1_propagator(dt_global);
    calculate_forces();
    U2_propagator(0.5 * dt_global);

    if (ensamble == NPT) {
        U2p_propagator(0.5 * dt_global);
    }

    if (ensamble == NVT or ensamble == NPT) {
        U_NHC_thermo_propagator(0.5 * dt_global);
    }

    if (ensamble == NPT) {
        U_NHC_baro_propagator(0.5 * dt_global);
    }

    calculate_thermodynamics();

    time += dt_global;
    ++cycle;
}

void MDDriver::set_kinetic_energies() {
    for (unsigned int i = 0; i < Nparticles; ++i) {
        energies_kin[i] = 0.5 * mass[types[i]] * velocities[i].length2();
    }
    ekin_tot = std::accumulate(energies_kin.begin(), energies_kin.end(), 0.);
    ekin_atom = ekin_tot / Nparticles;
}

void MDDriver::set_momentum_tot() {

    // Must reset each call; otherwise it accumulates across timesteps (wrong and will eventually
    // trip the check).
    momentum_tot = Vector3d(0., 0., 0.);
    for (unsigned int i = 0; i < Nparticles; ++i) {
        momentum_tot += velocities[i] * mass[types[i]];
    }

    double const ptot_scale_2 =
        2. * mass[0] * ekin_tot; // I chose mass[0] as i assume that the mass of all the particles
                                 // are of the same magnitude
    if (momentum_tot.length2() / ptot_scale_2 > 1e-14) {
        printf("FATAL: total momentum = (%g, %g, %g) scale=%g\n", momentum_tot.x, momentum_tot.y,
               momentum_tot.z, std::sqrt(ptot_scale_2));
        exit(1);
    }
}

void MDDriver::calculate_configurational_temperature() {
    // The Laplacian is current only on sampled steps; otherwise keep the last sampled T_conf.
    if (not sample_observables)
        return;

    total_force2 = 0.;

    for (auto const &f : forces) {
        total_force2 += f.length2();
    }

    if (total_force2 != 0. and total_laplacian != 0.) {
        temperature_conf_inverse = Units::k_boltz * total_laplacian / total_force2;
    } else {
        temperature_conf_inverse = 0.;
    }
}

void MDDriver::calculate_thermodynamics() {

    set_momentum_tot();

    epot_atom = std::accumulate(energies_pot.begin(), energies_pot.end(), 0.) / Nparticles;
    epot_eam_atom = std::accumulate(F_eam.begin(), F_eam.end(), 0.) / Nparticles;
    etot_atom = ekin_atom + epot_atom;

    double const etot_cons_old = etot_cons;
    etot_cons = etot_atom * Nparticles;

    //----- Eq. 9 in: T.-Q. Yu et al. / Chemical Physics 370 (2010) 294305
    if (ensamble == NVT or ensamble == NPT) {
        etot_cons += thermostat.calculate_etot();
    }

    if (ensamble == NPT) {
        etot_cons += box.volume() * barostat.pressure + barostat.calculate_etot();
    }

    etot_cons_err = std::abs(etot_cons / etot_cons_initial - 1.);
    double const etot_cons_err_current = std::abs(etot_cons / etot_cons_old - 1.);

    temperature_kin = 2. * ekin_tot / (Units::k_boltz * Nf);
    calculate_configurational_temperature();

    pressure = (ekin_tot * 2. / 3. + virial) / box.volume();

    stress_tensor = -born_stress_tensor;
    for (unsigned int i = 0; i < Nparticles; ++i) {
        Vector3d const v = velocities[i];
        double const m = mass[types[i]];

        for (unsigned int iv = 0; iv < 6; ++iv) {
            unsigned int const ii = voigt_first_index[iv];
            unsigned int const jj = voigt_second_index[iv];
            stress_tensor(iv) += m * v[ii] * v[jj];
        }
    }
    stress_tensor /= box.volume();
    born_stress_tensor /= box.volume();
    if (sample_observables) {
        born_elasticity_tensor /= box.volume();
    }

    double const p_trace = stress_tensor[0] + stress_tensor[1] + stress_tensor[2];
    double p_err = std::abs(pressure / (p_trace / 3.) - 1.);
    if (p_err > 1e-5) {
        printf("FATAL PRESS err=%g %1.16E %1.16E\n", p_err, pressure, p_trace / 3.);
        exit(1);
    }

    bool const bad_conservation =
        (etot_cons_err > fatal_conservation or etot_cons_err_current > fatal_conservation) and
        cycle > 0;
    if (cycle % 100 == 0 or bad_conservation) {

        printf("%lld t=%1.3E dt=%g: "
               "e=%g ep=%g eEMB=%g ek=%g[eV] "
               "er0=%1.1E erc=%1.1E "
               "Tk=%g[K] "
               "Tc=%g[K] "
               "Tth=%g[K] "
               "V/V0=%g "
               "L=%g[A] "
               "P=%g "
               "Pb=%g[GPa]\n",
               cycle, time, dt_global, etot_atom / Units::ev, epot_atom / Units::ev,
               epot_eam_atom / Units::ev, ekin_atom / Units::ev, etot_cons_err,
               etot_cons_err_current, temperature_kin, 1. / temperature_conf_inverse,
               thermostat.temperature, box.volume() / box_initial.volume(),
               box.length() / Units::angstrom, pressure / Units::GPa,
               barostat.pressure / Units::GPa);

        if (bad_conservation) {
            printf("FATAL: bad energy conservation...\n");
            exit(1);
        }
    }

    // Pair-correlation (g(r)) is expensive (O(N^2)) and is better sampled explicitly from Python
    // at the desired cadence (e.g. only after thermalization). See python `molecular_dynamics.py`.
}

void MDDriver::calculate_pair_correlation() {
    std::fill(gr.begin(), gr.end(), 0.);

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto const &ci = coordinates[i];
        for (unsigned int j = i + 1; j < Nparticles; ++j) {
            Vector3d const dp = get_position(ci - coordinates[j]);
            double const dr = dp.length();

            if (dr < rmaxg) {
                unsigned int const k = std::floor(dr / drg);
                assert(k < gr.size());
                gr[k] += 1.;
            }
        }
    }

    for (unsigned int k = 0; k < rg.size(); ++k) {
        gr[k] *= box.volume() / (2. * M_PI * Nparticles * Nparticles * rg[k] * rg[k] * drg);
    }
}

void MDDriver::calculate_pair_correlations() {
    std::fill(gr.begin(), gr.end(), 0.);
    for (auto &g : gr_partial) {
        std::fill(g.begin(), g.end(), 0.);
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto const &ci = coordinates[i];
        unsigned int const type_i = types[i];
        for (unsigned int j = i + 1; j < Nparticles; ++j) {
            unsigned int const type_j = types[j];
            Vector3d const dp = get_position(ci - coordinates[j]);
            double const dr = dp.length();

            if (dr < rmaxg) {
                unsigned int const k = std::floor(dr / drg);
                assert(k < rg.size());
                gr[k] += 1.;

                unsigned int const tmax = std::max(type_i, type_j);
                unsigned int const tmin = std::min(type_i, type_j);
                unsigned int const pair_index = tmax * (tmax + 1) / 2 + tmin;
                assert(pair_index < gr_partial.size());
                gr_partial[pair_index][k] += 1.;
            }
        }
    }

    for (unsigned int k = 0; k < rg.size(); ++k) {
        gr[k] *= box.volume() / (2. * M_PI * Nparticles * Nparticles * rg[k] * rg[k] * drg);
    }

    unsigned int pair_index = 0;
    for (unsigned int b = 0; b < Ntypes; ++b) {
        for (unsigned int a = 0; a <= b; ++a) {
            double const Na = static_cast<double>(Nparticles_per_type[a]);
            double const Nb = static_cast<double>(Nparticles_per_type[b]);
            double pair_count = 0.;
            if (a == b) {
                if (Na > 1) {
                    pair_count = Na * (Na - 1.) / 2.;
                }
            } else {
                pair_count = Na * Nb;
            }

            if (pair_count > 0.) {
                double const prefactor = box.volume() / (4. * M_PI * pair_count * drg);
                for (unsigned int k = 0; k < rg.size(); ++k) {
                    gr_partial[pair_index][k] *= prefactor / (rg[k] * rg[k]);
                }
            } else {
                std::fill(gr_partial[pair_index].begin(), gr_partial[pair_index].end(), 0.);
            }
            pair_index += 1;
        }
    }
}

void MDDriver::calculate_partial_pair_correlation() {
    for (auto &g : gr_partial) {
        std::fill(g.begin(), g.end(), 0.);
    }

    for (unsigned int i = 0; i < Nparticles; ++i) {
        auto const &ci = coordinates[i];
        unsigned int const type_i = types[i];
        for (unsigned int j = i + 1; j < Nparticles; ++j) {
            unsigned int const type_j = types[j];
            Vector3d const dp = get_position(ci - coordinates[j]);
            double const dr = dp.length();

            if (dr < rmaxg) {
                unsigned int const k = std::floor(dr / drg);
                assert(k < rg.size());
                unsigned int const tmax = std::max(type_i, type_j);
                unsigned int const tmin = std::min(type_i, type_j);
                unsigned int const pair_index = tmax * (tmax + 1) / 2 + tmin;
                assert(pair_index < gr_partial.size());
                gr_partial[pair_index][k] += 1.;
            }
        }
    }

    unsigned int pair_index = 0;
    for (unsigned int b = 0; b < Ntypes; ++b) {
        for (unsigned int a = 0; a <= b; ++a) {
            double const Na = static_cast<double>(Nparticles_per_type[a]);
            double const Nb = static_cast<double>(Nparticles_per_type[b]);
            double pair_count = 0.;
            if (a == b) {
                if (Na > 1) {
                    pair_count = Na * (Na - 1.) / 2.;
                }
            } else {
                pair_count = Na * Nb;
            }

            if (pair_count > 0.) {
                double const prefactor = box.volume() / (4. * M_PI * pair_count * drg);
                for (unsigned int k = 0; k < rg.size(); ++k) {
                    gr_partial[pair_index][k] *= prefactor / (rg[k] * rg[k]);
                }
            } else {
                std::fill(gr_partial[pair_index].begin(), gr_partial[pair_index].end(), 0.);
            }
            pair_index += 1;
        }
    }
}
