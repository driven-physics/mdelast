#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <memory>
#include "md_driver.hpp"
#include "units.hpp"
#include "potential.hpp"
#include "eam.hpp"
#include "vector3d.hpp"
#include <Eigen/Dense>

#include <pybind11/stl_bind.h>

void bind_module(pybind11::module &m) {

    /***************************************
     * Vector3d
     **************************************/
    pybind11::class_<Vector3d, std::shared_ptr<Vector3d>>(m, "Vector3d")
        .def(pybind11::init<double const, double const, double const>())
        .def_readwrite("x", &Vector3d::x)
        .def_readwrite("y", &Vector3d::y)
        .def_readwrite("z", &Vector3d::z)
        .def("length", &Vector3d::length)
        .def("length2", &Vector3d::length2)
        .def("__repr__", [](const Vector3d &v) {
            return "(" + std::to_string(v.x) + "," + std::to_string(v.y) + "," +
                   std::to_string(v.z) + ")";
        });

    /***************************************
     * Potentials
     **************************************/
    pybind11::class_<potential::values, std::shared_ptr<potential::values>>(m, "PotentialValue")
        .def_readonly("V", &potential::values::V)
        .def_readonly("dVdr", &potential::values::dVdr)
        .def_readonly("d2Vdr2", &potential::values::d2Vdr2);

    pybind11::class_<potential::Potential, std::shared_ptr<potential::Potential>>(m, "Potential")
        .def_readonly("rcut", &potential::Potential::rcut)
        .def("calculate", &potential::Potential::calculate);

    pybind11::class_<potential::Morse, potential::Potential, std::shared_ptr<potential::Morse>>(
        m, "Morse")
        .def(pybind11::init<double const, double const, double const, double const>())
        .def_readonly("D", &potential::Morse::D)
        .def_readonly("r0", &potential::Morse::r0)
        .def_readonly("alpha", &potential::Morse::alpha);

    pybind11::class_<potential::LennardJones, potential::Potential,
                     std::shared_ptr<potential::LennardJones>>(m, "LennardJones")
        .def(pybind11::init<double const, double const, double const>())
        .def_readonly("epsilon", &potential::LennardJones::epsilon)
        .def_readonly("sigma", &potential::LennardJones::sigma);

    pybind11::class_<potential::TablePotential, potential::Potential,
                     std::shared_ptr<potential::TablePotential>>(m, "TablePotential")
        .def(pybind11::init<double const, double const, double const, std::vector<double> const &,
                            std::vector<double> const &, std::vector<double> const &>());

    pybind11::class_<potential::EAM, std::shared_ptr<potential::EAM>>(m, "EAM")
        .def(pybind11::init<double const, double const, double const, double const,
                            std::vector<double> const &, std::vector<double> const &,
                            std::vector<double> const &, std::vector<double> const &,
                            std::vector<double> const &, std::vector<double> const &>())
        .def("rho", &potential::EAM::rho)
        .def("drhodr", &potential::EAM::drhodr)
        .def("d2rhodr2", &potential::EAM::d2rhodr2)
        .def("F", &potential::EAM::F)
        .def("dFdrho", &potential::EAM::dFdrho)
        .def("d2Fdrho2", &potential::EAM::d2Fdrho2);

    pybind11::bind_vector<std::vector<std::shared_ptr<potential::Potential>>>(m, "PotentialVector");
    pybind11::bind_vector<std::vector<std::shared_ptr<potential::EAM>>>(m, "EAMVector");

    /***************************************
     * MD Driver
     **************************************/
    pybind11::class_<MDDriver, std::shared_ptr<MDDriver>> md_driver(m, "MDDriver");
    md_driver
        .def(
            pybind11::init<std::vector<Vector3d> const &, std::vector<Vector3d> const &,
                           std::vector<potential::PotentialPtr> const &,
                           std::vector<potential::EAMPtr> const &, MDDriver::Ensamble const,
                           MDDriver::ThermostatCoupling const, Thermostat const &, Barostat const &,
                           MDDriver::Box const &, std::vector<double> const &, double const,
                           std::vector<double> const &, std::vector<unsigned int> const &,
                           double const, bool const, double const, double const>())
        .def_readonly("positions", &MDDriver::positions)
        .def_readonly("coordinates", &MDDriver::coordinates)
        .def_readonly("velocities", &MDDriver::velocities)
        .def_readonly("forces", &MDDriver::forces)
        .def_readonly("rho_eam", &MDDriver::rho_eam)
        .def_readonly("energies_pot", &MDDriver::energies_pot)
        .def_readonly("energies_kin", &MDDriver::energies_kin)
        .def_readonly("ensamble", &MDDriver::ensamble)
        .def_readonly("thermostat_coupling", &MDDriver::thermostat_coupling)
        .def_readonly("thermostat", &MDDriver::thermostat)
        .def_readonly("barostat", &MDDriver::barostat)
        .def_readonly("box_initial", &MDDriver::box_initial)
        .def_readonly("box", &MDDriver::box)
        .def_readonly("epot_eam_atom", &MDDriver::epot_eam_atom)
        .def_readonly("epot_atom", &MDDriver::epot_atom)
        .def_readonly("ekin_atom", &MDDriver::ekin_atom)
        .def_readonly("etot_atom", &MDDriver::etot_atom)
        .def_readonly("etot_cons", &MDDriver::etot_cons)
        .def_readonly("etot_cons_err", &MDDriver::etot_cons_err)
        .def_readonly("born_stress_tensor", &MDDriver::born_stress_tensor)
        .def_readonly("born_elasticity_tensor", &MDDriver::born_elasticity_tensor)
        .def_readonly("temperature_kin", &MDDriver::temperature_kin)
        .def_readonly("temperature_conf_inverse", &MDDriver::temperature_conf_inverse)
        .def_readonly("total_laplacian", &MDDriver::total_laplacian)
        .def_readonly("total_force2", &MDDriver::total_force2)
        .def_readonly("pressure", &MDDriver::pressure)
        .def_readonly("stress_tensor", &MDDriver::stress_tensor)
        .def_readonly("dt_global", &MDDriver::dt_global)
        .def_readonly("time", &MDDriver::time)
        .def_readonly("cycle", &MDDriver::cycle)
        .def_readonly("rg", &MDDriver::rg)
        .def_readonly("gr", &MDDriver::gr)
        .def_readonly("gr_partial", &MDDriver::gr_partial)
        .def_readonly("fatal_conservation", &MDDriver::fatal_conservation)
        .def_readonly("mass", &MDDriver::mass)
        .def_readonly("types", &MDDriver::types)
        .def_readonly("Nparticles_per_type", &MDDriver::Nparticles_per_type)
        .def("integrate", &MDDriver::integrate)
        .def("calculate_neighbors", &MDDriver::calculate_neighbors)
        .def("calculate_pair_correlation", &MDDriver::calculate_pair_correlation)
        .def("calculate_pair_correlations", &MDDriver::calculate_pair_correlations)
        .def("calculate_partial_pair_correlation", &MDDriver::calculate_partial_pair_correlation);

    pybind11::class_<MDDriver::Box, std::shared_ptr<MDDriver::Box>>(m, "Box")
        .def(pybind11::init<Eigen::Matrix3d const &>())
        .def_readonly("h", &MDDriver::Box::h)
        .def("volume", &MDDriver::Box::volume)
        .def("length", &MDDriver::Box::length);

    pybind11::class_<Thermostat, std::shared_ptr<Thermostat>>(m, "Thermostat")
        .def(pybind11::init<double const, std::vector<double> const &, std::vector<double> const &,
                            double const>())
        .def_readonly("temperature", &Thermostat::temperature)
        .def_readonly("q", &Thermostat::q)
        .def_readonly("p_eta", &Thermostat::p_eta)
        .def_readonly("eta", &Thermostat::eta)
        .def_readonly("nchain", &Thermostat::nchain)
        .def_readonly("Nf", &Thermostat::Nf);

    pybind11::class_<Barostat, std::shared_ptr<Barostat>>(m, "Barostat")
        .def(pybind11::init<double const, double const, double const, Thermostat const &>())
        .def_readonly("pressure", &Barostat::pressure)
        .def_readonly("w", &Barostat::w)
        .def_readonly("epsilon", &Barostat::epsilon)
        .def_readonly("p_epsilon", &Barostat::p_epsilon)
        .def_readonly("thermostat", &Barostat::thermostat);

    pybind11::class_<BarostatFlex, std::shared_ptr<BarostatFlex>>(m, "BarostatFlex")
        .def(
            pybind11::init<double const, double const, Eigen::Matrix3d const, Thermostat const &>())
        .def_readonly("pressure", &BarostatFlex::pressure)
        .def_readonly("w", &BarostatFlex::w)
        .def_readonly("pg", &BarostatFlex::pg)
        .def_readonly("thermostat", &BarostatFlex::thermostat);

    pybind11::enum_<MDDriver::Ensamble>(md_driver, "Ensamble")
        .value("NVE", MDDriver::Ensamble::NVE)
        .value("NVT", MDDriver::Ensamble::NVT)
        .value("NPT", MDDriver::Ensamble::NPT)
        .value("NPT_FLEX", MDDriver::Ensamble::NPT_FLEX)
        .export_values();

    pybind11::enum_<MDDriver::ThermostatCoupling>(md_driver, "ThermostatCoupling")
        .value("THERMOSTAT_KINETIC", MDDriver::ThermostatCoupling::THERMOSTAT_KINETIC)
        .value("THERMOSTAT_CONFIGURATIONAL",
               MDDriver::ThermostatCoupling::THERMOSTAT_CONFIGURATIONAL)
        .value("THERMOSTAT_AVERAGE", MDDriver::ThermostatCoupling::THERMOSTAT_AVERAGE)
        .export_values();

    /***************************************
     * Units
     **************************************/
    pybind11::class_<Units::BindingClass>(m, "Units")
        .def_readonly_static("k_boltz", &Units::k_boltz)
        .def_readonly_static("Navogadro", &Units::Navogadro)
        .def_readonly_static("ev", &Units::ev)
        .def_readonly_static("joule", &Units::joule)
        .def_readonly_static("angstrom", &Units::angstrom)
        .def_readonly_static("nm", &Units::nm)
        .def_readonly_static("bar", &Units::bar)
        .def_readonly_static("kbar", &Units::kbar)
        .def_readonly_static("atm", &Units::atm)
        .def_readonly_static("Pa", &Units::Pa)
        .def_readonly_static("GPa", &Units::GPa)
        .def_readonly_static("mus", &Units::mus)
        .def_readonly_static("ns", &Units::ns)
        .def_readonly_static("ps", &Units::ps)
        .def_readonly_static("fs", &Units::fs);
}

PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<potential::Potential>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<potential::EAM>>)

/*
 * the name 'md_driver_module' MUST also be the name of this file ('md_driver_module.cpp').
 * this same name should appear in the CMakeLists.txt
 */
PYBIND11_MODULE(md_driver_module, m) {
    m.doc() = "md driver c++ module";
    bind_module(m);
}
