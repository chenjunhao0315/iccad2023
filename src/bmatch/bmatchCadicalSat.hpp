#ifndef BMATCH_CADICALSAT_HPP
#define BMATCH_CADICALSAT_HPP

#include "cadical/cadical.hpp"

static inline int Bmatch_toLit(int lit) {
    return lit + 1;
}

static inline int Bmatch_toLitCond(int lit, int cond) {
    return (cond) ? -(lit + 1) : (lit + 1);
}

static CaDiCaL::Solver *Bmatch_sat_solver_new() {
    return new CaDiCaL::Solver();
}

static void Bmatch_sat_solver_delete(CaDiCaL::Solver *solver) {
    delete solver;
}

static void Bmatch_sat_solver_setnvars(CaDiCaL::Solver *solver, int vars) {
    solver->reserve(vars);
}

static void Bmatch_sat_solver_addvar(CaDiCaL::Solver *solver) {
    solver->reserve(solver->vars() + 1);
}

static int Bmatch_sat_solver_nvars(CaDiCaL::Solver *solver) {
    return solver->vars();
}

static void Bmatch_sat_solver_addclause(CaDiCaL::Solver *solver, int *lit_begin, int *lit_end) {
    for (int *i = lit_begin; i < lit_end; ++i) {
        solver->add(*i);
    }
    solver->add(0);
}

static int Bmatch_sat_solver_solve(CaDiCaL::Solver *solver, int *assum_begin, int *assum_end, int a, int b, int c, int d) {
    for (int *i = assum_begin; i < assum_end; ++i) {
        solver->assume(*i);
    }
    return solver->solve();
}

static int Bmatch_sat_solver_val(CaDiCaL::Solver *solver, int lit) {
    return solver->val(Bmatch_toLit(lit));
}



#endif // BMATCH_CADICALSAT_HPP
