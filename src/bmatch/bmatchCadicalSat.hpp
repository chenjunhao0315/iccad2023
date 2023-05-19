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

static int Bmatch_sat_solver_simplify(CaDiCaL::Solver *solver, int rounds = 3) {
    return solver->simplify(rounds);
}

static int Bmatch_sat_solver_var_value(CaDiCaL::Solver *solver, int var) {
    return solver->val(Bmatch_toLit(var));
}

static int *Bmatch_sat_solver_get_model(CaDiCaL::Solver *pSolver, int *pVars, int nVars) {
    int *pModel = (int*)malloc(sizeof(int) * (nVars + 1));
    for (int i = 0; i < nVars; i++)
        pModel[i] = Bmatch_sat_solver_var_value(pSolver, pVars[i]);
    return pModel;
}

#endif // BMATCH_CADICALSAT_HPP
