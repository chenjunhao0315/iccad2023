#ifndef BMATCH_CADICALSAT_HPP
#define BMATCH_CADICALSAT_HPP

#include "cadical/cadical.hpp"

static inline int Bmatch_toLit(int lit) {
    return lit + 1;
}

static inline int Bmatch_toLitCond(int lit, int cond) {
    return (cond) ? -(lit + 1) : (lit + 1);
}

static inline CaDiCaL::Solver *Bmatch_sat_solver_new() {
    CaDiCaL::Solver *s = new CaDiCaL::Solver();
    s->set("quiet", 1);
    return s;
}

static inline void Bmatch_sat_solver_delete(CaDiCaL::Solver *solver) {
    delete solver;
}

static inline void Bmatch_sat_solver_setnvars(CaDiCaL::Solver *solver, int vars) {
    solver->reserve(vars);
}

static inline void Bmatch_sat_solver_addvar(CaDiCaL::Solver *solver) {
    solver->reserve(solver->vars() + 1);
}

static inline int Bmatch_sat_solver_nvars(CaDiCaL::Solver *solver) {
    return solver->vars();
}

static inline void Bmatch_sat_solver_addclause(CaDiCaL::Solver *solver, int *lit_begin, int *lit_end) {
    for ( ; lit_begin < lit_end; ++lit_begin) {
        solver->add(*lit_begin);
    }
    solver->add(0);
}

static inline int Bmatch_sat_solver_solve(CaDiCaL::Solver *solver, int *assum_begin, int *assum_end, int a, int b, int c, int d) {
    for ( ; assum_begin < assum_end; ++assum_begin) {
        solver->assume(*assum_begin);
    }
    return solver->solve();
}

static inline int Bmatch_sat_solver_simplify(CaDiCaL::Solver *solver, int rounds = 0) {
    return solver->simplify(rounds);
}

static inline int Bmatch_sat_solver_var_value(CaDiCaL::Solver *solver, int var) {
    return solver->val(Bmatch_toLit(var)) > 0;
}

static int *Bmatch_sat_solver_get_model(CaDiCaL::Solver *pSolver, int *pVars, int nVars) {
    int *pModel = (int*)malloc(sizeof(int) * (nVars + 1));
    for (int i = 0; i < nVars; i++)
        pModel[i] = Bmatch_sat_solver_var_value(pSolver, pVars[i]) > 0;
    return pModel;
}

static inline int Bmatch_sat_solver_add_buffer(CaDiCaL::Solver *pSolver, int iVarA, int iVarB, int fCompl) {
    int Lits[2];
    assert( iVarA >= 0 && iVarB >= 0 );

    Lits[0] = Bmatch_toLitCond( iVarA, 0 );
    Lits[1] = Bmatch_toLitCond( iVarB, !fCompl );
    Bmatch_sat_solver_addclause( pSolver, Lits, Lits + 2 );

    Lits[0] = Bmatch_toLitCond( iVarA, 1 );
    Lits[1] = Bmatch_toLitCond( iVarB, fCompl );
    Bmatch_sat_solver_addclause( pSolver, Lits, Lits + 2 );

    return 2;
}

#endif // BMATCH_CADICALSAT_HPP
