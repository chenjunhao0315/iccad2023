#ifndef BMATCHQBF_HPP
#define BMATCHQBF_HPP

#include "sat/cnf/cnf.h"
#include "aig/gia/gia.h"

struct Bmatch_Qbf_Man_t {
    Gia_Man_t *     pGia;           // original miter
    int             nPars;          // parameter variables
    int             nVars;          // functional variables
    int             fVerbose;       // verbose flag
    // internal variables
    int             iParVarBeg;     // SAT var ID of the first par variable in the ver solver
    CaDiCaL::Solver *pSatVer;       // verification instance
    CaDiCaL::Solver *pSatSyn;       // synthesis instance
    Vec_Int_t *     vValues;        // variable values
    Vec_Int_t *     vParMap;        // parameter mapping
    Vec_Int_t *     vLits;          // literals for the SAT solver
    abctime         clkStart;       // global timeout
    abctime         clkSat;         // SAT solver time
};

#endif // BMATCHQBF_HPP