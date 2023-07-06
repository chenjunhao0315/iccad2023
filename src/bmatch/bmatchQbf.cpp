#include "bmatch.hpp"
#include "bmatchQbf.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

extern CaDiCaL::Solver *Bmatch_Cnf_DataWriteIntoSolver(CaDiCaL::Solver *pSolver, Cnf_Dat_t *p);
extern Gia_Man_t * Gia_QbfCofactor( Gia_Man_t * p, int nPars, Vec_Int_t * vValues, Vec_Int_t * vParMap );
extern void Cnf_SpecialDataLift( Cnf_Dat_t * p, int nVarsPlus, int firstPiVar, int lastPiVar);

Bmatch_Qbf_Man_t *Bmatch_Gia_QbfAlloc(Gia_Man_t *pGia, int nPars, int fVerbose);
int Bmatch_Gia_QbfVerify(Bmatch_Qbf_Man_t *p, Vec_Int_t *vValues);
void Bmatch_Gia_QbfOnePattern( Bmatch_Qbf_Man_t * p, Vec_Int_t * vValues );
void Bmatch_Gia_QbfFree( Bmatch_Qbf_Man_t * p );

int Bmatch_Gia_QbfSolveValue(Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);
int Bmatch_Gia_QbfSolveValueInt(Bmatch_Qbf_Man_t *p, Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);

#ifdef __cplusplus
}
#endif

Bmatch_Qbf_Man_t *Bmatch_Gia_QbfAlloc(Gia_Man_t *pGia, int nPars, int fVerbose) {
    Bmatch_Qbf_Man_t *p;
    Cnf_Dat_t * pCnf;
    Gia_ObjFlipFaninC0( Gia_ManPo(pGia, 0) );
    pCnf = (Cnf_Dat_t *)Mf_ManGenerateCnf( pGia, 8, 0, 1, 0, 0 );
    Gia_ObjFlipFaninC0( Gia_ManPo(pGia, 0) );
    p = ABC_CALLOC( Bmatch_Qbf_Man_t, 1 );
    p->clkStart   = Abc_Clock();
    p->pGia       = pGia;
    p->nPars      = nPars;
    p->nVars      = Gia_ManPiNum(pGia) - nPars;
    p->fVerbose   = fVerbose;
    p->iParVarBeg = pCnf->nVars - Gia_ManPiNum(pGia);// - 1;
    p->pSatVer    = Bmatch_Cnf_DataWriteIntoSolver(Bmatch_sat_solver_new(), pCnf);
    p->pSatSyn    = Bmatch_sat_solver_new();
    p->vValues    = Vec_IntAlloc( Gia_ManPiNum(pGia) );
    p->vParMap    = Vec_IntStartFull( nPars );
    p->vLits      = Vec_IntAlloc( nPars );
    Bmatch_sat_solver_setnvars( p->pSatSyn, nPars );
    Cnf_DataFree( pCnf );
    return p;
}

void Bmatch_Gia_QbfFree( Bmatch_Qbf_Man_t * p ) {
    Bmatch_sat_solver_delete( p->pSatVer );
    Bmatch_sat_solver_delete( p->pSatSyn );
    Vec_IntFree( p->vLits );
    Vec_IntFree( p->vValues );
    Vec_IntFree( p->vParMap );
    ABC_FREE( p );
}

int Bmatch_Gia_QbfVerify(Bmatch_Qbf_Man_t *p, Vec_Int_t *vValues) {
    int i, Entry, RetValue;
    assert( Vec_IntSize(vValues) == p->nPars );
    Vec_IntClear( p->vLits );
    Vec_IntForEachEntry( vValues, Entry, i )
        Vec_IntPush( p->vLits, Bmatch_toLitCond(p->iParVarBeg+i, !Entry) );
    RetValue = Bmatch_sat_solver_solve( p->pSatVer, Vec_IntArray(p->vLits), Vec_IntLimit(p->vLits), 0, 0, 0, 0 );
    if ( RetValue == 10 )
    {
        Vec_IntClear( vValues );
        for ( i = 0; i < p->nVars; i++ )
            Vec_IntPush( vValues, Bmatch_sat_solver_var_value(p->pSatVer, p->iParVarBeg+p->nPars+i) );
    }
    return RetValue == 10 ? 1 : 0;
}

int Gia_QbfAddCofactor( Bmatch_Qbf_Man_t * p, Gia_Man_t * pCof ) {
    Cnf_Dat_t * pCnf = (Cnf_Dat_t *)Mf_ManGenerateCnf( pCof, 8, 0, 1, 0, 0 );
    int i, useold = 0;
    int iFirstVar = useold ? Bmatch_sat_solver_nvars(p->pSatSyn) + pCnf->nVars - Gia_ManPiNum(pCof) : pCnf->nVars - Gia_ManPiNum(pCof); //-1   
    pCnf->pMan = NULL;
    
    if (useold)    
        Cnf_DataLift( pCnf, Bmatch_sat_solver_nvars(p->pSatSyn) );
    else
        Cnf_SpecialDataLift( pCnf, Bmatch_sat_solver_nvars(p->pSatSyn), iFirstVar, iFirstVar + Gia_ManPiNum(p->pGia) );

    for ( i = 0; i < pCnf->nClauses; i++ ) {
        AutoBuffer<int> pLits(pCnf->pClauses[i + 1] - pCnf->pClauses[i]);
        for (int j = 0, *k = pCnf->pClauses[i]; j < pLits.size(); ++j, ++k) {
            pLits[j] = Bmatch_toLitCond(((*k) >> 1), (*k) & 1);
        }
        Bmatch_sat_solver_addclause( p->pSatSyn, pLits, pLits + pLits.size() );
    }
    Cnf_DataFree( pCnf );
    // add connection clauses
    if (useold)
           for ( i = 0; i < Gia_ManPiNum(p->pGia); i++ )
            if ( !Bmatch_sat_solver_add_buffer( p->pSatSyn, i, iFirstVar+i, 0 ) )
                return 0;
    return 1;
}

void Bmatch_Gia_QbfOnePattern( Bmatch_Qbf_Man_t * p, Vec_Int_t * vValues ) {
    int i;
    Vec_IntClear( vValues );
    for ( i = 0; i < p->nPars; i++ )
        Vec_IntPush( vValues, Bmatch_sat_solver_var_value(p->pSatSyn, i) );
}

int Bmatch_Gia_QbfSolveValue(Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose) {
    Bmatch_Qbf_Man_t *p = Bmatch_Gia_QbfAlloc( pGia, nPars, fVerbose );
    int RetValue = Bmatch_Gia_QbfSolveValueInt(p, pGia, vValues, nPars, nIterLimit, nConfLimit, nTimeOut, fGlucose, fVerbose);
    Bmatch_Gia_QbfFree( p );
    return RetValue;
}

int Bmatch_Gia_QbfSolveValueInt(Bmatch_Qbf_Man_t *p, Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose) {
    Gia_Man_t *pCof;
    int i, status, RetValue = 0;
    abctime clk;
    if ( fVerbose )
        printf( "Solving QBF for \"%s\" with %d parameters, %d variables and %d AIG nodes.\n", Gia_ManName(pGia), p->nPars, p->nVars, Gia_ManAndNum(pGia) );
    if ( p->pSatVer == 0x0 )
        return 1;

    assert( Gia_ManRegNum(pGia) == 0 );
    Vec_IntFill( p->vValues, nPars, 0 );
    Vec_IntFill( vValues, nPars, 0 );
    for ( i = 0; Bmatch_Gia_QbfVerify(p, p->vValues); i++ ) {
        // generate next constraint
        assert( Vec_IntSize(p->vValues) == p->nVars );
        pCof = Gia_QbfCofactor( pGia, nPars, p->vValues, p->vParMap );
        status = Gia_QbfAddCofactor( p, pCof );
        Gia_ManStop( pCof );
        if ( status == 0 )       { RetValue =  1; break; }
        // synthesize next assignment
        clk = Abc_Clock();
        status = Bmatch_sat_solver_solve( p->pSatSyn, NULL, NULL, (ABC_INT64_T)nConfLimit, 0, 0, 0 );
        p->clkSat += Abc_Clock() - clk;
        if ( status == 20 ) { RetValue =  1; break; }
        if ( status == 0 ) { RetValue = -1; break; }
        // extract SAT assignment
        Bmatch_Gia_QbfOnePattern( p, p->vValues );
        Bmatch_Gia_QbfOnePattern( p, vValues );
        assert( Vec_IntSize(p->vValues) == p->nPars );
        if ( nIterLimit && i+1 == nIterLimit ) { RetValue = -1; break; }
        if ( nTimeOut && (Abc_Clock() - p->clkStart)/CLOCKS_PER_SEC >= nTimeOut ) { RetValue = -1; break; }
    }
    if ( RetValue == 0 && fVerbose) {
        int nZeros = Vec_IntCountZero( p->vValues );
        printf( "Parameters: " );
        assert( Vec_IntSize(p->vValues) == nPars );
        Vec_IntPrintBinary( p->vValues );
        printf( "  Statistics: 0=%d 1=%d\n", nZeros, Vec_IntSize(p->vValues) - nZeros );
        
    }
    if ( fVerbose ) {
        if ( RetValue == -1 && nTimeOut && (Abc_Clock() - p->clkStart)/CLOCKS_PER_SEC >= nTimeOut )
            printf( "The problem timed out after %d sec.  ", nTimeOut );
        else if ( RetValue == -1 && nConfLimit )
            printf( "The problem aborted after %d conflicts.  ", nConfLimit );
        else if ( RetValue == -1 && nIterLimit )
            printf( "The problem aborted after %d iterations.  ", nIterLimit );
        else if ( RetValue == 1 )
            printf( "The problem is UNSAT after %d iterations.  ", i );
        else 
            printf( "The problem is SAT after %d iterations.  ", i );
        printf( "\n" );
        Abc_PrintTime( 1, "SAT  ", p->clkSat );
        Abc_PrintTime( 1, "Other", Abc_Clock() - p->clkStart - p->clkSat );
        Abc_PrintTime( 1, "TOTAL", Abc_Clock() - p->clkStart );
    }
    return RetValue;
}

ABC_NAMESPACE_IMPL_END