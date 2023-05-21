#include "bmatch.hpp"

#include "proof/fra/fra.h"
#include "aig/aig/aig.h"
#include "proof/fraig/fraig.h"
#include "sat/cnf/cnf.h"
#include "base/io/ioAbc.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

extern Aig_Man_t *Abc_NtkToDar( Abc_Ntk_t *pNtk, int fExors, int fRegisters );
static void Bmatch_NtkVerifyReportError(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel, vMatch &MI, vMatch &MO);
int Bmatch_SatFraig(Abc_Ntk_t **ppNtk, int cadicalSat);
int Bmatch_FraigCadicalSat(Aig_Man_t *pMan, int fVerbose);
EcResult Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int cadicalSat, int fVerbose);

int Bmatch_Cnf_DataWriteOrClause(CaDiCaL::Solver *pSolver, Cnf_Dat_t *pCnf);
CaDiCaL::Solver *Bmatch_Cnf_DataWriteIntoSolver(CaDiCaL::Solver *pSolver, Cnf_Dat_t *p, int offset = 0);
CaDiCaL::Solver *Bmatch_ControllableInputOutputSat(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int &controlPiOffset);
CaDiCaL::Solver *Bmatch_ControllableInputSat(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, int &controlPiOffset);
EcResult Bmatch_NtkControllableInputEcFraig(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI);
EcResult Bmatch_NtkControllableInputOutputEcFraig(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);

#ifdef __cplusplus
}
#endif

EcResult Bmatch_NtkControllableInputOutputEcFraig(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO) {
    auto &pMiterSolver = pMan->pMiterSolverNew;
    auto controlOffset = pMan->controlOffset;

    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1))); // nPi + const (not include inv)
    int nControlPo = (int)(std::ceil(std::log2(Abc_NtkPoNum(pNtk1) + 1))); // nPo + nonmap (not include inv)
    AutoBuffer<int> ioControl((nControlPi + 1) * Abc_NtkPiNum(pNtk2) + (nControlPo + 1) * Abc_NtkPoNum(pNtk2), 0);

    int count = 0;
    for (int xi = 0; xi < MI.size(); ++xi) {
        for (auto y : MI[xi]) {
            int x = xi;
            int yi = y.var();
            for (int k = 0; k < nControlPi; ++k, x >>= 1) {
                ioControl[count++] = Bmatch_toLitCond(yi * (nControlPi + 1) + k + controlOffset, !(x & 1));
            }
            ioControl[count++] = Bmatch_toLitCond(yi * (nControlPi + 1) + nControlPi + controlOffset, !(y.sign()));
        }
    }

    int controlPiOffset = (nControlPi + 1) * Abc_NtkPiNum(pNtk2) + controlOffset;
    AutoBuffer<int> mappedPo(Abc_NtkPoNum(pNtk2), 0);
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto g : MO[fi]) {
            int f = fi;
            int gi = g.var();
            mappedPo[gi] = 1;
            for (int k = 0; k < nControlPo; ++k, f >>= 1) {
                ioControl[count++] = Bmatch_toLitCond(gi * (nControlPo + 1) + k + controlPiOffset, !(f & 1));
            }
            ioControl[count++] = Bmatch_toLitCond(gi * (nControlPo + 1) + nControlPo + controlPiOffset, !(g.sign()));
        }
    }
    for (int gi = 0; gi < mappedPo.size(); ++gi) {
        if (!mappedPo[gi]) {
            int f = Abc_NtkPoNum(pNtk1);
            for (int k = 0; k < nControlPo; ++k, f >>= 1) {
                ioControl[count++] = Bmatch_toLitCond(gi * (nControlPo + 1) + k + controlPiOffset, !(f & 1));
            }
            ioControl[count++] = Bmatch_toLitCond(gi * (nControlPo + 1) + nControlPo + controlPiOffset, 0);
        }
    }

    assert(count == ioControl.size());
    
    int status = Bmatch_sat_solver_solve(pMiterSolver, ioControl, ioControl + ioControl.size(), 0, 0, 0, 0);
    
    if (status == 10) {
        int *pModel1 = ABC_ALLOC(int, Abc_NtkPiNum(pNtk1));
        for (int i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
            pModel1[i] = Bmatch_sat_solver_var_value(pMiterSolver, (nControlPi + 1) * Abc_NtkPiNum(pNtk2) + (nControlPo + 1) * Abc_NtkPoNum(pNtk2) + controlOffset + i);
        }
        return {NON_EQUIVALENT, pModel1};
    }

    return {EQUIVALENT, NULL};
}

EcResult Bmatch_NtkControllableInputEcFraig(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI) {
    auto &pMiterSolver = pMan->pMiterSolver;
    auto &controlOffset = pMan->controlOffset;

    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));
    AutoBuffer<int> inputControl((nControlPi + 1) * Abc_NtkPiNum(pNtk2), 0);

    int count = 0;
    for (int xi = 0; xi < MI.size(); ++xi) {
        for (auto y : MI[xi]) {
            int x = xi;
            int yi = y.var();
            for (int k = 0; k < nControlPi; ++k, x >>= 1) {
                inputControl[count++] = Bmatch_toLitCond(yi * (nControlPi + 1) + k + controlOffset, !(x & 1));
            }
            inputControl[count++] = Bmatch_toLitCond(yi * (nControlPi + 1) + nControlPi + controlOffset, !(y.sign()));
        }
    }
    
    int status = Bmatch_sat_solver_solve(pMiterSolver, inputControl, inputControl + inputControl.size(), 0, 0, 0, 0);
    
    if (status == 10) {
        int *pModel1 = ABC_ALLOC(int, Abc_NtkPiNum(pNtk1));
        for (int i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
            pModel1[i] = Bmatch_sat_solver_var_value(pMiterSolver, (nControlPi + 1) * Abc_NtkPiNum(pNtk2) + controlOffset + i);
        }
        return {NON_EQUIVALENT, pModel1};
    }

    return {EQUIVALENT, NULL};
}

CaDiCaL::Solver *Bmatch_ConvertNtk2Sat(Abc_Ntk_t *pNtk, int &controlPiOffset) {
    Aig_Man_t *pMan = Abc_NtkToDar(pNtk, 0, 0);
    Cnf_Dat_t *pCnf;
    CaDiCaL::Solver *pSat;

    assert(Aig_ManRegNum(pMan) == 0);
    pMan->pData = NULL;

    // derive CNF
    pCnf = Cnf_Derive(pMan, Aig_ManCoNum(pMan));
    controlPiOffset = pCnf->pVarNums[Abc_NtkPi(pNtk, 0)->Id];

    // FlipBits
    Cnf_DataTranformPolarity(pCnf, 0);

    // write CNF
    pSat = Bmatch_Cnf_DataWriteIntoSolver(Bmatch_sat_solver_new(), pCnf, 0);
    if (pSat == NULL) {
        Cnf_DataFree( pCnf );
        return NULL;
    }
    Bmatch_Cnf_DataWriteOrClause(pSat, pCnf);

    Cnf_DataFree(pCnf);

    int status = Bmatch_sat_solver_simplify(pSat);

    if (status == 20) {
        Bmatch_sat_solver_delete(pSat);
        return NULL;
    }

    return pSat;
}

CaDiCaL::Solver *Bmatch_ControllableInputSat(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, int &controlPiOffset) {
    return Bmatch_ConvertNtk2Sat(Bmatch_NtkControllableInputMiter(pNtk1, pNtk2, MO), controlPiOffset);
}

CaDiCaL::Solver *Bmatch_ControllableInputOutputSat(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int &controlPiOffset) {
    return Bmatch_ConvertNtk2Sat(Bmatch_NtkControllableInputOutputMiter(pNtk1, pNtk2), controlPiOffset);
}

EcResult Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int cadicalSat, int fVerbose) {
    abctime clk = Abc_Clock();
    Abc_Ntk_t *pNtkMiter;
    int RetValue, Status;
    int *model = NULL;

    pNtkMiter = Bmatch_NtkMiter(pNtk1, pNtk2, MI, MO);

    if (pNtkMiter == NULL) {
        if (fVerbose) Abc_Print(1, "Miter computation has failed.\n");
        return {MITER_FAIL, NULL};
    }

    RetValue = Abc_NtkMiterIsConstant(pNtkMiter);
    if (RetValue == 0) {
        if (fVerbose) Abc_Print(1, "Networks are NOT EQUIVALENT after structural hashing.\n");
        // report the error
        pNtkMiter->pModel = Abc_NtkVerifyGetCleanModel(pNtkMiter, 1);
        if (fVerbose) Bmatch_NtkVerifyReportError(pNtk1, pNtk2, pNtkMiter->pModel, MI, MO);
        // ABC_FREE(pNtkMiter->pModel);
        model = pNtkMiter->pModel; pNtkMiter->pModel = NULL;
        if (fVerbose) Abc_PrintTime(1, "Time", Abc_Clock() - clk);
        Abc_NtkDelete(pNtkMiter);
        return {NON_EQUIVALENT, model};
    } else if (RetValue == 1) {
        if (fVerbose) Abc_Print(1, "Networks are equivalent after structural hashing.\n");
        if (fVerbose) Abc_PrintTime(1, "Time", Abc_Clock() - clk);
        Abc_NtkDelete(pNtkMiter);
        return {EQUIVALENT, NULL};
    }

    RetValue = Bmatch_SatFraig(&pNtkMiter, cadicalSat);
    if (RetValue == -1) {
        if (fVerbose) Abc_Print(1, "Networks are undecided (resource limits is reached).\n");
        Status = RESOURCE_LIMIT;
    } else if (RetValue == 0) {
        int *pSimInfo = Abc_NtkVerifySimulatePattern(pNtkMiter, pNtkMiter->pModel);
        if (pSimInfo[0] != 1) {
            if (fVerbose) Abc_Print(1, "ERROR in Abc_NtkMiterProve(): Generated counter-example is invalid.\n");
            Status = PROVE_ERROR;
        } else {
            if (fVerbose) Abc_Print(1, "Networks are NOT EQUIVALENT.\n");
            Status = NON_EQUIVALENT;
        }
        ABC_FREE(pSimInfo);
    } else {
        if (fVerbose) Abc_Print(1, "Networks are equivalent.\n");
        Status = EQUIVALENT;
    }
    if (fVerbose) Abc_PrintTime(1, "Time", Abc_Clock() - clk);
    if (pNtkMiter->pModel) {
        if (fVerbose) Bmatch_NtkVerifyReportError(pNtk1, pNtk2, pNtkMiter->pModel, MI, MO);
    }
    model = pNtkMiter->pModel; pNtkMiter->pModel = NULL;
    Abc_NtkDelete(pNtkMiter);

    return {Status, model};
}

int Bmatch_SatFraig(Abc_Ntk_t **ppNtk, int cadicalSat) {
    Abc_Ntk_t *pNtk = *ppNtk;
    Abc_Obj_t *pObj, *pFanin;
    Aig_Man_t *pMan;
    int RetValue;

    pObj = Abc_NtkPo(pNtk, 0);
    pFanin = Abc_ObjFanin0(pObj);
    if (Abc_ObjFanin0(pObj)->fPhase != (unsigned)Abc_ObjFaninC0(pObj)) {
        pNtk->pModel = ABC_CALLOC(int, Abc_NtkCiNum(pNtk));
        return 0;
    }

    pMan = Abc_NtkToDar(pNtk, 0, 0);
    RetValue = (cadicalSat) ? Bmatch_FraigCadicalSat(pMan, 0) : Fra_FraigSat(pMan, (ABC_INT64_T)5000, (ABC_INT64_T)0, 0, 0, 0, 1, 0, 0, 0);
    pNtk->pModel = (int *)pMan->pData, pMan->pData = NULL;
    Aig_ManStop(pMan);

    return RetValue;
}

CaDiCaL::Solver *Bmatch_Cnf_DataWriteIntoSolver(CaDiCaL::Solver *pSolver, Cnf_Dat_t *p, int offset) {
    CaDiCaL::Solver *pSat = pSolver;
    int i, f, status;
    assert(pSat);

    Bmatch_sat_solver_setnvars(pSat, p->nVars + offset);
    for (i = 0; i < p->nClauses; i++) {
        AutoBuffer<int> pLits(p->pClauses[i + 1] - p->pClauses[i]);
        for (int j = 0, *k = p->pClauses[i]; j < pLits.size(); ++j, ++k) {
            pLits[j] = Bmatch_toLitCond(((*k) >> 1) + offset, (*k) & 1);
        }
        Bmatch_sat_solver_addclause(pSat, pLits, pLits + pLits.size());
    }
    status = Bmatch_sat_solver_simplify(pSat);
    if (status == 20) {
        Bmatch_sat_solver_delete(pSat);
        return NULL;
    }
    return pSat;
}

int Bmatch_Cnf_DataWriteOrClause(CaDiCaL::Solver *pSolver, Cnf_Dat_t *pCnf) {
    CaDiCaL::Solver *pSat = pSolver;
    Aig_Obj_t *pObj;
    int i;
    AutoBuffer<int> pLits(Aig_ManCoNum(pCnf->pMan));
    Aig_ManForEachCo(pCnf->pMan, pObj, i)
        pLits[i] = Bmatch_toLitCond(pCnf->pVarNums[pObj->Id], 0);
    Bmatch_sat_solver_addclause( pSat, pLits, pLits + Aig_ManCoNum(pCnf->pMan));

    return 1;
}

int Bmatch_FraigCadicalSat(Aig_Man_t *pMan, int fVerbose) {
    CaDiCaL::Solver *pSat;
    Cnf_Dat_t *pCnf;
    int status, RetValue = 0;
    Vec_Int_t *vCiIds;
    abctime clk = Abc_Clock();

    assert(Aig_ManRegNum(pMan) == 0);
    pMan->pData = NULL;

    // derive CNF
    pCnf = Cnf_Derive(pMan, Aig_ManCoNum(pMan));

    // FlipBits
    Cnf_DataTranformPolarity(pCnf, 0);

    if (fVerbose) {
        printf( "CNF stats: Vars = %6d. Clauses = %7d. Literals = %8d. ", pCnf->nVars, pCnf->nClauses, pCnf->nLiterals );
        Abc_PrintTime(1, "Time", Abc_Clock() - clk);
    }

    // write CNF
    pSat = Bmatch_Cnf_DataWriteIntoSolver(Bmatch_sat_solver_new(), pCnf);
    if (pSat == NULL) {
        Cnf_DataFree( pCnf );
        return 1;
    }
    Bmatch_Cnf_DataWriteOrClause(pSat, pCnf);

    vCiIds = Cnf_DataCollectPiSatNums(pCnf, pMan);
    Cnf_DataFree(pCnf);

    clk = Abc_Clock();
    status = Bmatch_sat_solver_simplify(pSat);

    if (status == 20) {
        Vec_IntFree( vCiIds );
        Bmatch_sat_solver_delete( pSat );
        return 1;
    }

    clk = Abc_Clock();

    status = Bmatch_sat_solver_solve(pSat, NULL, NULL, 0, 0, 0, 0);
    if (status == 0) {
        RetValue = -1;
    } else if (status == 10) {
        RetValue = 0;
    } else if (status == 20) {
        RetValue = 1;
    } else
        assert(false);

    // if the problem is SAT, get the counterexample
    if (status == 10) {
        pMan->pData = Bmatch_sat_solver_get_model(pSat, vCiIds->pArray, vCiIds->nSize);
    }

    Bmatch_sat_solver_delete(pSat);
    Vec_IntFree( vCiIds );
    return RetValue;
}

void Bmatch_NtkVerifyReportError(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch &MI, vMatch &MO) {
    Vec_Ptr_t *vNodes;
    Abc_Obj_t *pNode;
    int *pValues1, *pValues2;
    int nErrors, nPrinted, i, iNode = -1;
    int *pModel2 = ABC_ALLOC(int, Abc_NtkPiNum(pNtk2));

    // get the input pattern of pNtk2
    for (i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
        for (auto &p : MI[i]) {
            pModel2[p.var()] = (p.sign()) ? !pModel1[i] : pModel1[i];
        }
    }
    for (auto &p : MI[Abc_NtkPiNum(pNtk1)]) {
        pModel2[p.var()] = (p.sign()) ? 0 : 1;
    }

    // get the CO values under this model
    pValues1 = Abc_NtkVerifySimulatePattern(pNtk1, pModel1);
    pValues2 = Abc_NtkVerifySimulatePattern(pNtk2, pModel2);
    // count the mismatches
    nErrors = 0;
    for (i = 0; i < Abc_NtkCoNum(pNtk1); i++) {
        for (auto &p : MO[i])
            nErrors += (int)(pValues1[i] != (p.sign() ? !pValues2[p.var()] : pValues2[p.var()]));
    }
    printf("Verification failed for at least %d outputs: ", nErrors);
    // print the first 3 outputs
    nPrinted = 0;
    for (i = 0; i < Abc_NtkCoNum(pNtk1); i++) {
        for (auto &p : MO[i]) {
            if (pValues1[i] != (p.sign() ? !pValues2[p.var()] : pValues2[p.var()])) {
                if (iNode == -1)
                    iNode = i;
                printf(" %s", Abc_ObjName(Abc_NtkCo(pNtk1, i)));
                if (++nPrinted == 3)
                    break;
            }
        }
    }
    if (nPrinted != nErrors)
        printf( " ..." );
    printf( "\n" );
    // report mismatch for the first output
    if (iNode >= 0 && MO[iNode].size() == 1) {
        printf( "Output %s: Value in Network1 = %d. Value in Network2 = %d.\n", 
            Abc_ObjName(Abc_NtkCo(pNtk1, iNode)), pValues1[iNode], (MO[iNode][0].sign()) ? !pValues2[MO[iNode][0].var()] : pValues2[MO[iNode][0].var()]);
        printf( "Input pattern: " );
        // collect PIs in the cone
        pNode = Abc_NtkCo(pNtk1, iNode);
        vNodes = Abc_NtkNodeSupport(pNtk1, &pNode, 1);
        // set the PI numbers
        Abc_NtkForEachCi(pNtk1, pNode, i)
            pNode->pCopy = (Abc_Obj_t *)(ABC_PTRINT_T)i;
        // print the model
        if (Vec_PtrSize(vNodes)) {
            pNode = (Abc_Obj_t *)Vec_PtrEntry(vNodes, 0);
            if (Abc_ObjIsCi(pNode)) {
                Vec_PtrForEachEntry(Abc_Obj_t *, vNodes, pNode, i) {
                    assert(Abc_ObjIsCi(pNode));
                    printf(" %s=%d", Abc_ObjName(pNode), pModel1[(int)(ABC_PTRINT_T)pNode->pCopy]);
                }
            }
        }
        printf("\n");
        Vec_PtrFree(vNodes);
    }
    ABC_FREE(pValues1);
    ABC_FREE(pValues2);
    ABC_FREE(pModel2);
}

ABC_NAMESPACE_IMPL_END