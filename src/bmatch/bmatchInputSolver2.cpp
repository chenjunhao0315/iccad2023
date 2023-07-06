#include "bmatch.hpp"
#include "bmatchQbf.hpp"
// #include "print.hpp"

#include "base/io/ioAbc.h"
#include "aig/gia/giaAig.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

/*=== base/abci/abcDar.c ==============================================*/
extern Aig_Man_t * Abc_NtkToDar( Abc_Ntk_t * pNtk, int fExors, int fRegisters );
extern Abc_Ntk_t * Abc_NtkDC2( Abc_Ntk_t * pNtk, int fBalance, int fUpdateLevel, int fFanout, int fPower, int fVerbose );

/*=== aig/gia/giaAig.c ================================================*/
extern Gia_Man_t * Gia_ManFromAig( Aig_Man_t * p );

/*=== aig/gia/giaQbf.c ================================================*/
extern int Gia_QbfSolveValue( Gia_Man_t * pGia, Vec_Int_t * vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose );

extern void Bmatch_EncodeControlSignal(int row, int col, int nControlPi, int fCompl, AutoBuffer<int> &pLits);

// bmatchQbf.cpp
extern Bmatch_Qbf_Man_t *Bmatch_Gia_QbfAlloc(Gia_Man_t *pGia, int nPars, int fVerbose);
extern void Bmatch_Gia_QbfFree( Bmatch_Qbf_Man_t * p );
extern int Bmatch_Gia_QbfSolveValue(Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);
extern int Bmatch_Gia_QbfSolveValueInt(Bmatch_Qbf_Man_t *p, Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);

int Bmatch_LegalMI(Bmatch_Man_t *pMan, int n, int m);
void Bmatch_InitQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_FillPossibleMIbyStrSupp(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
void Bmatch_ReducePossibleMIbyUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
void Bmatch_ReducePossibleMIbySymmetry(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
Abc_Ntk_t *Bmatch_NtkQbfMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
InputMapping Bmatch_SolveQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

InputMapping Bmatch_SolveQbfInputSolver2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
void Bmatch_CnfReduce(std::vector<AutoBuffer<int> > &cnf);

#ifdef __cplusplus
}
#endif

Abc_Obj_t *Bmatch_NtkCreateMultiplexer2(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal) {
    int j;
    std::vector<Abc_Obj_t *> prevSignal = pSignal;
    std::vector<Abc_Obj_t *> currSignal;
    
    for (int i = 0; i < pControl.size(); ++i) {
        Abc_Obj_t *pC = pControl[i];
        for (j = 0; j < prevSignal.size(); j += 2) {
            Abc_Obj_t *p0 = prevSignal[j];
            Abc_Obj_t *p1 = ((j + 1) == prevSignal.size()) ? prevSignal[j] : prevSignal[j + 1];
            currSignal.push_back(Abc_AigMux(pMan, pC, p1, p0));
        }
        prevSignal = currSignal;
        currSignal.clear();
    }

    return prevSignal[0];
}

Abc_Ntk_t *Bmatch_NtkQbfMiterReduced(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, std::vector<int> &forceYi2Xi, vMatch &MO) {
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    Abc_Obj_t *pObj, *pObjTemp, *pObjNew;
    auto &MapReduceMI = pMan->MapReduceMI;
    int nControlPi = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1)))); // normal matrix
    int i;
    
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // Control PI
    std::vector<std::vector<Abc_Obj_t *> > controlPis;
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPi;
        for (int j = 0; j < nControlPi; ++j) {
            pObj = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "controlPi_%d_%d", i, j);
            Abc_ObjAssignName(pObj, Buffer, NULL);
            controlPi.emplace_back(std::move(pObj));
        }

        controlPis.emplace_back(std::move(controlPi));
    }
    assert(controlPis.size() == Abc_NtkPiNum(pNtk2));

    // Ntk1 PI
    std::vector<Abc_Obj_t *> Ntk1_Pis;
    Abc_NtkForEachPi(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "_ntk1");
        Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), Buffer);
        Ntk1_Pis.emplace_back(pObjNew);
        Ntk1_Pis.emplace_back(Abc_ObjNot(pObjNew));
        pObj->pCopy = pObjNew;
    }
    Ntk1_Pis.emplace_back(Abc_AigConst1(pNtkMiter));
    Ntk1_Pis.emplace_back(Abc_ObjNot(Abc_AigConst1(pNtkMiter)));

    // PO
    pObjNew = Abc_NtkCreatePo(pNtkMiter);
    Abc_ObjAssignName(pObjNew, "miter", NULL);

    // Create control input of Ntk2
    // Init map of reduce MI
    MapReduceMI = Mat(Abc_NtkPiNum(pNtk2), std::vector<int>());
    Abc_NtkForEachPi(pNtk2, pObj, i) {
        std::vector<Abc_Obj_t *> &pControl = controlPis[i];

        if (forceYi2Xi[i] < 0) { // non forced
            std::vector<Abc_Obj_t *> Reduced_Ntk1_Pis;
            for (int j = 0; j < 2 * (Abc_NtkPiNum(pNtk1) + 1); ++j) {
                if (Bmatch_LegalMI(pMan, i, j)) {
                    MapReduceMI[i].push_back(j);
                    Reduced_Ntk1_Pis.push_back(Ntk1_Pis[j]);
                }
            }
            if (Reduced_Ntk1_Pis.empty())
                pObjNew = Abc_AigConst1(pNtkMiter);
            else
                pObjNew = Bmatch_NtkCreateMultiplexer2((Abc_Aig_t *)pNtkMiter->pManFunc, pControl, Reduced_Ntk1_Pis);
        } else { // force connection
            pObjNew = Ntk1_Pis[forceYi2Xi[i]];
        }
        pObj->pCopy = pObjNew;
    }

    // Construct over Ntk1 and Ntk2
    std::vector<int> Ntk1PoUsed;
    std::vector<int> Ntk2PoUsed;
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i) {
        if (!MO[i].empty()) {
            Ntk1PoUsed.push_back(i);
            for (auto &gi : MO[i]) {
                Ntk2PoUsed.push_back(gi.var());
            }
        }
    }
    Bmatch_NtkBuildWithCone(pNtk1, pNtkMiter, Ntk1PoUsed);
    Bmatch_NtkBuildWithCone(pNtk2, pNtkMiter, Ntk2PoUsed);

    std::vector<Abc_Obj_t *> Xnors;
    Abc_NtkForEachCo(pNtk1, pObj, i) {
        if (!MO[i].empty()) {
            pObjTemp = Abc_ObjChild0Copy(pObj);
            for (int j = 0; j < MO[i].size(); ++j) {
                pObj = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
                pObj = (MO[i][j].sign()) ? Abc_ObjNot(pObj) : pObj;
                pObjNew = Abc_ObjNot(Abc_AigXor((Abc_Aig_t *)pNtkMiter->pManFunc, pObjTemp, pObj));
                Xnors.push_back(pObjNew);
            }
        }
    }
    pObj = Bmatch_NtkCreateAnd((Abc_Aig_t *)pNtkMiter->pManFunc, Xnors);
    Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pObj);

    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    // DC2
    // if (Abc_NtkNodeNum(pNtkMiter) > 500) {
        Abc_Ntk_t *pNtkTemp;
        pNtkMiter = Abc_NtkDC2(pNtkTemp = pNtkMiter, 0, 0, 1, 0, 0);
        Abc_NtkDelete(pNtkTemp);
    // }

    if (!Abc_NtkCheck(pNtkMiter)) {
        Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
        Abc_NtkDelete(pNtkMiter);

        return NULL;
    }

    // printf("Miter Node: %d -> ", Abc_NtkNodeNum(pNtkMiter));

    Io_Write(pNtkMiter, "miter_reduced.v", IO_FILE_VERILOG);

    return pNtkMiter;
}

int Bmatch_SolveQbfInputInt(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, std::vector<int> &forceYi2Xi, Vec_Int_t *vControl, vMatch &MO) {
    int nControlPi = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1)))); // normal matrix
    int nPars = nControlPi * Abc_NtkPiNum(pNtk2);

    pMan->possibleMI.resize(pMan->ni * pMan->mi);
    pMan->possibleMI.fill(0);
    Bmatch_FillPossibleMIbyStrSupp(pMan, pNtk1, pNtk2, MO);
    Bmatch_ReducePossibleMIbyUnate(pMan, pNtk1, pNtk2, MO);
    Bmatch_ReducePossibleMIbySymmetry(pMan, pNtk1, pNtk2, MO);
    Abc_Ntk_t *pNtkMiter = Bmatch_NtkQbfMiterReduced(pMan, pNtk1, pNtk2, forceYi2Xi, MO);
    Aig_Man_t *pAig = Abc_NtkToDar(pNtkMiter, 0, 0);
    Gia_Man_t *pGia = Gia_ManFromAig(pAig);

    Bmatch_Qbf_Man_t *pQbfMan = Bmatch_Gia_QbfAlloc(pGia, nPars, 0);
    CaDiCaL::Solver  *pSatSyn = pQbfMan->pSatSyn;

    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        // Set floating port to const
        if (pMan->MapReduceMI[i].empty()) {
            for (int j = 0; j < nControlPi; ++j) {
                int Lit = Bmatch_toLitCond(i * nControlPi + j, 1);
                Bmatch_sat_solver_addclause(pSatSyn, &Lit, &Lit + 1);
            }
        } else { // Forbid unsed index
            std::vector<AutoBuffer<int> > cnfs;
            for (int j = pMan->MapReduceMI[i].size(); j < (1 << nControlPi); ++j) {
                AutoBuffer<int> cnf(nControlPi, 0);
                for (int code = j, k = 0; k < nControlPi; code >>= 1, ++k) {
                    cnf[k] = Bmatch_toLitCond(i * nControlPi + k, code & 1);
                }
                cnfs.emplace_back(cnf);
            }
            Bmatch_CnfReduce(cnfs);
            for (auto &pLits : cnfs) {
                Bmatch_sat_solver_addclause(pSatSyn, pLits, pLits + pLits.size());
            }
        }
    }

    int result = Bmatch_Gia_QbfSolveValueInt(pQbfMan, pGia, vControl, nPars, 1024, 0, 100, 0, 1);
    // int result = Gia_QbfSolveValue(pGia, vControl, nPars, 1024, 0, 100, 0, 1);
    Bmatch_Gia_QbfFree(pQbfMan);

    Gia_ManStop(pGia);
    Aig_ManStop(pAig);
    Abc_NtkDelete(pNtkMiter);

    auto &ostrFunc2 = pMan->oStrSupp2;
    std::set<int> allStrInput;
    for (auto &fi : MO) {
        for (auto &gi : fi) {
            allStrInput.insert(ostrFunc2[gi.var()].begin(), ostrFunc2[gi.var()].end());
        }
    }

    int RetValue = 1;
    if (result == 1 || result == -1) {
        RetValue = 0;
    } else if (result == 0) {
        for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
            if (allStrInput.count(i) == 0 || pMan->MapReduceMI[i].empty()) continue;
            int decode = 0;
            for (int j = nControlPi - 1; j >= 0; --j) {
                int value = (1 << j) * Vec_IntEntry(vControl, i * nControlPi + j);
                decode += value;
                if (decode >= pMan->MapReduceMI[i].size())
                    decode -= value;
            }
            decode = pMan->MapReduceMI[i][decode];
            for (int j = 0; j < pMan->mi; ++j) {
                if (j == decode) {
                    forceYi2Xi[i] = decode;
                }
            }
        }
    }

    return RetValue;
}

// Iterative solving
InputMapping Bmatch_SolveQbfInputSolver2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    std::vector<int> forceYi2Xi(Abc_NtkPiNum(pNtk2), -1);

    // Convert to gi base
    std::vector<Literal> gi2fi(Abc_NtkPoNum(pNtk2), -1);
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            gi2fi[gi.var()] = Literal(fi, gi.sign());
        }
    }

    int nControlPi = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1)))); // normal matrix
    int nPars = nControlPi * Abc_NtkPiNum(pNtk2);
    Vec_Int_t *vControl = Vec_IntAlloc(nPars);
        for (int i = 0; i < nPars; ++i)
            Vec_IntSetEntry(vControl, i, i);

    vMatch MO_success(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
    int RetValue = 0;
    for (int i = 0; i < gi2fi.size(); ++i) { // Iterative matching
        if (gi2fi[i] < 0) continue;

        // Stage 1 (Fix previous input and try to add new one output port)
        printf("Step: %d -> Stage1 -> ", i);
        vMatch MO_stage1(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
        MO_stage1[gi2fi[i].var()].push_back(Literal(i, gi2fi[i].sign()));

        int stage1 = Bmatch_SolveQbfInputInt(pMan, pNtk1, pNtk2, forceYi2Xi, vControl, MO_stage1);
        RetValue = stage1;
        printf("stage1: %s", stage1 ? "Equivalence\n" : "Non-Equivalence -> Relax and retry...\n");

        if (stage1 == 1) {
            MO_success[gi2fi[i].var()].push_back(Literal(i, gi2fi[i].sign()));
            continue;
        }

        // // Stage 2 (Relax all and try new one output port)
        // printf("Step: %d -> Stage2 -> ", i);
        // std::vector<int> freeYi2Xi(Abc_NtkPiNum(pNtk2), -1);
        // int stage2 = Bmatch_SolveQbfInputInt(pMan, pNtk1, pNtk2, freeYi2Xi, vControl, MO_stage1);
        // printf("stage2: %s", stage2 ? "Possible solution\n" : "Impossible solution\n");

        // if (stage2 == 0) break; // Impossible to find the match

        // Stage 3 (Relax and try to add add new one output port)
        vMatch MO_stage3 = MO_success;
        MO_stage3[gi2fi[i].var()].push_back(Literal(i, gi2fi[i].sign()));
        auto &relaxInput = pMan->oStrSupp2[i];
        std::set<int> affectedOutput;
        for (auto &yi : relaxInput) {
            for (auto &gi : pMan->iFuncSupp2[yi]) {
                for (auto &yip : pMan->oStrSupp2[gi]) {
                    forceYi2Xi[yip] = -1;
                }
            }
        }
        printf("Step: %d -> Stage3 -> ", i);
        int stage3 = Bmatch_SolveQbfInputInt(pMan, pNtk1, pNtk2, forceYi2Xi, vControl, MO_stage3);
        RetValue = stage3;
        printf("stage3: %s", stage3 ? "Equivalence\n" : "Non-Equivalence\n");

        if (stage3 == 0) break;

        MO_success = MO_stage3;
    }
    if (RetValue == 1) { // Find Match
        std::replace(forceYi2Xi.begin(), forceYi2Xi.end(), -1, 2 * Abc_NtkPiNum(pNtk1));
        for (int i = 0; i < forceYi2Xi.size(); ++i) {
            MI[(forceYi2Xi[i] % pMan->mi) / 2].push_back(Literal(i, (forceYi2Xi[i] % pMan->mi) & 1));
        }
    }

    Vec_IntFree(vControl);

    return {RetValue, MI};
}

inline int Bmatch_LegalMI(Bmatch_Man_t *pMan, int n, int m) {
    return pMan->possibleMI[n * pMan->mi + m];
}

void Bmatch_InitQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int nControlPi = pMan->nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1))); // nPi + const (not include inv)
    int mi = pMan->mi = (Abc_NtkPiNum(pNtk1) + 1) * 2;
    int ni = pMan->ni = Abc_NtkPiNum(pNtk2);
    auto &possibleMI = pMan->possibleMI;

    possibleMI.resize(ni * mi);
    possibleMI.fill(0);
}

Abc_Ntk_t *Bmatch_NtkQbfMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    Abc_Obj_t *pObj, *pObjTemp, *pObjNew;
    auto &MapReduceMI = pMan->MapReduceMI;
    int nControlPi = pMan->nControlPi;
    int i;
    
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // Control PI
    std::vector<std::vector<Abc_Obj_t *> > controlPis;
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPi;
        for (int j = 0; j < nControlPi; ++j) {
            pObj = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "controlPi_%d_%d", i, j);
            Abc_ObjAssignName(pObj, Buffer, NULL);
            controlPi.emplace_back(std::move(pObj));
        }

        controlPis.emplace_back(std::move(controlPi));
    }
    assert(controlPis.size() == Abc_NtkPiNum(pNtk2));

    // Ntk1 PI
    std::vector<Abc_Obj_t *> Ntk1_Pis;
    Abc_NtkForEachPi(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "_ntk1");
        Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), Buffer);
        Ntk1_Pis.emplace_back(pObjNew);
        Ntk1_Pis.emplace_back(Abc_ObjNot(pObjNew));
        pObj->pCopy = pObjNew;
    }
    Ntk1_Pis.emplace_back(Abc_AigConst1(pNtkMiter));
    Ntk1_Pis.emplace_back(Abc_ObjNot(Abc_AigConst1(pNtkMiter)));

    // PO
    pObjNew = Abc_NtkCreatePo(pNtkMiter);
    Abc_ObjAssignName(pObjNew, "miter", NULL);

    // Create control input of Ntk2
    // Init map of reduce MI
    MapReduceMI = Mat(Abc_NtkPiNum(pNtk2), std::vector<int>());
    Abc_NtkForEachPi(pNtk2, pObj, i) {
        std::vector<Abc_Obj_t *> &pControl = controlPis[i];
        std::vector<Abc_Obj_t *> Reduced_Ntk1_Pis;
        for (int j = 0; j < 2 * (Abc_NtkPiNum(pNtk1) + 1); ++j) {
            if (Bmatch_LegalMI(pMan, i, j)) {
                MapReduceMI[i].push_back(j);
                Reduced_Ntk1_Pis.push_back(Ntk1_Pis[j]);
            }
        }
        pObjNew = Bmatch_NtkCreateMultiplexer2((Abc_Aig_t *)pNtkMiter->pManFunc, pControl, Reduced_Ntk1_Pis);

        pObj->pCopy = pObjNew;
    }

    // Construct over Ntk1 and Ntk2
    std::vector<int> Ntk1PoUsed;
    std::vector<int> Ntk2PoUsed;
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i) {
        if (!MO[i].empty()) {
            Ntk1PoUsed.push_back(i);
            for (auto &gi : MO[i]) {
                Ntk2PoUsed.push_back(gi.var());
            }
        }
    }
    Bmatch_NtkBuildWithCone(pNtk1, pNtkMiter, Ntk1PoUsed);
    Bmatch_NtkBuildWithCone(pNtk2, pNtkMiter, Ntk2PoUsed);

    std::vector<Abc_Obj_t *> Xnors;
    Abc_NtkForEachCo(pNtk1, pObj, i) {
        if (!MO[i].empty()) {
            pObjTemp = Abc_ObjChild0Copy(pObj);
            for (int j = 0; j < MO[i].size(); ++j) {
                pObj = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
                pObj = (MO[i][j].sign()) ? Abc_ObjNot(pObj) : pObj;
                pObjNew = Abc_ObjNot(Abc_AigXor((Abc_Aig_t *)pNtkMiter->pManFunc, pObjTemp, pObj));
                Xnors.push_back(pObjNew);
            }
        }
    }
    pObj = Bmatch_NtkCreateAnd((Abc_Aig_t *)pNtkMiter->pManFunc, Xnors);
    Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pObj);

    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    // DC2
    if (Abc_NtkNodeNum(pNtkMiter) > 500) {
        Abc_Ntk_t *pNtkTemp;
        pNtkMiter = Abc_NtkDC2(pNtkTemp = pNtkMiter, 1, 0, 1, 0, 0);
        Abc_NtkDelete(pNtkTemp);
    }

    if (!Abc_NtkCheck(pNtkMiter)) {
        Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
        Abc_NtkDelete(pNtkMiter);

        return NULL;
    }

    // Io_Write(pNtkMiter, "miter_reduced.v", IO_FILE_VERILOG);

    return pNtkMiter;
}

void Bmatch_CnfReduce(std::vector<AutoBuffer<int> > &cnf) {
    for (auto iter = cnf.begin(); iter != cnf.end(); ) {
        auto &c1 = *iter;
        bool reduce = false;
        for (auto &c2 : cnf) {
            if (c1.size() != c2.size()) continue;

            int index = -1;
            for (int i = 0; i < c1.size(); ++i) {
                if (c1[i] != c2[i] && c1[i] == -c2[i]) {
                    if (true == reduce) {
                        reduce = false;
                        break;
                    } else {
                        index = i;
                        reduce = true;
                    }
                } else if (c1[i] != c2[i]) {
                    break;
                }
            }
            if (reduce) {
                for (int i = index; i < c2.size() - 1; ++i) {
                    c2[i] = c2[i + 1];
                }
                c2.resize(c2.size() - 1);
                break;
            }
        }
        if (reduce) {
            iter = cnf.erase(iter);
        } else {
            ++iter;
        }
    }
}

int Bmatch_PruneQbfSynthesizerByPossibleMI(CaDiCaL::Solver *pSatSyn, Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    auto &possibleMI = pMan->possibleMI;
    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));

    std::vector<AutoBuffer<int> > cnfs;
    for (int i = 0; i < pMan->ni; ++i) {
        for (int j = pMan->MapReduceMI[i].size() * 2; j < (1 << (nControlPi + 1)); ++j) {
        // for (int j = pMan->MapReduceMI[i].size() * 2; j < pMan->mi; ++j) {
            AutoBuffer<int> pLits(nControlPi + 1);
            Bmatch_EncodeControlSignal(i, j, nControlPi, 1, pLits);
            cnfs.emplace_back(pLits);
        }
    }

    Bmatch_CnfReduce(cnfs);

    for (auto &pLits : cnfs) {
        Bmatch_sat_solver_addclause(pSatSyn, pLits, pLits + pLits.size());
    }

    return 1;
}

InputMapping Bmatch_SolveQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    int nControlPi = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1)))); // normal matrix
    int nPars = nControlPi * Abc_NtkPiNum(pNtk2);
    Vec_Int_t *vControl = Vec_IntAlloc(0);
    std::vector<int> forceYi2Xi(Abc_NtkPiNum(pNtk2), -1);

    auto &ostrFunc2 = pMan->oStrSupp2;
    std::set<int> allStrInput;
    for (auto &fi : MO) {
        for (auto &gi : fi) {
            allStrInput.insert(ostrFunc2[gi.var()].begin(), ostrFunc2[gi.var()].end());
        }
    }
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        if (allStrInput.count(i) == 0)
            forceYi2Xi[i] = 2 * Abc_NtkPiNum(pNtk1);
    }

    int RetValue = Bmatch_SolveQbfInputInt(pMan, pNtk1, pNtk2, forceYi2Xi, vControl, MO);

    if (RetValue == 1) { // Find Match
        std::replace(forceYi2Xi.begin(), forceYi2Xi.end(), -1, 2 * Abc_NtkPiNum(pNtk1));
        for (int i = 0; i < forceYi2Xi.size(); ++i) {
            MI[(forceYi2Xi[i] % pMan->mi) / 2].push_back(Literal(i, (forceYi2Xi[i] % pMan->mi) & 1));
        }
    }

    Vec_IntFree(vControl);

    return {RetValue, MI};
}

void Bmatch_ReducePossibleMIbySymmetry(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    if (pMan->vSymm1.empty() || pMan->vSymm2.empty()) return;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &vSymm1 = pMan->vSymm1;
    auto &vSymm2 = pMan->vSymm2;
    auto &possibleMI = pMan->possibleMI;

    for (int fi = 0; fi < Abc_NtkPoNum(pNtk1); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (auto &symm1 : vSymm1[fi]) {
                for (auto &symm2 : vSymm2[gi]) {
                    if (symm1.size() > symm2.size()) {
                        for (auto &xi : symm1) {
                            for (auto &yi : symm2) {
                                possibleMI[yi * mi + 2 * xi + 0] = 0;
                                possibleMI[yi * mi + 2 * xi + 1] = 0;
                            }
                        }
                    }
                }
            }
        }
    }
}

void Bmatch_FillPossibleMIbyStrSupp(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;
    auto &possibleMI = pMan->possibleMI;
    
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    possibleMI[n * mi + m * 2] = 1;
                    possibleMI[n * mi + m * 2 + 1] = 1;
                }
                possibleMI[n * mi + mi - 1] = 1;
                possibleMI[n * mi + mi - 2] = 1;
            }
        }
    }
}

void Bmatch_ReducePossibleMIbyUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    auto &unateMat1 = pMan->unateMat1;
    auto &unateMat2 = pMan->unateMat2;
    auto &possibleMI = pMan->possibleMI;

    int mi = pMan->mi;
    int ni = pMan->ni;

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (int xi = 0; xi < mi / 2 - 1; ++xi) {
                int unateness1 = unateMat1[fi][xi];

                for (int yi = 0; yi < ni; ++yi) {
                    int unateness2 = unateMat2[gi][yi];

                    if (unateness1 == -1 || unateness2 == -1 || unateness1 == 3 || unateness2 == 3) continue;

                    if (unateness1 == unateness2) {
                        possibleMI[yi * mi + xi * 2 + (1 - g.sign())] = 0;
                    } else if (unateness1 != unateness2) {
                        possibleMI[yi * mi + xi * 2 + g.sign()] = 0;
                    }
                }
            }
        }
    }
}

ABC_NAMESPACE_IMPL_END