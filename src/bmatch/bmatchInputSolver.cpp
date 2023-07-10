#include "bmatch.hpp"
#include "bmatchQbf.hpp"
// #include "print.hpp"

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

// bmatchQbf.cpp
extern Bmatch_Qbf_Man_t *Bmatch_Gia_QbfAlloc(Gia_Man_t *pGia, int nPars, int fVerbose);
extern void Bmatch_Gia_QbfFree( Bmatch_Qbf_Man_t * p );
extern int Bmatch_Gia_QbfSolveValue(Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);
extern int Bmatch_Gia_QbfSolveValueInt(Bmatch_Qbf_Man_t *p, Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);

void Bmatch_InitControllableInputMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
void Bmatch_InitControllableInputOutputMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_InitInputControl(Bmatch_Man_t *pMan, int offset);

int Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_ApplyInputSolverRowConstraint(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_PruneInputSolverByBusOrdered(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverBySymmetryProperty(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverByUnate(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverBySymmetry(Bmatch_Man_t *pMan, vMatch &MI);
int Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI, vMatch& MO);
InputMapping Bmatch_HeuristicSolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
InputMapping Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);
InputMapping Bmatch_SolveInputQbf(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

int Bmatch_PruneSynthesizerByImpossibleMI(CaDiCaL::Solver *pSatSyn, Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
int Bmatch_PruneSynthesizerByFunSupport(CaDiCaL::Solver *pSatSyn, Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

void Bmatch_EncodeControlSignal(int row, int col, int nControlPi, int fCompl, AutoBuffer<int> &pLits);

#ifdef __cplusplus
}
#endif

void Bmatch_InitControllableInputOutputMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    CaDiCaL::Solver *pMiterSolver = pMan->pMiterSolverNew;
    if (pMiterSolver) Bmatch_sat_solver_delete(pMiterSolver);
    pMiterSolver = Bmatch_ControllableInputOutputSat(pNtk1, pNtk2, pMan->controlOffset);

    pMan->pMiterSolverNew = pMiterSolver;
}

void Bmatch_InitControllableInputMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    int &controlOffset = pMan->controlOffset;
    CaDiCaL::Solver *pMiterSolver = pMan->pMiterSolver;
    if (pMiterSolver) Bmatch_sat_solver_delete(pMiterSolver);
    pMiterSolver = Bmatch_ControllableInputSat(pNtk1, pNtk2, MO, pMan->controlOffset);

    Bmatch_InitInputControl(pMan, controlOffset);

    pMan->pMiterSolver = pMiterSolver;
}

void Bmatch_InitInputControl(Bmatch_Man_t *pMan, int offset) {
    auto &inputControl = pMan->inputControl;
    inputControl.resize(pMan->mi * pMan->ni);

    for (int i = 0; i < pMan->ni; ++i) {
        for (int j = 0; j < pMan->mi; ++j) {
            inputControl[i * pMan->mi + j] = Bmatch_toLitCond(i * pMan->mi + j + offset, 1);
        }
    }
}

InputMapping Bmatch_SolveInputQbf(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    Abc_Ntk_t *pNtkMiter = Bmatch_NtkControllableInputMiter(pNtk1, pNtk2, MO, 1);
    printf("Solving miter nodes: %d\n", Abc_NtkNodeNum(pNtkMiter));
    Aig_Man_t *pAig = Abc_NtkToDar(pNtkMiter, 0, 0);
    Gia_Man_t *pGia = Gia_ManFromAig(pAig);

    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1))); // nPi + const (not include inv)
    int nPars = (nControlPi + 1) * Abc_NtkPiNum(pNtk2);

    Vec_Int_t *vControl = Vec_IntAlloc(nPars);
    for (int i = 0; i < nPars; ++i)
        Vec_IntSetEntry(vControl, i, i);
    // int result = Gia_QbfSolveValue(pGia, vControl, nPars, 1024, 0, 100, 0, 0);
    // int result = Bmatch_Gia_QbfSolveValue(pGia, vControl, nPars, 1024, 0, 100, 0, 0);
    Bmatch_Qbf_Man_t *pQbfMan = Bmatch_Gia_QbfAlloc(pGia, nPars, 0);
    CaDiCaL::Solver  *pSatSyn = pQbfMan->pSatSyn;

    Bmatch_PruneSynthesizerByImpossibleMI(pSatSyn, pMan, pNtk1, pNtk2, MO);
    // Bmatch_PruneSynthesizerByFunSupport(pSatSyn, pMan, pNtk1, pNtk2, MO);

    int result = Bmatch_Gia_QbfSolveValueInt(pQbfMan, pGia, vControl, nPars, 1024, 0, 100, 0, 0);
    Bmatch_Gia_QbfFree(pQbfMan);

    Gia_ManStop(pGia);
    Aig_ManStop(pAig);
    Abc_NtkDelete(pNtkMiter);

    int RetValue = 1;
    if (result == 1 || result == -1) {
        RetValue = 0;
    } else if (result == 0) {
        for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
            int decode = 0;
            for (int j = nControlPi - 1; j >= 0; --j) {
                int value = (1 << j) * Vec_IntEntry(vControl, i * (nControlPi + 1) + j);
                decode += value;
                if (decode >= Abc_NtkPiNum(pNtk2) + 1)
                    decode -= value;
                // printf("decode + value: %d Abc_NtkPiNum(pNtk1): %d\n", decode + value, Abc_NtkPiNum(pNtk1));
            }
            decode = decode * 2 + Vec_IntEntry(vControl, i * (nControlPi + 1) + nControlPi);
            // printf("Decode: %d %d\n", decode, pMan->mi);
            for (int j = 0; j < pMan->mi; ++j) {
                if (j == decode) {
                    // printf("Resize: %d\n", resize);
                    MI[decode / 2].push_back(Literal(i, decode & 1));
                }
            }
        }
    }

    Vec_IntFree(vControl);

    return {RetValue, MI};
}

InputMapping Bmatch_HeuristicSolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());

}

InputMapping Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    auto *pSolver = pMan->pInputSolver;
    // pSolver->set("verbose", 3);

    int status = Bmatch_sat_solver_simplify(pSolver, 0);
    if (status == 20) return {0, MI};
    status = Bmatch_sat_solver_solve(pSolver, bLits, eLits, 0, 0, 0, 0);
    if (status == 20 || status == 0) return {0, MI};

    int n = pMan->ni;
    int m = pMan->mi;

    if (fVerbose) { printf("       "); for (int j = 0; j < m - 2; ++j) { printf("%c%-3s", ((j & 1) != 0) ? '~' : ' ', Abc_ObjName(Abc_NtkPi(pNtk1, j / 2))); } printf(" 1   0\n"); }
    for (int i = 0; i < n; ++i) {
        if (fVerbose) printf("%3s: ", Abc_ObjName(Abc_NtkPi(pNtk2, i)));
        for (int j = 0; j < m; ++j) {
            if (fVerbose) printf(" %3d", Bmatch_sat_solver_var_value(pSolver, i * m + j));
            if (Bmatch_sat_solver_var_value(pSolver, i * m + j)) {
                MI[j / 2].emplace_back(i, (int)((j & 1) != 0));
            }
        }
        if (fVerbose) printf("\n");
    }
    if (fVerbose) printf("\n");

    return {1, MI};
}

void Bmatch_EncodeControlSignal(int row, int col, int nControlPi, int fCompl, AutoBuffer<int> &pLits) {
    int code = col / 2;
    int sign = col & 1;
    for (int k = 0; k < nControlPi; ++k, code >>= 1) {
        pLits[k] = Bmatch_toLitCond(row * (nControlPi + 1) + k, (fCompl) ? code & 1 : !(code & 1));
    }
    pLits[nControlPi] = Bmatch_toLitCond(row * (nControlPi + 1) + nControlPi, (fCompl) ? sign : !sign);
}

int Bmatch_PruneSynthesizerByImpossibleMI(CaDiCaL::Solver *pSatSyn, Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    auto &impossibleMI = pMan->impossibleMI;
    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));
    AutoBuffer<int> pLits(nControlPi + 1);

    for (int i = 0; i < pMan->ni; ++i) {
        for (int j = 0; j < pMan->mi; ++j) {
            if (impossibleMI[i * pMan->mi + j] == 0) continue;
            Bmatch_EncodeControlSignal(i, j, nControlPi, 1, pLits);
            Bmatch_sat_solver_addclause(pSatSyn, pLits, pLits + pLits.size());
        }
    }
    // for (int i = 0; i < pMan->ni; ++i) {
    //     for (int j = pMan->mi; j < (1 << (nControlPi + 1)); ++j) {
    //         Bmatch_EncodeControlSignal(i, j, nControlPi, 1, pLits);
    //         Bmatch_sat_solver_addclause(pSatSyn, pLits, pLits + pLits.size());
    //     }
    // }

    return 1;
}

void Bmatch_sat_solver_add_dnf_clause(CaDiCaL::Solver *pSolver, std::vector<AutoBuffer<int> > &dnf) {
    AutoBuffer<int> ORs(dnf.size());
    for (int i = 0; i < dnf.size(); ++i) {
        auto &d = dnf[i];
        assert(d.size() > 0);
        // Tseytin transformation
        int GateC = d[0];
        if (d.size() > 1) {
            AutoBuffer<int> pLits(d.size() + 1);
            GateC = Bmatch_toLit(Bmatch_sat_solver_nvars(pSolver));
            for (int i = 0; i < d.size(); ++i) {
                pLits[i] = -d[i];
            }
            pLits[d.size()] = GateC;
            Bmatch_sat_solver_addclause(pSolver, pLits, pLits + pLits.size());
            pLits[0] = -GateC;
            for (int i = 0; i < d.size(); ++i) {
                pLits[1] = d[i];
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
            }
        }
        // for (int i = 1; i < d.size(); ++i) {
        //     int GateB = d[i];
        //     int GateC = Bmatch_toLit(Bmatch_sat_solver_nvars(pSolver));
        //     pLits[0] = -GateC;
        //     pLits[1] = GateA;
        //     Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
        //     pLits[1] = GateB;
        //     Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
        //     pLits[0] = GateC;
        //     pLits[1] = -GateA;
        //     pLits[2] = -GateB;
        //     Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 3);
        //     GateA = GateC;
        // }
        ORs[i] = GateC;
    }
    Bmatch_sat_solver_addclause(pSolver, ORs, ORs + ORs.size());
}

int Bmatch_PruneSynthesizerByFunSupport(CaDiCaL::Solver *pSatSyn, Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    int ni = pMan->ni;
    int mi = pMan->mi;
    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));

    auto &oSupp1 = pMan->oFuncSupp1;
    auto &oSupp2 = pMan->oFuncSupp2;

    std::vector<AutoBuffer<int> > dnf;

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];       // functional support of fi
            auto &cond2 = oSupp2[gi.var()]; // functional support of gi

            for (auto &m : cond1) {         // every support(fi)
                int i = 0;
                for (auto &n : cond2) {     // every support(gi)
                    AutoBuffer<int> pLitsA(nControlPi + 1);
                    Bmatch_EncodeControlSignal(n, m * 2, nControlPi, 0, pLitsA);
                    AutoBuffer<int> pLitsB(nControlPi + 1);
                    Bmatch_EncodeControlSignal(n, m * 2 + 1, nControlPi, 0, pLitsB);
                    dnf.emplace_back(std::move(pLitsA));
                    dnf.emplace_back(std::move(pLitsB));
                }
                Bmatch_sat_solver_add_dnf_clause(pSatSyn, dnf);
            }
        }
    }

    return 1;
}

int Bmatch_PruneInputSolverByBusOrdered(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    if (pMan->BI1.empty() || pMan->BI2.empty()) return 1;

    auto *pSolver = pMan->pInputSolver;

    auto &BI1 = pMan->BI1;
    auto &BI2 = pMan->BI2;

    int mi = pMan->mi;

    for (auto &bi1 : BI1) {
        for (auto &bi2 : BI2) {
            for (int i = 0; i < bi1.size(); ++i) {
                int xi = bi1[i];
                for (int j = 0; j < bi2.size(); ++j) {
                    int yi = bi2[j];

                    AutoBuffer<int, 3> pLits;
                    for (int invA = 0; invA < 2; ++invA) {
                        int match1 = yi * mi + xi * 2 + invA;
                        pLits[0] = Bmatch_toLitCond(match1, 1);

                        for (int k = 0; k < bi1.size(); ++k) {
                            if (k == i) continue;
                            for (int l = 0; l < bi2.size(); ++l) {
                                if (l == j) continue;

                                for (int invB = 0; invB < 2; ++invB) {
                                    int match2 = l * mi + k * 2 + invB;
                                    pLits[1] = Bmatch_toLitCond(match2, 1);

                                    for (int m = 0; m < bi1.size(); ++m) {
                                        if (m == i || m == k) continue;
                                        for (int n = 0; n < bi2.size(); ++n) {
                                            if (n == j || n == l) continue;

                                            if ((l > j && (m > k || m < i) && (n < l && n > j)) || (l < j && (m < k || m > i) && (n > l))) {
                                                for (int invC = 0; invC < 2; ++invC) {
                                                    int disable = n * mi + m * 2 + invC;
                                                    pLits[2] = Bmatch_toLitCond(disable, 1);

                                                    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 3);
                                                }
                                            } 
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 1;
}

int Bmatch_PruneInputSolverByBusExactMap(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    if (pMan->BI1.empty() || pMan->BI2.empty()) return 1;

    auto *pSolver = pMan->pInputSolver;
    auto &BI1 = pMan->BI1;
    auto &BI2 = pMan->BI2;

    int mi = pMan->mi;

    for (auto &bi1 : BI1) {
        for (auto &bi2 : BI2) {
            if (bi1.size() != bi2.size()) continue;

            int size = bi1.size();
            AutoBuffer<int, 3> pLits;
            for (int i = 0; i < size; ++i) {
                int xi = bi1[i], yi = bi2[i];
               
                for (int k = 0; k < 2; ++k) {
                    pLits[0] = Bmatch_toLitCond(yi * mi + xi * 2 + k, 1);
                    for (int j = 0; j < size; ++j) {
                        if (i == j) continue;
                        int xj = bi1[j], yj = bi2[j];

                        //printf("(%d %d) -> (%d %d) (%d %d)\n", yi, xi * 2 + k, yj, xj * 2, yj, xj * 2 + 1);
                        pLits[1] = Bmatch_toLit(yj * mi + xj * 2);
                        pLits[2] = Bmatch_toLit(yj * mi + xj * 2 + 1);
                        Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 3);
                    }
                }
            }
            for (int i = 0; i < size; ++i) {
                int xi = bi1[i], yi = bi2[size - 1 - i];
               
                for (int k = 0; k < 2; ++k) {
                    pLits[0] = Bmatch_toLitCond(yi * mi + xi * 2 + k, 1);
                    for (int j = 0; j < size; ++j) {
                        if (i == size - 1 - j) continue;
                        int xj = bi1[j], yj = bi2[size - 1 - j];

                        //printf("(%d %d) -> (%d %d) (%d %d)\n", yi, xi * 2 + k, yj, xj * 2, yj, xj * 2 + 1);
                        pLits[1] = Bmatch_toLit(yj * mi + xj * 2);
                        pLits[2] = Bmatch_toLit(yj * mi + xj * 2 + 1);
                        Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 3);
                    }
                }
            }
        }
    }

    return 1;
}

int Bmatch_PruneInputSolverByUnate(Bmatch_Man_t *pMan, vMatch &MO) {
    auto *pSolver = pMan->pInputSolver;

    Mat &unateMat1 = pMan->unateMat1;
    Mat &unateMat2 = pMan->unateMat2;
    auto &impossibleMI = pMan->impossibleMI;

    int ret = 1;

    int mi = pMan->mi;
    int ni = pMan->ni;

    int Lit;

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (int xi = 0; xi < mi / 2 - 1; ++xi) {
                int unateness1 = unateMat1[fi][xi];

                for (int yi = 0; yi < ni; ++yi) {
                    int unateness2 = unateMat2[gi][yi];

                    if (unateness1 == -1 || unateness2 == -1 || unateness1 == 3 || unateness2 == 3) continue;

                    if (unateness1 == unateness2) {
                        Lit = Bmatch_toLitCond(yi * mi + xi * 2 + (1 - g.sign()), 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        impossibleMI[yi * mi + xi * 2 + (1 - g.sign())] = 1;
                        // printf("(%d, %d) ", yi, xi * 2 + 1);
                    } else if (unateness1 != unateness2) {
                        Lit = Bmatch_toLitCond(yi * mi + xi * 2 + g.sign(), 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        impossibleMI[yi * mi + xi * 2 + g.sign()] = 1;
                        // printf("(%d, %d) ", yi, xi * 2);
                    }
                }
            }
            // printf("\n");
        }
    }

    return ret;
}

int Bmatch_PruneInputSolverBySymmetryProperty(Bmatch_Man_t *pMan, vMatch &MO) {
    if (pMan->vSymm1.empty() || pMan->vSymm2.empty()) return 1;
    auto *pSolver = pMan->pInputSolver;

    vSymm &symm1 = pMan->vSymm1;
    vSymm &symm2 = pMan->vSymm2;

    int ret = 1;
    int mi = pMan->mi;

    AutoBuffer<int> pLits(pMan->mi * pMan->ni);

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (auto &vSymm1 : symm1[fi]) {      // vector of fi input symmetry group
                for (auto &vSymm2 : symm2[gi]) {  // vector of gi input symmetry group
                    for (auto &xi : vSymm1) {
                        for (auto &yi : vSymm2) {
                            for (int j = 0; j < 2; ++j) {
                                int start = 0;
                                pLits[start++] = Bmatch_toLitCond(yi * mi + xi * 2 + j, (1 - j));
                                pLits[start++] = Bmatch_toLitCond(yi * mi + xi * 2 + (1 - j), j);
                                // printf("%c(%d, %d) ", (1 - j) ? '!' : ' ', yi, xi * 2 + j);
                                // printf("%c(%d, %d) ", (j) ? '!' : ' ', yi, xi * 2 + (1 - j));
                                for (auto &xir : vSymm1) {
                                    if (xi == xir) continue;
                                    for (auto &yia : vSymm2) {
                                        if (yi == yia) continue;
                                        pLits[start++] = Bmatch_toLit(yia * mi + xir * 2);
                                        pLits[start++] = Bmatch_toLit(yia * mi + xir * 2 + 1);
                                        // printf("(%d, %d) ", yia, xir * 2);
                                        // printf("(%d, %d) ", yia, xir * 2 + 1);
                                    }
                                }
                                // printf("\n");
                                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + start);
                            }
                        }
                    }
                }
            }
            break;
        }
    }

    return ret;
}

int Bmatch_PruneInputSolverBySymmetry(Bmatch_Man_t *pMan, vMatch &MI) {
    if (pMan->vSymmPair1.empty() || pMan->vSymmPair2.empty()) return 1;
    if (MI.empty()) return 1;
    int ret = 1;

    auto *pSolver = pMan->pInputSolver;

    vSymmPair &symm2 = pMan->vSymmPair2;

    int mi = pMan->mi;

    int pLits[2];
    for (int i = 0; i < MI.size(); ++i) {
        for (auto &p : MI[i]) {
            for (auto &s : symm2) { // loop symmetry pairs
                if ((p.var() == s.first || p.var() == s.second)) {
                    int symm_port = (p.var() == s.first) ? s.second : s.first;
                    int find = 0;
                    for (int j = 0; j < MI.size(); ++j) {
                        for (auto &q : MI[j]) {
                            if (q.var() == symm_port) { // find the corresponding pair
                                find = 1;
                                pLits[0] = Bmatch_toLitCond(p.var() * mi + j * 2 + q.sign(), 1);
                                pLits[1] = Bmatch_toLitCond(q.var() * mi + i * 2 + p.sign(), 1);
                                // printf("(%d, %d) ", p.var(), j * 2 + q.sign());
                                // printf("(%d, %d)\n", q.var(), i * 2 + p.sign());
                                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
                                break;
                            }
                        }
                        if (find) break;
                    }
                }
            }
        }
    }

    return ret;
}

int Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI, vMatch& MO) {
    if (!pModel1) return 1;

    int ret = 1;
    auto *pSolver = pMan->pInputSolver;
    AutoBuffer<int> pModel2(Abc_NtkPiNum(pNtk2));
    // get the input pattern of pNtk2
    for (int i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
        for (auto &p : MI[i]) {
            pModel2[p.var()] = (p.sign()) ? !pModel1[i] : pModel1[i];
        }
    }
    for (auto &p : MI[Abc_NtkPiNum(pNtk1)]) {
        pModel2[p.var()] = (p.sign()) ? 0 : 1;
    }

    int mi = pMan->mi;
    int k = 0;
    AutoBuffer<int> pLits((Abc_NtkPiNum(pNtk1) + 2) * Abc_NtkPiNum(pNtk2));

    auto &sRedund1 = pMan->sRedund1;
    auto &sRedund2 = pMan->sRedund2;
    auto &impossibleMI = pMan->impossibleMI;

    // counter example
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        if (sRedund2[i]) continue;
        for (int j = 0; j < Abc_NtkPiNum(pNtk1); ++j) {
            if (sRedund1[j]) continue;
            else if (impossibleMI[i * mi + j * 2] || impossibleMI[i * mi + j * 2 + 1]) continue;
            pLits[k++] = Bmatch_toLit(i * mi + j * 2 + (pModel2[i] == pModel1[j]));
            // printf("(%d %d) ", i, j * 2 + (pModel2[i] == pModel1[j]));
        }
        pLits[k++] = Bmatch_toLit((i * mi) + mi - 1 - (pModel2[i] == 0)); // don't know how this work???
    }
    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + k);

    // discard current
    // k = 0;
    // for (int i = 0; i < Abc_NtkPiNum(pNtk1) + 1; ++i) {
    //     for (auto &p : MI[i]) {
    //         pLits[k++] = Bmatch_toLitCond(p.var() * mi + i * 2 + p.sign(), 1);
    //     }
    // }
    // Bmatch_sat_solver_addclause(pSolver, pLits, pLits + k);

    ABC_FREE(pModel1);

    return ret;
}

int Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oFuncSupp1;
    auto &oSupp2 = pMan->oFuncSupp2;
    auto &impossibleMI = pMan->impossibleMI;

    AutoBuffer<int> pLits(ni * 2);

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];       // functional support of fi
            auto &cond2 = oSupp2[gi.var()]; // functional support of gi

            for (auto &m : cond1) {         // every support(fi)
                int i = 0;
                for (auto &n : cond2) {     // every support(gi)
                    pLits[i++] = Bmatch_toLit(n * mi + m * 2);
                    pLits[i++] = Bmatch_toLit(n * mi + m * 2 + 1);
                    if (cond1.size() == cond2.size()) { // if functional support is the same, then the input of gi cannot be const
                        int Lit = Bmatch_toLitCond(n * mi + mi - 1, 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        impossibleMI[n * mi + mi - 1] = 1;
                        Lit = Bmatch_toLitCond(n * mi + mi - 2, 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        impossibleMI[n * mi + mi - 2] = 1;
                    }
                }
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + i);
            }
        }
    }

    return ret;
}

int Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;
    auto &impossibleMI = pMan->impossibleMI;

    AutoBuffer<int> valid(ni * mi, 0);
    
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    valid[n * mi + m * 2] = 1;
                    valid[n * mi + m * 2 + 1] = 1;
                    // printf("(%d %d) (%d %d) ", n, m * 2, n, m * 2 + 1);
                }
                // printf("\n");
            }
        }
    }

    for (int n = 0; n < ni; ++n) {
        for (int m = 0; m < mi - 2; ++m) {
            if (valid[n * mi + m] == 0) {
                int Lit = Bmatch_toLitCond(n * mi + m, 1);
                Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                impossibleMI[n * mi + m] = 1;
            }
        }
    }

    return ret;
}

int Bmatch_ApplyInputSolverRowConstraint(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;

    int start = 0;
    int n = pMan->ni;
    int m = pMan->mi;
    auto &sRedund1 = pMan->sRedund1;
    auto &sRedund2 = pMan->sRedund2;
    auto &impossibleMI = pMan->impossibleMI;

    AutoBuffer<int> pLits(m);

    for (int i = 0; i < n; ++i) {
        // if input of cir2 is redundant
        if (sRedund2[i]) {
            pLits[0] = Bmatch_toLit(i * m + m - 1); // set it to constant
            Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
            continue;
        }
        // V(aij v bij)
        start = 0;
        for (int j = 0; j < m; ++j) {
            if (j < m - 2 && sRedund1[j / 2]) continue; // do not consider cir1's redundant input
            else if (impossibleMI[i * m + j]) continue; // impossible MI
            pLits[start++] = Bmatch_toLit(i * m + j);
        }
        Bmatch_sat_solver_addclause(pSolver, pLits, pLits + start);

        // VC~(row, 2)
        for (int j = 0; j < m - 1; ++j) {
            if (j < m - 2 && sRedund1[j / 2]) continue; // do not consider cir1's redundant input
            else if (impossibleMI[i * m + j]) continue; // impossible MI
            pLits[0] = Bmatch_toLitCond(i * m + j, 1);
            for (int k = j + 1; k < m; ++k) {
                pLits[1] = Bmatch_toLitCond(i * m + k, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
            }
        }
    }

    return ret;
}

int Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;
    if (pSolver) Bmatch_sat_solver_delete(pSolver);
    pSolver = Bmatch_sat_solver_new();
    // pSolver->fVerbose = 1;
    // pSolver->fPrintClause = 1;
    // pSolver->verbosity = 2;

    int n = pMan->ni = Abc_NtkPiNum(pNtk2);
    int m = pMan->mi = 2 * (Abc_NtkPiNum(pNtk1) + 1);

    Bmatch_sat_solver_setnvars(pSolver, n * m);

    // pMan->heuristicStage = 0;
    pMan->impossibleMI.resize(n * m);
    pMan->impossibleMI.fill(0);
    pMan->pInputSolver = pSolver;

    return ret;
}


ABC_NAMESPACE_IMPL_END
