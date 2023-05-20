#include "bmatch.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_InitControllableMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
void Bmatch_InitInputControl(Bmatch_Man_t *pMan, int offset);

int Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_ApplyInputSolverRowConstraint(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverBySymmetryProperty(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverByUnate(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverBySymmetry(Bmatch_Man_t *pMan, vMatch &MI);
int Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI, vMatch& MO);
InputMapping Bmatch_HeuristicSolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
InputMapping Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);

#ifdef __cplusplus
}
#endif

void Bmatch_InitControllableMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    int &controlOffset = pMan->controlOffset;
    CaDiCaL::Solver *pMiterSolver = pMan->pMiterSolver;
    if (pMiterSolver) Bmatch_sat_solver_delete(pMiterSolver);
    pMiterSolver = Bmatch_ControlSat(pNtk1, pNtk2, MO, pMan->controlOffset);

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
                        // printf("(%d, %d) ", yi, xi * 2 + 1);
                    } else if (unateness1 != unateness2) {
                        Lit = Bmatch_toLitCond(yi * mi + xi * 2 + g.sign(), 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
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
                }
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
