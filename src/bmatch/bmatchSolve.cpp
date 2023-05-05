#include "bmatch.hpp"

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

void Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO);
void Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO);
void Bmatch_PruneINputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI);
vMatch Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);

vGroup Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan);

#ifdef __cplusplus
}
#endif

// vMatch MI = {{Literal(1, false), Literal(3, false)}, {Literal(2, true)}, {Literal(0, true)}, {}};
// vMatch MO = {{Literal(0, false)}};

// vMatch MI = {{Literal(2, false)}, {Literal(3, false)}, {Literal(0, false)}, {Literal(1, false)}, {Literal(4, false)}, {}};
// vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {}, {}};

// vMatch MI = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}, {Literal(4, false)}, {}};
// vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}};

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    int maxIter = 50, iter = 0;
    EcResult result;

    // if (option & VERBOSE_MASK) Bmatch_PrintBusInfo(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintInputSense(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintOutputSupport(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintSymm(pMan, pNtk1, pNtk2);

    Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);

    // testing flow
    auto groups = Bmatch_SolveOutputGroup(pMan);

    if (option & VERBOSE_MASK) Bmatch_PrintOutputGroup(pNtk1, pNtk2, groups);

    vMatch MI;
    vMatch MO = {{Literal(0, false)}};
    Bmatch_PruneInputSolverByStrSupport(pMan, MO);
    Bmatch_PruneInputSolverByFuncSupport(pMan, MO);

    while (result.status != EQUIVALENT && iter++ < maxIter) {
        Bmatch_PruneINputSolverByCounterPart(pMan, pNtk1, pNtk2, result.model, MI);
        MI = Bmatch_SolveInput(pMan, pNtk1, pNtk2, NULL, NULL, 1);

        if (option & VERBOSE_MASK) Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO);

        assert(MI.size() == Abc_NtkPiNum(pNtk1) + 1);
        assert(MO.size() == Abc_NtkPoNum(pNtk1));
        result = Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO, 0);
        print("Status:", (result.status == EQUIVALENT) ? "EQUIVALENT" : "NON-EQUIVALENT");
    }

    if (result.status == EQUIVALENT) printf("Find matching at iteration %d!!!\n", iter);
    else printf("Reach maximum iteration!\n");
}

vMatch Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    sat_solver *pSolver = pMan->pInputSolver;

    int status = sat_solver_solve(pSolver, bLits, eLits, 0, 0, 0, 0);
    print(status == l_True ? "InputSolver: SAT" : "InputSolver: UNSAT");
    int *model = pSolver->model;
    if (status == l_False) return MI;

    int n = pMan->ni;
    int m = pMan->mi;

    if (fVerbose) { printf("       "); for (int j = 0; j < m - 2; ++j) { printf("%c%-3s", ((j & 1) != 0) ? '~' : ' ', Abc_ObjName(Abc_NtkPi(pNtk1, j / 2))); } printf(" 1   0\n"); }
    for (int i = 0; i < n; ++i) {
        if (fVerbose) printf("%3s: ", Abc_ObjName(Abc_NtkPi(pNtk2, i)));
        for (int j = 0; j < m; ++j) {
            if (fVerbose) printf(" %3d", model[i * m + j] == l_True);
            if (model[i * m + j] == l_True) {
                MI[j / 2].emplace_back(i, (int)((j & 1) != 0));
            }
        }
        if (fVerbose) printf("\n");
    }
    if (fVerbose) printf("\n");

    return MI;
}

void Bmatch_PruneINputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI) {
    if (!pModel1) return;

    sat_solver *pSolver = pMan->pInputSolver;
    int *pModel2 = ABC_ALLOC(int, Abc_NtkPiNum(pNtk2));
    // get the input pattern of pNtk2
    for (int i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
        for (auto &p : MI[i]) {
            pModel2[p.var()] = (p.sign()) ? !pModel1[i] : pModel1[i];
        }
    }

    int mi = pMan->mi;
    int k = 0;
    int *pLits = ABC_ALLOC(int, Abc_NtkPiNum(pNtk1) * Abc_NtkPiNum(pNtk2));

    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        for (int j = 0; j < Abc_NtkPiNum(pNtk1); ++j) {
            pLits[k++] = toLit(i * mi + j * 2 + (pModel2[i] == pModel1[j]));
        }
    }

    sat_solver_addclause(pSolver, pLits, pLits + k);

    ABC_FREE(pModel1);
    ABC_FREE(pModel2);
    ABC_FREE(pLits);
}

void Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    sat_solver *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oFuncSupp1;
    auto &oSupp2 = pMan->oFuncSupp2;

    int *pLits = ABC_ALLOC(int, ni * 2);

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &m : cond1) {
                int i = 0;
                for (auto &n : cond2) {
                    pLits[i++] = toLit(n * mi + m * 2);
                    pLits[i++] = toLit(n * mi + m * 2 + 1);
                }
                sat_solver_addclause(pSolver, pLits, pLits + i);
            }
        }
    }

    ABC_FREE(pLits);
}

void Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    sat_solver *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;

    int *valid = ABC_CALLOC(int, ni * mi);
    
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    valid[n * mi + m * 2] = 1;
                    valid[n * mi + m * 2 + 1] = 1;
                    
                }
                valid[n * mi + mi - 2] = 1;
                valid[n * mi + mi - 1] = 1;
            }
        }
    }

    for (int n = 0; n < ni; ++n) {
        for (int m = 0; m < mi; ++m) {
            if (valid[n * mi + m] == 0) {
                int Lit = toLitCond(n * mi + m, 1);
                sat_solver_addclause(pSolver, &Lit, &Lit + 1);
            }
        }
    }

    ABC_FREE(valid);
}

void Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    sat_solver *pSolver = pMan->pInputSolver;
    if (pSolver) sat_solver_delete(pSolver);
    pSolver = sat_solver_new();
    // pSolver->fVerbose = 1;
    // pSolver->fPrintClause = 1;
    // pSolver->verbosity = 2;

    int n = pMan->ni = Abc_NtkPiNum(pNtk2);
    int m = pMan->mi = 2 * (Abc_NtkPiNum(pNtk1) + 1);

    int *pLits = ABC_ALLOC(int, m);

    sat_solver_setnvars(pSolver, n * m);
    for (int i = 0; i < n; ++i) {
        // V(aij v bij)
        for (int j = 0; j < m; ++j) {
            pLits[j] = toLit(i * m + j);
        }
        sat_solver_addclause(pSolver, pLits, pLits + m);

        // VC~(row, 2)
        for (int j = 0; j < m - 1; ++j) {
            pLits[0] = toLitCond(i * m + j, 1);
            for (int k = j + 1; k < m; ++k) {
                pLits[1] = toLitCond(i * m + k, 1);
                sat_solver_addclause(pSolver, pLits, pLits + 2);
            }
        }
    }

    ABC_FREE(pLits);
    pMan->pInputSolver = pSolver;
}

vGroup Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan) {
    vGroup groups;
    auto &supp1 = pMan->vSuppInfo1;
    auto &supp2 = pMan->vSuppInfo2;
    int suppFunc1 = 0, suppFunc2 = 0;
    int n1 = supp1.size() - 1;
    int n2 = supp2.size() - 1;

    groups.emplace_back(std::vector<int>(), std::vector<int>());
    for (int i = n1, j = n2; i >= 0 || j >= 0; --i, --j) {
        auto &group = groups.back();
        if (i >= 0) group.first.emplace_back(std::get<PO>(supp1[i]));
        if (j >= 0) group.second.emplace_back(std::get<PO>(supp2[j]));

        suppFunc1 = std::max(suppFunc1, std::get<SUPPFUNC>(supp1[i]));
        suppFunc2 = (j - 1 >= 0) ? std::get<SUPPFUNC>(supp2[j - 1]) : std::get<SUPPFUNC>(supp2[0]);

        if (suppFunc1 > suppFunc2 && (i > 0 || j > 0)) {
            suppFunc1 = 0;
            groups.emplace_back(std::vector<int>(), std::vector<int>());
        }
    }

    return groups;
}

ABC_NAMESPACE_IMPL_END