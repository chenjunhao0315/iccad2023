#include <algorithm>

#include "bmatch.hpp"
#include "opt/sim/sim.h"

ABC_NAMESPACE_IMPL_START

#define AIG_OBJ_NAME2INDEX(TYPE, pObj)

#ifdef __cplusplus
extern "C" {
#endif

#define LINEAR_SEARCH_FUNC(TYPE) int Bmatch_LinearSearch##TYPE##Name2Index(Abc_Ntk_t *pNtk, Abc_Obj_t *pObj);
LINEAR_SEARCH_FUNC(Pi)
LINEAR_SEARCH_FUNC(Po)
#undef LINEAR_SEARCH_FUNC

void Bmatch_BusNameMaping(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_CalSuppAndSymm(Abc_Ntk_t *pNtk, vSupp &iFuncSupp, vSupp &oFuncSupp, vSymm &vSymm, vSymmPair &vSymmPair, int calSymm);
void Bmatch_CalStructSupp(vSupp &oStrSupp, Abc_Ntk_t *pNtk);
void Bmatch_CalRedundSupp(vSupp &rSupp, vSupp &oStrSupp, vSupp &oFuncSupp);
void Bmatch_CalSuppInfo(vSuppInfo& vSuppInfo, vSupp &oFuncSupp, vSupp &oStrSupp);
void Bmatch_CalCirRedund(Abc_Ntk_t *pNtk, vSupp &oStrSupp, AutoBuffer<int> &sRedund);
void Bmatch_CalCir2RedundWithGivenMapping(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, std::set<int> &sRedund);
void Bmatch_CalUnate(Abc_Ntk_t *pNtk, Mat &unateMat);
void Bmatch_CalEqual(vEqual &oEqual1, Abc_Ntk_t *pNtk);
void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

#ifdef __cplusplus
}
#endif

void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    // Setup network (Resythn or something)
    Bmatch_NtkSetup(pNtk1, pNtk2, option);

    // Convert bus information from string to index
    Bmatch_BusNameMaping(pMan, pNtk1, pNtk2);

    // Functional support and Symmetry information
    int calSymm = Abc_NtkPiNum(pNtk1) <= 20 && Abc_NtkPiNum(pNtk2) <= 20;
    printf("Calculate Functional support");
    if (calSymm) printf(" and Symm");
    printf(" information...");
    Bmatch_CalSuppAndSymm(pNtk1, pMan->iFuncSupp1, pMan->oFuncSupp1, pMan->vSymm1, pMan->vSymmPair1, calSymm);
    Bmatch_CalSuppAndSymm(pNtk2, pMan->iFuncSupp2, pMan->oFuncSupp2, pMan->vSymm2, pMan->vSymmPair2, calSymm);
    printf("Done!\n");

    // Structural support information
    printf("Calculate Structural support information...");
    Bmatch_CalStructSupp(pMan->oStrSupp1, pNtk1);
    Bmatch_CalStructSupp(pMan->oStrSupp2, pNtk2);
    printf("Done!\n");

    // Redundant support information
    printf("Calculate Redundant Information...");
    Bmatch_CalRedundSupp(pMan->oRedundSupp1, pMan->oStrSupp1, pMan->oFuncSupp1);
    Bmatch_CalRedundSupp(pMan->oRedundSupp2, pMan->oStrSupp2, pMan->oFuncSupp2);
    Bmatch_CalCirRedund(pNtk1, pMan->oStrSupp1, pMan->sRedund1);
    Bmatch_CalCirRedund(pNtk2, pMan->oStrSupp2, pMan->sRedund2);

    Bmatch_CalSuppInfo(pMan->vSuppInfo1, pMan->oFuncSupp1, pMan->oStrSupp1);
    Bmatch_CalSuppInfo(pMan->vSuppInfo2, pMan->oFuncSupp2, pMan->oStrSupp2);
    printf("Done!\n");

    // calculate unateness
    printf("Calculate Unateness Informatin...");
    Bmatch_CalUnate(pNtk1, pMan->unateMat1);
    Bmatch_CalUnate(pNtk2, pMan->unateMat2);
    printf("Done!\n");

    //equality
    printf("Calculate Equality Information...");
    Bmatch_CalEqual(pMan->oEqual1, pNtk1);
    Bmatch_CalEqual(pMan->oEqual2, pNtk2);
    printf("Done!\n");

    // Sensitivity or others
}

void Bmatch_CalUnate(Abc_Ntk_t *pNtk, Mat &unateMat) {
    Bmatch_RandomSimUnate(pNtk, unateMat, 10);
    for (int i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
        Bmatch_SatUnate(pNtk, unateMat, i, 1);
    }
}

void Bmatch_CalCirRedund(Abc_Ntk_t *pNtk, vSupp &oStrSupp, AutoBuffer<int> &sRedund) {
    sRedund.resize(Abc_NtkPiNum(pNtk));
    sRedund.fill(1);

    for (auto &s : oStrSupp)
        for (auto &p : s)
            sRedund[p] = 0;
}

void Bmatch_BusNameMaping(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int i;
    Abc_Obj_t *pObj;

    #define LINEAR_MAP(IO, sBIO, pNtk, BIO)                          \
    do {                                                             \
        for (auto &b : sBIO) {                                       \
            std::vector<int> bus;                                    \
            for (auto &p : b) {                                      \
                Abc_NtkForEach##IO(pNtk, pObj, i) {                  \
                    if (strcmp(p.c_str(), Abc_ObjName(pObj)) == 0) { \
                        bus.push_back(i);                            \
                    }                                                \
                }                                                    \
            }                                                        \
            if (!bus.empty()) BIO.emplace_back(std::move(bus));      \
        }                                                            \
    } while (0)

    LINEAR_MAP(Pi, pMan->sBIO1, pNtk1, pMan->BI1);
    LINEAR_MAP(Po, pMan->sBIO1, pNtk1, pMan->BO1);
    LINEAR_MAP(Pi, pMan->sBIO2, pNtk2, pMan->BI2);
    LINEAR_MAP(Po, pMan->sBIO2, pNtk2, pMan->BO2);

    #undef LINEAR_MAP

    pMan->sBIO1.clear();
    pMan->sBIO2.clear();
}

void Bmatch_CalSuppAndSymm(Abc_Ntk_t *pNtk, vSupp &iFuncSupp, vSupp &oFuncSupp, vSymm &vSymm, vSymmPair &vSymmPair, int calSymm) {
    Sym_Man_t *pSymMan;

    srand(0xABC);

    pSymMan = Sym_ManStart(pNtk, 0);

    if (calSymm) {
        pSymMan->nPairsTotal = pSymMan->nPairsRem = Sim_UtilCountAllPairs(pSymMan->vSuppFun, pSymMan->nSimWords, pSymMan->vPairsTotal);

        // detect symmetries using circuit structure
        Sim_SymmsStructCompute(pNtk, pSymMan->vMatrSymms, pSymMan->vSuppFun);
        Sim_UtilCountPairsAll(pSymMan);
        pSymMan->nPairsSymmStr = pSymMan->nPairsSymm;

        // detect symmetries using simulation
        for (int i = 1; i <= 1000; i++) {
            // simulate this pattern
            Sim_UtilSetRandom(pSymMan->uPatRand, pSymMan->nSimWords);
            Sim_SymmsSimulate(pSymMan, pSymMan->uPatRand, pSymMan->vMatrNonSymms);
            if (i % 50 != 0)
                continue;
            // check disjointness
            assert(Sim_UtilMatrsAreDisjoint(pSymMan));
            // count the number of pairs
            // Sim_UtilCountPairsAll(pSymMan);
        }

        // detect symmetries using SAT
        for (int i = 1; Sim_SymmsGetPatternUsingSat(pSymMan, pSymMan->uPatRand); i++) {
            // simulate this pattern in four polarities
            Sim_SymmsSimulate(pSymMan, pSymMan->uPatRand, pSymMan->vMatrNonSymms);
            Sim_XorBit(pSymMan->uPatRand, pSymMan->iVar1);
            Sim_SymmsSimulate(pSymMan, pSymMan->uPatRand, pSymMan->vMatrNonSymms);
            Sim_XorBit(pSymMan->uPatRand, pSymMan->iVar2);
            Sim_SymmsSimulate(pSymMan, pSymMan->uPatRand, pSymMan->vMatrNonSymms);
            Sim_XorBit(pSymMan->uPatRand, pSymMan->iVar1);
            Sim_SymmsSimulate(pSymMan, pSymMan->uPatRand, pSymMan->vMatrNonSymms);
            Sim_XorBit(pSymMan->uPatRand, pSymMan->iVar2);

            if ( i % 10 != 0 )
                continue;
            // check disjointness
            assert(Sim_UtilMatrsAreDisjoint(pSymMan));
            // count the number of pairs
            // Sim_UtilCountPairsAll(pSymMan);
        }

        // count the number of pairs
        // Sim_UtilCountPairsAll(pSymMan);

        auto insert_symm = [](std::vector<std::set<int> > &symm_groups, int i, int j) -> void {
            for (auto &g : symm_groups) {
                if (g.count(i)) { g.insert(j); return; }
                else if (g.count(j)) { g.insert(i); return; }
            }
            symm_groups.emplace_back(std::set<int>{i, j});
        };

        for (int i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
            int j, k, Index1, Index2;
            Extra_BitMat_t *pMat = (Extra_BitMat_t *)Vec_PtrEntry(pSymMan->vMatrSymms, i);
            Vec_Int_t *vSupport = Vec_VecEntryInt(pSymMan->vSupports, i);
            std::vector<std::set<int> > symm_groups;

            // symmetry
            Vec_IntForEachEntry(vSupport, j, Index1) {
                Vec_IntForEachEntryStart(vSupport, k, Index2, Index1 + 1) {
                    if (Extra_BitMatrixLookup1(pMat, j, k)) {
                        insert_symm(symm_groups, j, k);
                    }
                }
            }
            vSymm.emplace_back(std::move(symm_groups));
        }
    }

    iFuncSupp.resize(Abc_NtkPiNum(pNtk));
    for (int i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
        int j, k, Index1, Index2;
        Vec_Int_t *vSupport = Vec_VecEntryInt(pSymMan->vSupports, i);
        std::set<int> ofuncSupp;

        // functional support
        for (j = 0; j < Abc_NtkPiNum(pNtk); ++j) {
            if (Sim_SuppFunHasVar(pSymMan->vSuppFun, i, j)) {
                ofuncSupp.insert(j);
                iFuncSupp[j].insert(i);
            }
        }
        oFuncSupp.emplace_back(std::move(ofuncSupp));
    }

    if (calSymm) {
        for (auto &symm : vSymm) {
            for (auto &s : symm) {
                int supp_check = 1;
                std::set<int> supp_prev;
                for (auto &p : s) {
                    auto &supp = iFuncSupp[p];
                    if (supp_prev.empty()) {
                        supp_prev = supp;
                    } else if (supp_prev == supp || supp.empty()) {
                        
                    } else {
                        supp_check = 0;
                        break;
                    }
                }
                if (supp_check == 1) {
                    std::vector<int> symm_set(s.begin(), s.end());
                    for (int i = 0; i < symm_set.size() - 1; ++i) {
                        for (int j = i + 1; j < symm_set.size(); ++j) {
                            vSymmPair.emplace_back(symm_set[i], symm_set[j]);
                        }
                    }
                }
            }
        }
    }

    Sym_ManStop(pSymMan);
}

void Bmatch_CalSuppInfo(vSuppInfo& vSuppInfo, vSupp &oFuncSupp, vSupp &oStrSupp) {
    assert(oFuncSupp.size() == oStrSupp.size());

    for (int i = 0; i < oFuncSupp.size(); ++i) {
        vSuppInfo.emplace_back(std::make_tuple(i, oFuncSupp[i].size(), oStrSupp[i].size()));
    }

    std::sort(vSuppInfo.begin(), vSuppInfo.end(), [](std::tuple<int, int, int> const &t1, std::tuple<int, int, int> const &t2) { return std::get<SUPPFUNC>(t1) < std::get<SUPPFUNC>(t2); });
}

void Bmatch_CalEqual(vEqual &oEqual, Abc_Ntk_t *pNtk){
    int i, j, index1, index2;
    
    Abc_Obj_t *pNode1, *pNode2;
    std::vector<int> temp;
    std::vector<int> cal;
    Abc_NtkForEachPo( pNtk, pNode1, i ){
        index1 = Bmatch_LinearSearchPoName2Index(pNtk, pNode1);
        if(std::find(cal.begin(), cal.end(), index1) != cal.end()) continue;
        
        temp.emplace_back(index1);
        bool add = false;
        Abc_NtkForEachPo(pNtk, pNode2, j){
            if(pNode1 != pNode2){
                if((Abc_ObjFanin0(pNode1) == Abc_ObjFanin0(pNode2)) & (pNode1->fCompl0 == pNode2->fCompl0)){
                    add = true;
                    index2 = Bmatch_LinearSearchPoName2Index(pNtk, pNode2);
                    temp.emplace_back(index2);
                    cal.emplace_back(index2);
            }
            }
        }

        if (add) oEqual.emplace_back(temp);
        if (add) cal.emplace_back(index1);
        temp.clear();
    }
}

void Bmatch_CalStructSupp(vSupp &oStrSupp, Abc_Ntk_t *pNtk) {
    Vec_Ptr_t *vSupp, *vNodes;
    Abc_Obj_t *pObj, *pObjTemp;
    int i, j, index;

    oStrSupp.resize(Abc_NtkPoNum(pNtk));
    Abc_NtkForEachPo(pNtk, pObj, i) {
        vSupp = Abc_NtkNodeSupport(pNtk, &pObj, 1);
        // vNodes = Abc_NtkDfsNodes(pNtk, &pObj, 1);
        Vec_PtrForEachEntry(Abc_Obj_t *, vSupp, pObjTemp, j) {
            if (Abc_ObjIsPi(pObjTemp)) {
                index = Bmatch_LinearSearchPiName2Index(pNtk, pObjTemp);
                if (index != -1) oStrSupp[i].insert(index);
            }
        }
        Vec_PtrFree(vSupp);
        // Vec_PtrFree(vNodes);
    }
    Abc_NtkCleanMarkA(pNtk);
}

void Bmatch_CalRedundSupp(vSupp &rSupp, vSupp &oStrSupp, vSupp &oFuncSupp) {
    // StrSupp - FuncSupp
    assert(oStrSupp.size() == oFuncSupp.size());
    rSupp.resize(oStrSupp.size());

    for (int i = 0; i < oStrSupp.size(); ++i) {
        for (auto &p : oStrSupp[i]) {
            if (oFuncSupp[i].count(p) == 0) rSupp[i].insert(p);
        }
    }
}

#define LINEAR_SEARCH_FUNC_IMPL(TYPE)                                         \
int Bmatch_LinearSearch##TYPE##Name2Index(Abc_Ntk_t *pNtk, Abc_Obj_t *pObj) { \
    int i;                                                                    \
    Abc_Obj_t *pNode;                                                         \
    Abc_NtkForEach##TYPE(pNtk, pNode, i) {                                    \
        if (pNode == pObj) return i;                                          \
    }                                                                         \
    return -1;                                                                \
}
LINEAR_SEARCH_FUNC_IMPL(Pi)
LINEAR_SEARCH_FUNC_IMPL(Po)
#undef LINEAR_SEARCH_IMPL

ABC_NAMESPACE_IMPL_END
