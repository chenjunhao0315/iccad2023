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
void Bmatch_CalSupp(Abc_Ntk_t *pNtk, vSupp &vSupp, vSense &iFuncSupp, vSense &oFuncSupp, vSense &oStrSupp, vSense &oRedundSupp, int option);
void Bmatch_CalStructSupp(vSense &oStrSupp, Abc_Ntk_t *pNtk);
void Bmatch_CalFuncSupp(vSupp &pVec, vSense &iSen, vSense &oFuncSupp, Abc_Ntk_t *pNtk, int option);
void Bmatch_CalRedundSupp(vSense &rSupp, vSense &oStrSupp, vSense &oFuncSupp);
void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

#ifdef __cplusplus
}
#endif

void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    // Setup network (Resythn or something)
    Bmatch_NtkSetup(pNtk1, pNtk2, option);

    // Convert bus information from string to index
    Bmatch_BusNameMaping(pMan, pNtk1, pNtk2);

    // Bmatch_CalSupp(pNtk1, pMan->suppFunc1, pMan->FI1, pMan->FO1, pMan->SO1, pMan->RO1, option);
    // Bmatch_CalSupp(pNtk2, pMan->suppFunc2, pMan->FI2, pMan->FO2, pMan->SO2, pMan->RO2, option);

    // Structural support information
    Bmatch_CalStructSupp(pMan->SO1, pNtk1);
    Bmatch_CalStructSupp(pMan->SO2, pNtk2);

    // Functional support information
    Bmatch_CalFuncSupp(pMan->suppFunc1, pMan->FI1, pMan->FO1, pNtk1, option);
    Bmatch_CalFuncSupp(pMan->suppFunc2, pMan->FI2, pMan->FO2, pNtk2, option);

    // Redundant support information
    Bmatch_CalRedundSupp(pMan->RO1, pMan->SO1, pMan->FO1);
    Bmatch_CalRedundSupp(pMan->RO2, pMan->SO2, pMan->FO2);

    // Sensitivity or others
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

void Bmatch_CalFuncSupp(vSupp &pVec, vSense &iSen, vSense &oFuncSupp, Abc_Ntk_t *pNtk, int option) {
    int i, j, k;
    int suppFunc;
    Vec_Ptr_t *vResult;

    pVec.reserve(Abc_NtkPoNum(pNtk));
    iSen.resize(Abc_NtkPiNum(pNtk));
    oFuncSupp.reserve(Abc_NtkPoNum(pNtk));
    vResult = Sim_ComputeFunSupp(pNtk, (option & VERBOSE_DETAIL_MASK) != 0);
    for (i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
        std::set<int> sen;
        suppFunc = 0;
        for (j = 0, k = (Abc_NtkPiNum(pNtk) / 32) + (Abc_NtkPiNum(pNtk) % 32 != 0); j < k; ++j) {
            suppFunc += __builtin_popcount(((int*)vResult->pArray[i])[j]);

            if (option & VERBOSE_MASK) Abc_Print(1, "%s(%d): ", Abc_ObjName(Abc_NtkPo(pNtk, i)), i);
            unsigned n = ((unsigned*)vResult->pArray[i])[j];
            for(int l = 0, num = (j == (Abc_NtkPiNum(pNtk) / 32)) ? ((Abc_NtkPiNum(pNtk) % 32 == 0) ? 32 : (Abc_NtkPiNum(pNtk) % 32)) : 32; l < num; ++l) {
                if (n & 1) {
                    if (option & VERBOSE_MASK) Abc_Print(1, "%s, ", Abc_ObjName(Abc_NtkPi(pNtk, j * 32 + l)));
                    iSen[j * 32 + l].insert(i);
                    sen.insert(j * 32 + l);
                }
                n >>= 1;
            }
        }
        if (option & VERBOSE_MASK) Abc_Print(1, "suppFunc = %d\n", suppFunc);
        pVec.emplace_back(i, suppFunc);
        oFuncSupp.emplace_back(std::move(sen));
    }

    std::sort(pVec.begin(), pVec.end(), [](std::pair<int, int> &t1, std::pair<int, int> &t2) { return t1.second < t2.second; });

    Vec_PtrFree(vResult);
}

void Bmatch_CalStructSupp(vSense &oStrSupp, Abc_Ntk_t *pNtk) {
    Vec_Ptr_t *vSupp, *vNodes;
    Abc_Obj_t *pObj, *pObjTemp;
    int i, j, index;

    oStrSupp.resize(Abc_NtkPoNum(pNtk));
    Abc_NtkForEachPo(pNtk, pObj, i) {
        vSupp = Abc_NtkNodeSupport(pNtk, &pObj, 1);
        vNodes = Abc_NtkDfsNodes(pNtk, &pObj, 1);
        Vec_PtrForEachEntry(Abc_Obj_t *, vSupp, pObjTemp, j) {
            if (Abc_ObjIsPi(pObjTemp)) {
                index = Bmatch_LinearSearchPiName2Index(pNtk, pObjTemp);
                if (index != -1) oStrSupp[i].insert(index);
            }
        }
        Vec_PtrFree(vSupp);
        Vec_PtrFree(vNodes);
    }
    Abc_NtkCleanMarkA(pNtk);
}

void Bmatch_CalRedundSupp(vSense &rSupp, vSense &oStrSupp, vSense &oFuncSupp) {
    // StrSupp - FuncSupp
    assert(oStrSupp.size() == oFuncSupp.size());
    rSupp.resize(oStrSupp.size());

    for (int i = 0; i < oStrSupp.size(); ++i) {
        for (auto &p : oStrSupp[i]) {
            if (oFuncSupp[i].count(p) == 0) rSupp[i].insert(p);
        }
    }
}

void Bmatch_CalSupp(Abc_Ntk_t *pNtk, vSupp &vSupp, vSense &iFuncSupp, vSense &oFuncSupp, vSense &oStrSupp, vSense &oRedundSupp, int option) {
    int i, j, k;
    int suppFunc;
    Vec_Ptr_t *vResultFun, *vResultStr;

    iFuncSupp.resize(Abc_NtkPiNum(pNtk));
    oFuncSupp.reserve(Abc_NtkPoNum(pNtk));
    oStrSupp.resize(Abc_NtkPoNum(pNtk));
    oRedundSupp.resize(Abc_NtkPoNum(pNtk));

    vResultFun = Sim_ComputeFunSupp(pNtk, (option & VERBOSE_DETAIL_MASK) != 0);
    vResultStr = Sim_ComputeStrSupp(pNtk);
    for (i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
        std::set<int> oFunc;
        suppFunc = 0;
        for (j = 0, k = (Abc_NtkPiNum(pNtk) / 32) + (Abc_NtkPiNum(pNtk) % 32 != 0); j < k; ++j) {
            suppFunc += __builtin_popcount(((int*)vResultFun->pArray[i])[j]);

            if (option & VERBOSE_MASK) Abc_Print(1, "%s(%d): ", Abc_ObjName(Abc_NtkPo(pNtk, i)), i);
            unsigned f = ((unsigned*)vResultFun->pArray[i])[j];
            unsigned s = ((unsigned*)vResultStr->pArray[i])[j];
            printf("hahahah:::  %u", s);
            unsigned r = s & ~f;
            for(int l = 0, num = (j == (Abc_NtkPiNum(pNtk) / 32)) ? ((Abc_NtkPiNum(pNtk) % 32 == 0) ? 32 : (Abc_NtkPiNum(pNtk) % 32)) : 32; l < num; ++l) {
                if (f & 1) {
                    if (option & VERBOSE_MASK) Abc_Print(1, "%s", Abc_ObjName(Abc_NtkPi(pNtk, j * 32 + l)));
                    iFuncSupp[j * 32 + l].insert(i);
                    oFunc.insert(j * 32 + l);
                }
                if (s & 1) {
                    oStrSupp[i].insert(j * 32 + l);
                }
                if (r & 1) {
                    if (option & VERBOSE_MASK) Abc_Print(1, "(X), ");
                    oRedundSupp[i].insert(j * 32 + l);
                } else {
                    if (option & VERBOSE_MASK) Abc_Print(1, "(V), ");
                }
                f >>= 1;
                s >>= 1;
                r >>= 1;
            }
        }
        if (option & VERBOSE_MASK) Abc_Print(1, "suppFunc = %d\n", suppFunc);
        vSupp.emplace_back(i, suppFunc);
        oFuncSupp.emplace_back(std::move(oFunc));
    }
    
    std::sort(vSupp.begin(), vSupp.end(), [](std::pair<int, int> &t1, std::pair<int, int> &t2) { return t1.second < t2.second; });

    Vec_PtrFree(vResultFun);
    Vec_PtrFree(vResultStr);
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