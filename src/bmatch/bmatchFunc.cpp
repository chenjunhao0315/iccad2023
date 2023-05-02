#include <algorithm>

#include "bmatch.hpp"
#include "opt/sim/sim.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_BusNameMaping(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_CalStructSupp(vSense &oSupp, Abc_Ntk_t *pNtk);
void Bmatch_CalFuncSupp(vSupp &pVec, vSense &iSen, vSense &oSen, Abc_Ntk_t *pNtk, int option);
void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

#ifdef __cplusplus
}
#endif

void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    // Setup network (Resythn or something)
    Bmatch_NtkSetup(pNtk1, pNtk2, option);

    // Convert bus information from string to index
    Bmatch_BusNameMaping(pMan, pNtk1, pNtk2);

    // Functional support information
    Bmatch_CalFuncSupp(pMan->suppFunc1, pMan->FI1, pMan->FO1, pNtk1, option);
    Bmatch_CalFuncSupp(pMan->suppFunc2, pMan->FI2, pMan->FO2, pNtk2, option);

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

void Bmatch_CalFuncSupp(vSupp &pVec, vSense &iSen, vSense &oSen, Abc_Ntk_t *pNtk, int option) {
    int i, j, k;
    int suppFunc;
    Vec_Ptr_t *vResult;

    pVec.reserve(Abc_NtkPoNum(pNtk));
    iSen.resize(Abc_NtkPiNum(pNtk));
    oSen.reserve(Abc_NtkPoNum(pNtk));
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
        oSen.emplace_back(std::move(sen));
    }
    std::sort(pVec.begin(), pVec.end(), [](std::pair<int, int> &t1, std::pair<int, int> &t2) { return t1.second < t2.second; });
}

void Bmatch_CalStructSupp(vSense &oSupp, Abc_Ntk_t *pNtk) {
    Vec_Ptr_t *vSupp, *vNodes;
    Abc_Obj_t *pObj, *pObjTemp;
    int i, j;

    oSupp.resize(Abc_NtkPoNum(pNtk));
    Abc_NtkForEachPo(pNtk, pObj, i) {
        vSupp = Abc_NtkNodeSupport(pNtk, &pObj, 1);
        vNodes = Abc_NtkDfsNodes(pNtk, &pObj, 1);
        Vec_PtrForEachEntry(Abc_Obj_t *, vSupp, pObjTemp, j) {
            if (Abc_ObjIsPi(pObjTemp)) {
                
            }
        }
    }
}

ABC_NAMESPACE_IMPL_END