#include <stdio.h>

#include "bmatch.hpp"
#include "opt/sim/sim.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_CalFuncSupp(vSupp &pVec, Abc_Ntk_t *pNtk, int option);
void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

#ifdef __cplusplus
}
#endif

void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    Bmatch_NtkSetup(pNtk1, pNtk2, option);

    // Functional support information
    Bmatch_CalFuncSupp(pMan->suppFunc1, pNtk1, option);
    Bmatch_CalFuncSupp(pMan->suppFunc2, pNtk2, option);

    // Sensitivity or others
}

void Bmatch_CalFuncSupp(vSupp &pVec, Abc_Ntk_t *pNtk, int option) {
    int i, j, k;
    int suppFunc;
    Vec_Ptr_t *vResult;

    pVec.reserve(Abc_NtkPoNum(pNtk));
    vResult = Sim_ComputeFunSupp(pNtk, (option & VERBOSE_DETAIL_MASK) != 0);
    for (i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
        suppFunc = 0;
        for (j = 0, k = (Abc_NtkPiNum(pNtk) / 32) + (Abc_NtkPiNum(pNtk) % 32 != 0); j < k; ++j) {
            suppFunc += __builtin_popcount(((int*)vResult->pArray[i])[j]);

            if (option & VERBOSE_MASK) {
                Abc_Print(1, "%s(%d): ", Abc_ObjName(Abc_NtkPo(pNtk, i)), i);
                unsigned n = ((unsigned*)vResult->pArray[i])[j];
                for(int l = 0, num = (j == (Abc_NtkPiNum(pNtk) / 32)) ? ((Abc_NtkPiNum(pNtk) % 32 == 0) ? 32 : (Abc_NtkPiNum(pNtk) % 32)) : 32; l < num; ++l) {
                    if (n & 1) {
                        Abc_Print(1, "%s, ", Abc_ObjName(Abc_NtkPi(pNtk, j * 32 + l)));
                    }
                    n >>= 1;
                }
            }
        }
        if (option & VERBOSE_MASK) Abc_Print(1, "suppFunc = %d\n", suppFunc);
        pVec.emplace_back(i, suppFunc);
    }
    std::sort(pVec.begin(), pVec.end(), [](std::pair<int, int> &t1, std::pair<int, int> &t2) { return t1.second < t2.second; });
}

ABC_NAMESPACE_IMPL_END