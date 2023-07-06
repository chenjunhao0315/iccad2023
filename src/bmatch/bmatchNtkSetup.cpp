#include "bmatch.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_NtkSetup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);
Abc_Ntk_t* Bmatch_NtkResynth(Abc_Ntk_t *pNtk);

extern Abc_Ntk_t * Abc_NtkDC2( Abc_Ntk_t * pNtk, int fBalance, int fUpdateLevel, int fFanout, int fPower, int fVerbose );

#ifdef __cplusplus
}
#endif

void Bmatch_NtkSetup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    Abc_NtkSetName(pNtk1, Extra_UtilStrsav("cir1"));
    Abc_NtkSetName(pNtk2, Extra_UtilStrsav("cir2"));

    if (option & RESYNTH_MASK) {
        pNtk1 = Bmatch_NtkResynth(pNtk1);
        pNtk2 = Bmatch_NtkResynth(pNtk2);
    }

    Abc_NtkOrderObjsByName(pNtk1, 0);
    Abc_NtkOrderObjsByName(pNtk2, 0);

    if (option & VERBOSE_MASK) Bmatch_NtkPrint(pNtk1);
    if (option & VERBOSE_DETAIL_MASK) Bmatch_NtkPrintIO(pNtk1);
    if (option & VERBOSE_MASK) Bmatch_NtkPrint(pNtk2);
    if (option & VERBOSE_DETAIL_MASK) Bmatch_NtkPrintIO(pNtk2);
}

Abc_Ntk_t* Bmatch_NtkResynth(Abc_Ntk_t *pNtk) {
    Abc_Ntk_t *pNtkTemp;

    // pNtk = Abc_NtkBalance(pNtkTemp = pNtk, 0, 0, 1);  Abc_NtkDelete(pNtkTemp);

    // Abc_NtkRewrite(pNtk, 1, 0, 0, 0, 0);
    // Abc_NtkRewrite(pNtk, 1, 1, 0, 0, 0);

    // pNtk = Abc_NtkBalance(pNtkTemp = pNtk, 0, 0, 1);  Abc_NtkDelete(pNtkTemp);

    // Abc_NtkRewrite(pNtk, 1, 1, 0, 0, 0);

    // pNtk = Abc_NtkBalance(pNtkTemp = pNtk, 0, 0, 1);  Abc_NtkDelete(pNtkTemp);

    pNtk = Abc_NtkStrash( pNtk, 0, 1, 0 );

    pNtk = Abc_NtkDC2( pNtk, 0, 0, 1, 0, 0 );

    return pNtk;
}

ABC_NAMESPACE_IMPL_END