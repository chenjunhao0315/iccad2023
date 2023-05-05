#include "bmatch.hpp"

#include "proof/fra/fra.h"
#include "aig/aig/aig.h"
#include "proof/fraig/fraig.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

extern Aig_Man_t *Abc_NtkToDar( Abc_Ntk_t *pNtk, int fExors, int fRegisters );
static void Bmatch_NtkVerifyReportError(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel, vMatch &MI, vMatch &MO);
int Bmatch_SatFraig(Abc_Ntk_t **ppNtk);
EcResult Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int fVerbose);

#ifdef __cplusplus
}
#endif

EcResult Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int fVerbose) {
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
        Bmatch_NtkVerifyReportError(pNtk1, pNtk2, pNtkMiter->pModel, MI, MO);
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

    RetValue = Bmatch_SatFraig(&pNtkMiter);
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
        Bmatch_NtkVerifyReportError(pNtk1, pNtk2, pNtkMiter->pModel, MI, MO);
    }
    model = pNtkMiter->pModel; pNtkMiter->pModel = NULL;
    Abc_NtkDelete(pNtkMiter);

    return {Status, model};
}

int Bmatch_SatFraig(Abc_Ntk_t **ppNtk) {
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
    RetValue = Fra_FraigSat(pMan, (ABC_INT64_T)5000, (ABC_INT64_T)0, 0, 0, 0, 1, 0, 0, 0); 
    pNtk->pModel = (int *)pMan->pData, pMan->pData = NULL;
    Aig_ManStop(pMan);

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