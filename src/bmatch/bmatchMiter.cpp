#include "bmatch.hpp"

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

// MI/MO
// ckt1[0] = {ckt2[?], ckt2[?]}
// ckt1[1] = {ckt2[?]}

Abc_Ntk_t *Bmatch_NtkMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
int  Bmatch_NtkMiterCheck(vMatch &MI, vMatch &MO, Abc_Ntk_t *pNtk2);
void Bmatch_NtkMiterPrepare(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MI);
void Bmatch_NtkMiterAddOne(Abc_Ntk_t *pNtk, Abc_Ntk_t *pNtkMiter);
void Bmatch_NtkMiterFinalize(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MO);

#ifdef __cplusplus
}
#endif

Abc_Ntk_t *Bmatch_NtkMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO) {
    char Buffer[1000];
    Abc_Ntk_t *pNtkMiter;

    if (!Bmatch_NtkMiterCheck(MI, MO, pNtk2)) return NULL;

    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    sprintf(Buffer, "%s_%s_miter", Abc_NtkName(pNtk1), Abc_NtkName(pNtk2));
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));

    // Prepare
    Bmatch_NtkMiterPrepare(pNtk1, pNtk2, pNtkMiter, MI);
    // Construct
    Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
    Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);
    // Finalize
    Bmatch_NtkMiterFinalize(pNtk1, pNtk2, pNtkMiter, MO);
    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    if (!Abc_NtkCheck(pNtkMiter)) {
        Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
        Abc_NtkDelete(pNtkMiter);

        return NULL;
    }

    return pNtkMiter;
}

int  Bmatch_NtkMiterCheck(vMatch &MI, vMatch &MO, Abc_Ntk_t *pNtk2) {
    int MIs = 0, MOs = 0;
    for (auto &i : MI)
        MIs += i.size();
    for (auto &o : MO)
        if ((MOs = o.size()) != 0) break;

    return MIs && MOs && MIs == Abc_NtkPiNum(pNtk2);
}

void Bmatch_NtkMiterPrepare(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MI) {
    int i, start = 0;
    Abc_Obj_t *pObj, *pObjNew;
    char buffer[1000];

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    Abc_NtkForEachCi(pNtk1, pObj, i) {
        start = 0;
        memset(buffer, 0, sizeof(buffer));
        pObjNew = Abc_NtkCreatePi(pNtkMiter);

        pObj->pCopy = pObjNew;
        for (auto &p : MI[i]) {
            pObj = Abc_NtkCi(pNtk2, p.var());
            pObj->pCopy = (p.sign()) ? Abc_ObjNot(pObjNew) : pObjNew;
            if (p.sign()) strcat(buffer, "~");
            strcat(buffer, Abc_ObjName(pObj));
            if (start++ != MI[i].size() - 1) strcat(buffer, "_");
        }
        if (strlen(buffer) == 0) sprintf(buffer, "non_map_%d", i);
        Abc_ObjAssignName(pObjNew, buffer, NULL);
    }
    // Test for const input (not really sure if it works or not)
    for (auto &p : MI.back()) {
        pObj = Abc_NtkCi(pNtk2, p.var());
        pObj->pCopy = (p.sign()) ? Abc_ObjNot(Abc_AigConst1(pNtkMiter)) : Abc_AigConst1(pNtkMiter);
    }
    pObjNew = Abc_NtkCreatePo(pNtkMiter);
    Abc_ObjAssignName(pObjNew, "miter", Abc_ObjName(pObjNew));
}

void Bmatch_NtkMiterAddOne(Abc_Ntk_t *pNtk, Abc_Ntk_t *pNtkMiter) {
    Abc_Obj_t * pNode;
    int i;
    assert(Abc_NtkIsDfsOrdered(pNtk));
    Abc_AigForEachAnd(pNtk, pNode, i)
        pNode->pCopy = Abc_AigAnd((Abc_Aig_t *)pNtkMiter->pManFunc, Abc_ObjChild0Copy(pNode), Abc_ObjChild1Copy(pNode));
}

void Bmatch_NtkMiterFinalize(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MO) {
    Vec_Ptr_t *vPairs;
    Abc_Obj_t *pMiter, *pNode;
    int i, j;

    vPairs = Vec_PtrAlloc(100);

    Abc_NtkForEachCo(pNtk1, pNode, i) {
        if (!MO[i].empty()) {
            Vec_PtrPush(vPairs, Abc_ObjChild0Copy(pNode));
            for (j = 0; j < MO[i].size(); ++j) {
                pNode = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
                pNode = (MO[i][j].sign()) ? Abc_ObjNot(pNode) : pNode;
                Vec_PtrPush(vPairs, pNode);
            }
        }
    }

    pMiter = Abc_AigMiter((Abc_Aig_t *)pNtkMiter->pManFunc, vPairs, 0);
    Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pMiter);
    
    Vec_PtrFree(vPairs);
}

ABC_NAMESPACE_IMPL_END