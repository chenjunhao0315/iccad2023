#include "bmatch.hpp"

#include "print.hpp"

#include "base/io/ioAbc.h"
#include "aig/gia/giaAig.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

// MI/MO
// ckt1[0] = {ckt2[?], ckt2[?]}
// ckt1[1] = {ckt2[?]}

/*=== base/abci/abcDar.c ==============================================*/
extern Aig_Man_t * Abc_NtkToDar( Abc_Ntk_t * pNtk, int fExors, int fRegisters );
extern Abc_Ntk_t * Abc_NtkDC2( Abc_Ntk_t * pNtk, int fBalance, int fUpdateLevel, int fFanout, int fPower, int fVerbose );

/*=== aig/gia/giaAig.c ================================================*/
extern Gia_Man_t * Gia_ManFromAig( Aig_Man_t * p );

Abc_Ntk_t *Bmatch_NtkMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
Abc_Ntk_t *Bmatch_NtkMiterGia(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
Abc_Ntk_t *Bmatch_NtkControllableInputMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, int inv);
Abc_Ntk_t *Bmatch_NtkControllableInputOutputMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int  Bmatch_NtkMiterCheck(vMatch &MI, vMatch &MO, Abc_Ntk_t *pNtk2);
void Bmatch_NtkMiterPrepare(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MI, vMatch &MO);
void Bmatch_NtkMiterAddOne(Abc_Ntk_t *pNtk, Abc_Ntk_t *pNtkMiter);
void Bmatch_NtkMiterFinalize(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MO);

Abc_Obj_t *Bmatch_NtkCreateOr(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal);
Abc_Obj_t *Bmatch_NtkCreateAnd(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal);
Abc_Obj_t *Bmatch_NtkCreateParallelCase(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal);
Abc_Obj_t *Bmatch_NtkCreateMultiplexer(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal, int dontApplyNot);

#ifdef __cplusplus
}
#endif

Abc_Ntk_t *Bmatch_NtkMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO) {
    char Buffer[100];
    Abc_Ntk_t *pNtkMiter;

    if (!Bmatch_NtkMiterCheck(MI, MO, pNtk2)) return NULL;

    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    sprintf(Buffer, "%s_%s_miter", Abc_NtkName(pNtk1), Abc_NtkName(pNtk2));
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));

    // Prepare
    Bmatch_NtkMiterPrepare(pNtk1, pNtk2, pNtkMiter, MI, MO);
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

void Bmatch_NtkMiterPrepare(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Abc_Ntk_t *pNtkMiter, vMatch &MI, vMatch &MO) {
    int i, start = 0;
    Abc_Obj_t *pObj, *pObjNew;
    char buffer[100];

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
    Abc_ObjAssignName(pObjNew, "miter_", Abc_ObjName(pObjNew));
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
    Abc_Obj_t *pMiter, *pNode, *pNodeTemp;
    int i, j;

    vPairs = Vec_PtrAlloc(100);

    Abc_NtkForEachCo(pNtk1, pNode, i) {
        if (!MO[i].empty()) {
            pNodeTemp = pNode;
            for (j = 0; j < MO[i].size(); ++j) {
                // Move inside the loop (ensuring create a pair?)
                Vec_PtrPush(vPairs, Abc_ObjChild0Copy(pNodeTemp));
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

Abc_Ntk_t *Bmatch_NtkGiaMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO) {
    assert(Abc_NtkIsDfsOrdered(pNtk1));
    assert(Abc_NtkIsDfsOrdered(pNtk2));
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));
    Abc_Aig_t *pMan = (Abc_Aig_t *)pNtkMiter->pManFunc;
    int i;
    Abc_Obj_t *pObj, *pObjTemp, *pObjNew;

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // MI
    Abc_NtkForEachCi(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "_ntk1");
        Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), Buffer);
        pObj->pCopy = pObjNew;

        for (auto &p : MI[i]) {
            pObj = Abc_NtkPi(pNtk2, p.var());
            pObj->pCopy = (p.sign()) ? Abc_ObjNot(pObjNew) : pObjNew;
        }
    }
    for (auto &p : MI.back()) {
        pObj = Abc_NtkPi(pNtk2, p.var());
        pObj->pCopy = (p.sign()) ? Abc_ObjNot(Abc_AigConst1(pNtkMiter)) : Abc_AigConst1(pNtkMiter);
    }

    // MO
    int count = 0;
    for (auto &mo : MO) {
        for (int i = 0; i < mo.size(); ++i) {
            pObjNew = Abc_NtkCreatePo(pNtkMiter);
            sprintf(Buffer, "miter_cir1_");
            Abc_ObjAssignName(pObjNew, Buffer, (char*)std::to_string(count++).c_str());
            pObjNew = Abc_NtkCreatePo(pNtkMiter);
            sprintf(Buffer, "miter_cir2_");
            Abc_ObjAssignName(pObjNew, Buffer, (char*)std::to_string(count++).c_str());
        }
    }

    // Construct
    Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
    Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);

    // Miter
    count = 0;
    Abc_NtkForEachCo(pNtk1, pObj, i) {
        if (!MO[i].empty()) {
            pObjTemp = Abc_ObjChild0Copy(pObj);
            for (int j = 0; j < MO[i].size(); ++j) {
                pObj = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
                pObj = (MO[i][j].sign()) ? Abc_ObjNot(pObj) : pObj;
                Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, count++), pObjTemp);
                Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, count++), pObj);
            }
        }
    }
    assert(2 * count == Abc_NtkPoNum(pNtkMiter));

    // Cleanup
    Abc_AigCleanup(pMan);

    if (!Abc_NtkCheck(pNtkMiter)) {
        Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
        Abc_NtkDelete(pNtkMiter);

        return NULL;
    }

    return pNtkMiter;
}

Abc_Ntk_t *Bmatch_NtkControllableInputMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, int inv) {
    assert(Abc_NtkIsDfsOrdered(pNtk1));
    assert(Abc_NtkIsDfsOrdered(pNtk2));
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));
    Abc_Aig_t *pMan = (Abc_Aig_t *)pNtkMiter->pManFunc;
    int i;
    Abc_Obj_t *pObj, *pObjTemp, *pObjNew;

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // MI control
    std::vector<std::vector<Abc_Obj_t *> > controlPis;
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPi;
        int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));
        for (int j = 0; j < nControlPi; ++j) {
            pObj = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "controlPi_%d_%d", i, j);
            Abc_ObjAssignName(pObj, Buffer, NULL);
            controlPi.emplace_back(std::move(pObj));
        }
        pObj = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "controlPi_%d_%d_inv", i, Abc_NtkPiNum(pNtk1));
        Abc_ObjAssignName(pObj, Buffer, NULL);
        controlPi.emplace_back(std::move(pObj));

        controlPis.emplace_back(std::move(controlPi));
    }

    assert(controlPis.size() == Abc_NtkPiNum(pNtk2));

    // Ntk1 Pi
    std::vector<Abc_Obj_t *> Ntk1_Pis;
    Abc_NtkForEachPi(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "_ntk1");
        Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), Buffer);
        Ntk1_Pis.emplace_back(pObjNew);
        pObj->pCopy = pObjNew;
    }
    Ntk1_Pis.emplace_back(Abc_AigConst1(pNtkMiter));

    // PO
    pObjNew = Abc_NtkCreatePo(pNtkMiter);
    Abc_ObjAssignName(pObjNew, "miter", Abc_ObjName(pObjNew));

    // Create control input of Ntk2
    Abc_NtkForEachPi(pNtk2, pObj, i) {
        std::vector<Abc_Obj_t *> pControl = controlPis[i];
        pObjNew = Bmatch_NtkCreateMultiplexer(pMan, pControl, Ntk1_Pis, 0);
        pObj->pCopy = pObjNew;
    }

    // Construct over Ntk1 and Ntk2
    Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
    Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);

    // Finalize
    if (!inv) Bmatch_NtkMiterFinalize(pNtk1, pNtk2, pNtkMiter, MO);
    else {
        std::vector<Abc_Obj_t *> Xnors;
        Abc_NtkForEachCo(pNtk1, pObj, i) {
            if (!MO[i].empty()) {
                pObjTemp = Abc_ObjChild0Copy(pObj);
                for (int j = 0; j < MO[i].size(); ++j) {
                    pObj = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
                    pObj = (MO[i][j].sign()) ? Abc_ObjNot(pObj) : pObj;
                    pObjNew = Abc_ObjNot(Abc_AigXor(pMan, pObjTemp, pObj));
                    Xnors.push_back(pObjNew);
                }
            }
        }
        pObj = Bmatch_NtkCreateAnd(pMan, Xnors);
        Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pObj);
    }


    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    if (Abc_NtkNodeNum(pNtkMiter) > 500) {
        Abc_Ntk_t *pNtkTemp;
        pNtkMiter = Abc_NtkDC2(pNtkTemp = pNtkMiter, 1, 0, 1, 0, 0);
        Abc_NtkDelete(pNtkTemp);
    }

    if (!Abc_NtkCheck(pNtkMiter)) {
        Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
        Abc_NtkDelete(pNtkMiter);

        return NULL;
    }

    // Io_Write(pNtkMiter, "miter_rewrite.v", IO_FILE_VERILOG);

    return pNtkMiter;
}

Abc_Ntk_t *Bmatch_NtkControllableInputOutputMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    assert(Abc_NtkIsDfsOrdered(pNtk1));
    assert(Abc_NtkIsDfsOrdered(pNtk2));
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));
    Abc_Aig_t *pMan = (Abc_Aig_t *)pNtkMiter->pManFunc;
    int i, j;
    Abc_Obj_t *pObj, *pObjTemp, *pObjNew;
    Abc_Obj_t *pObjNtk1, *pObjNtk2;

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // MI control
    std::vector<std::vector<Abc_Obj_t *> > controlPis;
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPi;
        int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));
        for (int j = 0; j < nControlPi; ++j) {
            pObj = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "I_%d_%d", i, j);
            Abc_ObjAssignName(pObj, Buffer, NULL);
            controlPi.emplace_back(std::move(pObj));
        }
        pObj = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "I_%d_%d_inv", i, nControlPi);
        Abc_ObjAssignName(pObj, Buffer, NULL);
        controlPi.emplace_back(std::move(pObj));

        controlPis.emplace_back(std::move(controlPi));
    }

    assert(controlPis.size() == Abc_NtkPiNum(pNtk2));

    // MO control
    std::vector<std::vector<Abc_Obj_t *> > controlPos;
    for (int i = 0; i < Abc_NtkPoNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPo;
        int nControlPo = (int)(std::ceil(std::log2(Abc_NtkPoNum(pNtk1) + 1)));
        for (int j = 0; j < nControlPo; ++j) {
            pObj = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "O_%d_%d", i, j);
            Abc_ObjAssignName(pObj, Buffer, NULL);
            controlPo.emplace_back(std::move(pObj));
        }
        pObj = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "O_%d_inv", i);
        Abc_ObjAssignName(pObj, Buffer, NULL);
        controlPo.emplace_back(std::move(pObj));

        controlPos.emplace_back(std::move(controlPo));
    }

    // PO
    pObjNew = Abc_NtkCreatePo(pNtkMiter);
    Abc_ObjAssignName(pObjNew, "miter", Abc_ObjName(pObjNew));

    // Ntk1 Pi
    std::vector<Abc_Obj_t *> Ntk1_Pis;
    Abc_NtkForEachPi(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        sprintf(Buffer, "_ntk1");
        Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), Buffer);
        Ntk1_Pis.emplace_back(pObjNew);
        pObj->pCopy = pObjNew;
    }
    Ntk1_Pis.emplace_back(Abc_AigConst1(pNtkMiter));

    // Create control input of Ntk2
    Abc_NtkForEachPi(pNtk2, pObj, i) {
        std::vector<Abc_Obj_t *> pControl = controlPis[i];
        pObjNew = Bmatch_NtkCreateMultiplexer(pMan, pControl, Ntk1_Pis, 0);
        pObj->pCopy = pObjNew;
    }

    // Construct over Ntk1 and Ntk2
    Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
    Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);

    // Finalize
    // All possible miter
    std::vector<Abc_Obj_t *> Ntk_Miters;
    Abc_NtkForEachPo(pNtk2, pObjNtk2, i) {
        int nControlPo = (int)(std::ceil(std::log2(Abc_NtkPoNum(pNtk1) + 1)));
        std::vector<Abc_Obj_t *> Ntk_Xors;
        auto &controlPo = controlPos[i];
        Abc_Obj_t *pNtk2Po = Abc_ObjChild0Copy(pObjNtk2);
        Abc_Obj_t *pNtk2PoCInv = Abc_AigXor(pMan, controlPo[nControlPo], pNtk2Po);

        Abc_NtkForEachPo(pNtk1, pObjNtk1, j) {
            Abc_Obj_t *pNtk1Po = Abc_ObjChild0Copy(pObjNtk1);

            Ntk_Xors.push_back(Abc_AigXor(pMan, pNtk2PoCInv, pNtk1Po));
        }
        Ntk_Xors.push_back(Abc_ObjNot(Abc_AigConst1(pNtkMiter))); // nonmap always equivalent

        Abc_Obj_t *pFinal = Bmatch_NtkCreateMultiplexer(pMan, controlPo, Ntk_Xors, 1);
        Ntk_Miters.push_back(pFinal);
    }
    Abc_Obj_t *pMiter = Bmatch_NtkCreateOr(pMan, Ntk_Miters);
    Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pMiter);

    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    if (Abc_NtkNodeNum(pNtkMiter) > 500) {
        Abc_Ntk_t *pNtkTemp;
        pNtkMiter = Abc_NtkDC2(pNtkTemp = pNtkMiter, 1, 0, 1, 0, 0);
        Abc_NtkDelete(pNtkTemp);
    }

    if (!Abc_NtkCheck(pNtkMiter)) {
        Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
        Abc_NtkDelete(pNtkMiter);

        return NULL;
    }

    // Io_Write(pNtkMiter, "miter_dc2_balance.v", IO_FILE_VERILOG);

    return pNtkMiter;
}

Abc_Obj_t *Bmatch_NtkCreateMultiplexer(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal, int dontApplyNot) {
    int j;
    std::vector<Abc_Obj_t *> prevSignal = pSignal;
    std::vector<Abc_Obj_t *> currSignal;
    
    for (int i = 0; i < pControl.size() - 1; ++i) {
        Abc_Obj_t *pC = pControl[i];
        for (j = 0; j < prevSignal.size(); j += 2) {
            Abc_Obj_t *p0 = prevSignal[j];
            Abc_Obj_t *p1 = ((j + 1) == prevSignal.size()) ? prevSignal[j] : prevSignal[j + 1];
            currSignal.push_back(Abc_AigMux(pMan, pC, p1, p0));
        }
        prevSignal = currSignal;
        currSignal.clear();
    }

    Abc_Obj_t *pMux = prevSignal[0];
    Abc_Obj_t *pFinal = (dontApplyNot) ? pMux : Abc_AigXor(pMan, pControl.back(), pMux);

    return pFinal;
}

Abc_Obj_t *Bmatch_NtkCreateParallelCase(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal) {
    assert(pControl.size() == pSignal.size());

    std::vector<Abc_Obj_t *> pAnds;
    for (int i = 0; i < pControl.size(); ++i) {
        Abc_Obj_t *pObj = Abc_AigAnd(pMan, pControl[i], pSignal[i]);
        pAnds.emplace_back(pObj);
    }

    return Bmatch_NtkCreateOr(pMan, pAnds);
}

Abc_Obj_t *Bmatch_NtkCreateAndRec(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal, int start, int end) {
    Abc_Obj_t *pNode = NULL;
    int mid = (start + end) / 2;
    if (start == end) {
        pNode = pSignal[start];
    } else if (start < end) {
        Abc_Obj_t *pNode1 = Bmatch_NtkCreateAndRec(pMan, pSignal, start, mid);
        Abc_Obj_t *pNode2 = Bmatch_NtkCreateAndRec(pMan, pSignal, mid + 1, end);
        pNode = Abc_AigAnd(pMan, pNode1, pNode2);
    }

    return pNode;
}

Abc_Obj_t *Bmatch_NtkCreateAnd(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal) {
    return Bmatch_NtkCreateAndRec(pMan, pSignal, 0, pSignal.size() - 1);
}

Abc_Obj_t *Bmatch_NtkCreateOrRec(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal, int start, int end) {
    Abc_Obj_t *pNode = NULL;
    int mid = (start + end) / 2;
    if (start == end) {
        pNode = pSignal[start];
    } else if (start < end) {
        Abc_Obj_t *pNode1 = Bmatch_NtkCreateOrRec(pMan, pSignal, start, mid);
        Abc_Obj_t *pNode2 = Bmatch_NtkCreateOrRec(pMan, pSignal, mid + 1, end);
        pNode = Abc_AigOr(pMan, pNode1, pNode2);
    }

    return pNode;
}

Abc_Obj_t *Bmatch_NtkCreateOr(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal) {
    return Bmatch_NtkCreateOrRec(pMan, pSignal, 0, pSignal.size() - 1);
}

ABC_NAMESPACE_IMPL_END
