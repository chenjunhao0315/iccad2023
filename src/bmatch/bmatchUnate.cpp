#include "bmatch.hpp"

#include <math.h>

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

int *Bmatch_SimPatternParallel(Abc_Ntk_t *pNtk, int *pModel);
void Bmatch_RandomPattern(int *pModel, int Pi);
void Bmatch_FlipSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int *pModel, int *pValueGolden, int Pi);
void Bmatch_MonoSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int *pModel);
void Bmatch_RandomSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int iter);

void Bmatch_RandomAlterPattern(int *pModel, int Pi);
void Bmatch_RandomExtendPattern(int *pModel, int Pi);
void Bmatch_CounterExampleSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int *pModel, int iter);
void Bmatch_SatUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int Po, int iter);
sat_solver *Bmatch_UnateCheckerSatSolver(Abc_Ntk_t *pNtk, int Po, int *inputs, int *inputControls, int *outputControl);
Abc_Ntk_t *Bmatch_UnateChecker(Abc_Ntk_t *pNtk, int Po, Abc_Obj_t **inputs, Abc_Obj_t **inputControls, Abc_Obj_t **outputContorls);

#ifdef __cplusplus
}
#endif

void Bmatch_SatUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int Po, int iter) {
    int trival_case = 0;
    int *inputs = ABC_ALLOC(int, Abc_NtkPiNum(pNtk));
    int *inputControls = ABC_ALLOC(int, Abc_NtkPiNum(pNtk));
    int outputControl;
    lit *controls = ABC_ALLOC(lit, Abc_NtkPiNum(pNtk) + 5);
    int *pModel = ABC_ALLOC(int, Abc_NtkPiNum(pNtk));

    sat_solver *pSolver = Bmatch_UnateCheckerSatSolver(pNtk, Po, inputs, inputControls, &outputControl);
    lbool status;

    if (!pSolver) trival_case = 1;
    if (pSolver && sat_solver_simplify(pSolver) == 0) trival_case = 1;

    if (!trival_case) {
        for (int i = 0; i < Abc_NtkPiNum(pNtk); ++i) {
            controls[i] = toLitCond(inputControls[i], 1);
        }

        for (int i = 0; i < Abc_NtkPiNum(pNtk); ++i) {
            // not structural support
            if (inputs[i] < 0) {
                unateMat[Po][i] = 0;
                continue;
            }

            int unateness = unateMat[Po][i];
            if (unateness == 3 || unateness == -1) continue;

            int index = Abc_NtkPiNum(pNtk);
            // change the i-th input into constant 0, 1
            controls[i] = lit_neg(controls[i]);
            controls[index++] = toLitCond(inputs[i], 1);

            if (!unateness) {
                status = sat_solver_solve(pSolver, controls, controls + index, 0, 0, 0, 0);
                if (status == l_True) {
                    // SAT -> dependent
                    // check the value of the output control
                    // (0 -> negative unate, 1 -> positive unate)
                    int oc = sat_solver_var_value(pSolver, outputControl);
                    unateness = 1 << (1 - oc);

                    // get the SAT model
                    for (int j = 0; j < Abc_NtkPiNum(pNtk); ++j) {
                        // if inputs[j] < 0 (notice that inputs[j] != 0)
                        // the input is not a structural support
                        if (inputs[j] > 0)
                            pModel[j] = sat_solver_var_value(pSolver, inputs[j]);
                    }
                    Bmatch_CounterExampleSimUnate(pNtk, unateMat, pModel, 10);
                }
            }

            if (unateness > 0 && unateness < 3) {
                // get the value of the output control
                // (0 -> negative unate, 1 -> positive unate)
                int oc = unateness & 1;
                // set the output control to the opposite
                // and check the sat again
                controls[index++] = toLit(outputControl) ^ oc;
                status = sat_solver_solve(pSolver, controls, controls + index, 0, 0, 0, 0);

                if (status == l_True) {
                    // SAT -> binate
                    unateness = 3;
                    // get the SAT model
                    for (int j = 0; j < Abc_NtkPiNum(pNtk); ++j) {
                        // if inputs[j] is NULL
                        // the input is not a structural support
                        if (inputs[j] > 0)
                            pModel[j] = sat_solver_var_value(pSolver, inputs[j]);
                    }
                    Bmatch_CounterExampleSimUnate(pNtk, unateMat, pModel, iter);
                }
            }
            unateMat[Po][i] = unateness;
            controls[i] = lit_neg(controls[i]);
        }
    }
    
    if (pSolver) sat_solver_delete(pSolver);
    ABC_FREE(inputs);
    ABC_FREE(inputControls);
    ABC_FREE(controls);
    ABC_FREE(pModel);

    for (int i = 0; i < Abc_NtkPiNum(pNtk); ++i) {
        if (unateMat[Po][i] == 0)
            unateMat[Po][i] = -1;
    }
}

static inline int addClause(sat_solver * pSat, Vec_Int_t * vLits) {
    return sat_solver_addclause(pSat, vLits->pArray, vLits->pArray + vLits->nSize);
}

static inline int addAigClause(sat_solver * pSat, Abc_Obj_t * pObj, Vec_Int_t * vLits) {
    assert(!Abc_ObjIsComplement(pObj));
    // c = ab -> (a' + b' + c)(a + c')(b + c')
    lit a = toLitCond(Abc_ObjFanin0(pObj)->iTemp, Abc_ObjFaninC0(pObj));
    lit b = toLitCond(Abc_ObjFanin1(pObj)->iTemp, Abc_ObjFaninC1(pObj));
    lit c = toLit(pObj->iTemp);
    // (a' + b' + c)
    vLits->nSize = 0;
    Vec_IntPush(vLits, lit_neg(a));
    Vec_IntPush(vLits, lit_neg(b));
    Vec_IntPush(vLits, c);
    if (!addClause(pSat, vLits))
        return 0;
    // (b + c')
    vLits->nSize = 0;
    Vec_IntPush(vLits, b);
    Vec_IntPush(vLits, lit_neg(c));
    if (!addClause(pSat, vLits))
        return 0;
    // (a + c')
    vLits->nSize = 0;
    Vec_IntPush(vLits, a);
    Vec_IntPush(vLits, lit_neg(c));
    if (!addClause(pSat, vLits))
        return 0;
    return 1;
}

sat_solver *Bmatch_UnateCheckerSatSolver(Abc_Ntk_t *pNtk, int Po, int *inputs, int *inputControls, int *outputControl) {
    int i;
    Abc_Obj_t *pObj;
    
    Abc_Obj_t **inputsObj = ABC_ALLOC(Abc_Obj_t *, Abc_NtkPiNum(pNtk));
    Abc_Obj_t **inputControlsObj = ABC_ALLOC(Abc_Obj_t *, Abc_NtkPiNum(pNtk));
    Abc_Obj_t *outputControlObj;

    Abc_Ntk_t *pChecker = Bmatch_UnateChecker(pNtk, Po, inputsObj, inputControlsObj, &outputControlObj);

    sat_solver *pSolver = sat_solver_new();

    int cnt = 0;
    Abc_AigConst1(pChecker)->iTemp = cnt++;
    Abc_NtkForEachPi(pChecker, pObj, i) {
        pObj->iTemp = cnt++;
    }
    Abc_NtkForEachPo(pChecker, pObj, i) {
        pObj->iTemp = cnt++;
    }
    Abc_AigForEachAnd(pChecker, pObj, i) {
        pObj->iTemp = cnt++;
    }

    for (int i = 0; i < Abc_NtkPiNum(pNtk); ++i) {
        inputs[i]        = inputsObj[i] ? inputsObj[i]->iTemp : -1;
        inputControls[i] = inputControlsObj[i]->iTemp;
    }
    outputControl[0] = outputControlObj->iTemp;

    Vec_Int_t *vLits = Vec_IntAlloc(100);

    vLits->nSize = 0;
    // const
    Vec_IntPush(vLits, toLit(Abc_AigConst1(pChecker)->iTemp));
    if (!addClause(pSolver, vLits))
        return NULL;

    // and clause
    Abc_AigForEachAnd(pChecker, pObj, i) {
        if (!addAigClause(pSolver, pObj, vLits))
            return NULL;
    }

    // output clause
    vLits->nSize = 0;
    Vec_IntPush(vLits, toLitCond(Abc_ObjFanin0(Abc_NtkPo(pChecker, 0))->iTemp, Abc_ObjFaninC0(Abc_NtkPo(pChecker, 0))));
    if (!addClause(pSolver, vLits))
        return NULL;

    // Vec_IntFree(vLits);
    Abc_NtkDelete(pChecker);
    ABC_FREE(inputsObj);
    ABC_FREE(inputControlsObj);

    return pSolver;
}

Abc_Ntk_t *Bmatch_UnateChecker(Abc_Ntk_t *pNtk, int Po, Abc_Obj_t **inputs, Abc_Obj_t **inputControls, Abc_Obj_t **outputControls) {
    int i;
    Abc_Obj_t *pObj;

    int inv = Abc_ObjFaninC0(Abc_NtkPo(pNtk, Po));
    
    char Buffer[100];
    sprintf(Buffer, "%d_cone", Po);
    Abc_Ntk_t *pCone = Abc_NtkCreateCone(pNtk, Abc_ObjFanin0(Abc_NtkPo(pNtk, Po)), Buffer, 1);
    Abc_Ntk_t *pChecker =Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);

    Abc_AigConst1(pCone)->pCopy = Abc_AigConst1(pChecker);
    // create control signal
    Abc_NtkForEachPi(pCone, pObj, i) {
        inputControls[i] = Abc_NtkCreatePi(pChecker);
    }
    outputControls[0] = Abc_NtkCreatePi(pChecker);

    // golden path
    Abc_NtkForEachPi(pCone, pObj, i) {
        pObj->pCopy = inputs[i] = Abc_NtkCreatePi(pChecker);
    }
    Abc_AigForEachAnd(pCone, pObj, i) {
        pObj->pCopy = Abc_AigAnd((Abc_Aig_t *)(pChecker->pManFunc), Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj));
    }
    Abc_Obj_t *pGolden = Abc_AigXor((Abc_Aig_t *)(pChecker->pManFunc), Abc_ObjChild0Copy(Abc_NtkPo(pCone, 0)), Abc_ObjNotCond(outputControls[0], inv));

    // flip path
    Abc_NtkForEachPi(pCone, pObj, i) {
        pObj->pCopy = Abc_AigXor((Abc_Aig_t *)(pChecker->pManFunc), inputs[i], inputControls[i]);
    }
    Abc_AigForEachAnd(pCone, pObj, i) {
        pObj->pCopy = Abc_AigAnd((Abc_Aig_t *)(pChecker->pManFunc), Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj));
    }
    Abc_Obj_t *pFlip = Abc_AigXor((Abc_Aig_t *)(pChecker->pManFunc), Abc_ObjChild0Copy(Abc_NtkPo(pCone, 0)), Abc_ObjNotCond(outputControls[0], inv ^ 1));

    // checker
    Abc_Obj_t *pCheck = Abc_AigAnd((Abc_Aig_t *)(pChecker->pManFunc), pGolden, pFlip);
    Abc_ObjAddFanin(Abc_NtkCreatePo(pChecker), pCheck);

    Abc_NtkForEachPi(pCone, pObj, i) {
        if (Vec_IntSize(&(pObj->vFanouts)) == 0)
            inputs[i] = NULL;
    }

    Abc_NtkDelete(pCone);

    return pChecker;
}

int *Bmatch_SimPatternParallel(Abc_Ntk_t *pNtk, int *pModel) {
    assert(Abc_NtkIsStrash(pNtk));
    int i, v0, v1;
    Abc_Obj_t *pObj;

    Abc_AigConst1(pNtk)->iTemp = S_ONE;
    Abc_NtkForEachCi(pNtk, pObj, i) {
        pObj->iTemp = pModel[i];
    }
    Abc_NtkForEachNode(pNtk, pObj, i) {
        v0 = Abc_ObjFanin0(pObj)->iTemp ^ ((Abc_ObjFaninC0(pObj)) ? S_ONE : 0);
        v1 = Abc_ObjFanin1(pObj)->iTemp ^ ((Abc_ObjFaninC1(pObj)) ? S_ONE : 0);
        pObj->iTemp = v0 & v1;
    }
    int *pValue = ABC_ALLOC(int, Abc_NtkPoNum(pNtk));
    Abc_NtkForEachPo(pNtk, pObj, i) {
        pValue[i] = Abc_ObjFanin0(pObj)->iTemp ^ ((Abc_ObjFaninC0(pObj)) ? S_ONE : 0);
    }

    return pValue;
}

void Bmatch_CounterExampleSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int *pModel, int iter) {
    Bmatch_RandomExtendPattern(pModel, Abc_NtkPiNum(pNtk));
    Bmatch_MonoSimUnate(pNtk, unateMat, pModel);
    for (int i = 0; i < iter; ++i) {
        Bmatch_RandomAlterPattern(pModel, Abc_NtkPiNum(pNtk));
        Bmatch_MonoSimUnate(pNtk, unateMat, pModel);
    }
}

void Bmatch_FlipSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int *pModel, int *pValueGolden, int Pi) {
    pModel[Pi] ^= S_ONE;
    int *pValue = Bmatch_SimPatternParallel(pNtk, pModel);
    pModel[Pi] ^= S_ONE;

    int unateness;
    for (int i = 0; i < Abc_NtkPoNum(pNtk); ++i) {
        pValueGolden[i] ^= pModel[Pi];
        pValue[i]       ^= pModel[Pi];
        if ((~pValueGolden[i]) & (pValue[i])) { // postive unate
            unateMat[i][Pi] |= 1;
        }
        if ((pValueGolden[i]) & (~pValue[i])) { // negative unate
            unateMat[i][Pi] |= 2;
        }
        pValueGolden[i] ^= pModel[Pi];
        pValue[i]       ^= pModel[Pi];
    }

    ABC_FREE(pValue);
}

void Bmatch_MonoSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int *pModel) {
    int *pValueGolden = Bmatch_SimPatternParallel(pNtk, pModel);
    
    for (int i = 0; i < Abc_NtkPiNum(pNtk); ++i) {
        Bmatch_FlipSimUnate(pNtk, unateMat, pModel, pValueGolden, i);
    }

    ABC_FREE(pValueGolden);
}

void Bmatch_RandomSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int iter) {
    if (!unateMat.empty()) return;
    int *pModel = ABC_ALLOC(int, Abc_NtkPiNum(pNtk));

    unateMat = std::vector<std::vector<int> >(Abc_NtkPoNum(pNtk), std::vector<int>(Abc_NtkPiNum(pNtk), 0));
    for (int i = 0; i < iter; ++i) {
        Bmatch_RandomPattern(pModel, Abc_NtkPiNum(pNtk));
        Bmatch_MonoSimUnate(pNtk, unateMat, pModel);
    }

    ABC_FREE(pModel);
}

void Bmatch_RandomAlterPattern(int *pModel, int Pi) {
    assert(pModel);
    for (int i = 0; i < sqrt(Pi); ++i) {
        pModel[rand() % Pi] ^= rand();
    }
}

void Bmatch_RandomExtendPattern(int *pModel, int Pi) {
    assert(pModel);
    for (int i = 0; i < Pi; ++i) {
        pModel[i] = (pModel[i]) ? S_ONE : S_ZERO;
    }
    for (int i = 0; i < sqrt(Pi); ++i) {
        pModel[rand() % Pi] ^= ~(rand() | 1);
    }
}

void Bmatch_RandomPattern(int *pModel, int Pi) {
    assert(pModel);
    for (int i = 0; i < Pi; ++i) {
        pModel[i] = rand();
    }
}

ABC_NAMESPACE_IMPL_END