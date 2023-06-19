#include "bmatch.hpp"
#include "bmatchQbf.hpp"
// #include "print.hpp"

#include "base/io/ioAbc.h"
#include "aig/gia/giaAig.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

/*=== base/abci/abcDar.c ==============================================*/
extern Aig_Man_t * Abc_NtkToDar( Abc_Ntk_t * pNtk, int fExors, int fRegisters );
extern Abc_Ntk_t * Abc_NtkDC2( Abc_Ntk_t * pNtk, int fBalance, int fUpdateLevel, int fFanout, int fPower, int fVerbose );

/*=== aig/gia/giaAig.c ================================================*/
extern Gia_Man_t * Gia_ManFromAig( Aig_Man_t * p );

extern void Bmatch_EncodeControlSignal(int row, int col, int nControlPi, int fCompl, AutoBuffer<int> &pLits);

// bmatchQbf.cpp
extern Bmatch_Qbf_Man_t *Bmatch_Gia_QbfAlloc(Gia_Man_t *pGia, int nPars, int fVerbose);
extern void Bmatch_Gia_QbfFree( Bmatch_Qbf_Man_t * p );
extern int Bmatch_Gia_QbfSolveValue(Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);
extern int Bmatch_Gia_QbfSolveValueInt(Bmatch_Qbf_Man_t *p, Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);

void Bmatch_InitQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_FillPossibleMIbyStrSupp(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
void Bmatch_ReducePossibleMIbyUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
Abc_Ntk_t *Bmatch_NtkQbfMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
InputMapping Bmatch_SolveQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

InputMapping Bmatch_SolveQbfInputSolver2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

#ifdef __cplusplus
}
#endif

// Abc_Obj_t *Bmatch_NtkCreateMultiplexer2(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal) {
//     int j;
//     std::vector<Abc_Obj_t *> prevSignal = pSignal;
//     std::vector<Abc_Obj_t *> currSignal;
    
//     for (int i = 0; i < pControl.size(); ++i) {
//         Abc_Obj_t *pC = pControl[i];
//         for (j = 0; j < prevSignal.size(); j += 2) {
//             Abc_Obj_t *p0 = prevSignal[j];
//             Abc_Obj_t *p1 = ((j + 1) == prevSignal.size()) ? prevSignal[j] : prevSignal[j + 1];
//             currSignal.push_back(Abc_AigMux(pMan, pC, p1, p0));
//         }
//         prevSignal = currSignal;
//         currSignal.clear();
//     }

//     return prevSignal[0];
// }

// void Bmatch_NtkBuildWithCone(Abc_Ntk_t *pNtk, Abc_Ntk_t *pNtkMiter, std::vector<int> &outs) {
//     Abc_Obj_t *pObj;
//     int i;

//     std::set<Abc_Obj_t *> cone;
//     for (auto &out : outs) {
//         Vec_Ptr_t * vNodes;
//         Abc_Obj_t *pOut = Abc_NtkPo(pNtk, out);
//         vNodes = Abc_NtkDfsNodes( pNtk, &pOut, 1 );
//         Vec_PtrForEachEntry(Abc_Obj_t *, vNodes, pObj, i) {
//             cone.insert(pObj);
//         }
//     }
// }

// Abc_Ntk_t *Bmatch_NtkQbfMiterReduced(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, std::vector<int> &forceYi2Xi, vMatch &MO) {
//     char Buffer[1000];
//     Abc_Ntk_t * pNtkMiter;
//     pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
//     Abc_Obj_t *pObj, *pObjTemp, *pObjNew;
//     auto &MapReduceMI = pMan->MapReduceMI;
//     int nControlPi = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1)))); // normal matrix
//     int i;
    
//     sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
//     Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));

//     Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
//     Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

//     // Control PI
//     std::vector<std::vector<Abc_Obj_t *> > controlPis;
//     for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
//         std::vector<Abc_Obj_t *> controlPi;
//         for (int j = 0; j < nControlPi; ++j) {
//             pObj = Abc_NtkCreatePi(pNtkMiter);
//             sprintf(Buffer, "controlPi_%d_%d", i, j);
//             Abc_ObjAssignName(pObj, Buffer, NULL);
//             controlPi.emplace_back(std::move(pObj));
//         }
//         pObj = Abc_NtkCreatePi(pNtkMiter);
//         sprintf(Buffer, "controlPi_%d_%d_inv", i, Abc_NtkPiNum(pNtk1));
//         Abc_ObjAssignName(pObj, Buffer, NULL);
//         controlPi.emplace_back(std::move(pObj));

//         controlPis.emplace_back(std::move(controlPi));
//     }
//     assert(controlPis.size() == Abc_NtkPiNum(pNtk2));

//     // Ntk1 PI
//     std::vector<Abc_Obj_t *> Ntk1_Pis;
//     Abc_NtkForEachPi(pNtk1, pObj, i) {
//         pObjNew = Abc_NtkCreatePi(pNtkMiter);
//         sprintf(Buffer, "_ntk1");
//         Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), Buffer);
//         Ntk1_Pis.emplace_back(pObjNew);
//         Ntk1_Pis.emplace_back(Abc_ObjNot(pObjNew));
//         pObj->pCopy = pObjNew;
//     }
//     Ntk1_Pis.emplace_back(Abc_AigConst1(pNtkMiter));
//     Ntk1_Pis.emplace_back(Abc_ObjNot(Abc_AigConst1(pNtkMiter)));

//     // PO
//     pObjNew = Abc_NtkCreatePo(pNtkMiter);
//     Abc_ObjAssignName(pObjNew, "miter", NULL);

//     // Create control input of Ntk2
//     // Init map of reduce MI
//     MapReduceMI = Mat(Abc_NtkPiNum(pNtk2), std::vector<int>());
//     Abc_NtkForEachPi(pNtk2, pObj, i) {
//         std::vector<Abc_Obj_t *> &pControl = controlPis[i];

//         if (forceYi2Xi[i] < 0) { // non forced
//             std::vector<Abc_Obj_t *> Reduced_Ntk1_Pis;
//             for (int j = 0; j < 2 * (Abc_NtkPiNum(pNtk1) + 1); ++j) {
//                 if (Bmatch_LegalMI(pMan, i, j)) {
//                     MapReduceMI[i].push_back(j);
//                     Reduced_Ntk1_Pis.push_back(Ntk1_Pis[j]);
//                 }
//             }
//             pObjNew = Bmatch_NtkCreateMultiplexer2((Abc_Aig_t *)pNtkMiter->pManFunc, pControl, Reduced_Ntk1_Pis);
//         } else { // force connection
//             pObjNew = Ntk1_Pis[forceYi2Xi[i]];
//         }
//         pObj->pCopy = pObjNew;
//     }

//     // Construct over Ntk1 and Ntk2
//     Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
//     Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);

//     std::vector<Abc_Obj_t *> Xnors;
//     Abc_NtkForEachCo(pNtk1, pObj, i) {
//         if (!MO[i].empty()) {
//             pObjTemp = Abc_ObjChild0Copy(pObj);
//             for (int j = 0; j < MO[i].size(); ++j) {
//                 pObj = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
//                 pObj = (MO[i][j].sign()) ? Abc_ObjNot(pObj) : pObj;
//                 pObjNew = Abc_ObjNot(Abc_AigXor((Abc_Aig_t *)pNtkMiter->pManFunc, pObjTemp, pObj));
//                 Xnors.push_back(pObjNew);
//             }
//         }
//     }
//     pObj = Bmatch_NtkCreateAnd((Abc_Aig_t *)pNtkMiter->pManFunc, Xnors);
//     Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pObj);

//     // Cleanup
//     Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

//     // DC2
//     if (Abc_NtkNodeNum(pNtkMiter) > 500) {
//         Abc_Ntk_t *pNtkTemp;
//         pNtkMiter = Abc_NtkDC2(pNtkTemp = pNtkMiter, 1, 0, 1, 0, 0);
//         Abc_NtkDelete(pNtkTemp);
//     }

//     if (!Abc_NtkCheck(pNtkMiter)) {
//         Abc_Print(-1, "Bmatch_NtkMiter: The network check has failed.\n");
//         Abc_NtkDelete(pNtkMiter);

//         return NULL;
//     }

//     // Io_Write(pNtkMiter, "miter_reduced.v", IO_FILE_VERILOG);

//     return pNtkMiter;
// }

// // Recursive solving
// InputMapping Bmatch_SolveQbfInputSolver2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
//     for (int fi = 0; fi < MO.size(); ++fi) {
//         for (auto &gi : MO[fi]) {

//         }
//     }
// }

inline int Bmatch_LegalMI(Bmatch_Man_t *pMan, int n, int m) {
    return pMan->possibleMI[n * pMan->mi + m];
}

void Bmatch_InitQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int nControlPi = pMan->nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1))); // nPi + const (not include inv)
    int mi = pMan->mi = (Abc_NtkPiNum(pNtk1) + 1) * 2;
    int ni = pMan->ni = Abc_NtkPiNum(pNtk2);
    auto &possibleMI = pMan->possibleMI;

    possibleMI.resize(ni * mi);
    possibleMI.fill(0);
}

Abc_Ntk_t *Bmatch_NtkQbfMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    Abc_Obj_t *pObj, *pObjTemp, *pObjNew;
    auto &MapReduceMI = pMan->MapReduceMI;
    int nControlPi = pMan->nControlPi;
    int i;
    
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    Abc_NtkSetName(pNtkMiter, Extra_UtilStrsav(Buffer));

    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // Control PI
    std::vector<std::vector<Abc_Obj_t *> > controlPis;
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPi;
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

    // Ntk1 PI
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
    Abc_ObjAssignName(pObjNew, "miter", NULL);

    // Create control input of Ntk2
    // Init map of reduce MI
    MapReduceMI = Mat(Abc_NtkPiNum(pNtk2), std::vector<int>());
    Abc_NtkForEachPi(pNtk2, pObj, i) {
        std::vector<Abc_Obj_t *> &pControl = controlPis[i];

        std::vector<Abc_Obj_t *> Reduced_Ntk1_Pis;
        for (int j = 0; j < Abc_NtkPiNum(pNtk1) + 1; ++j) {
            if (Bmatch_LegalMI(pMan, i, j * 2) || Bmatch_LegalMI(pMan, i, j * 2 + 1)) {
                MapReduceMI[i].push_back(j);
                Reduced_Ntk1_Pis.push_back(Ntk1_Pis[j]);
            }
        }

        pObjNew = Bmatch_NtkCreateMultiplexer((Abc_Aig_t *)pNtkMiter->pManFunc, pControl, Reduced_Ntk1_Pis, 0);
        pObj->pCopy = pObjNew;
    }

    // Construct over Ntk1 and Ntk2
    Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
    Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);

    std::vector<Abc_Obj_t *> Xnors;
    Abc_NtkForEachCo(pNtk1, pObj, i) {
        if (!MO[i].empty()) {
            pObjTemp = Abc_ObjChild0Copy(pObj);
            for (int j = 0; j < MO[i].size(); ++j) {
                pObj = Abc_ObjChild0Copy(Abc_NtkPo(pNtk2, MO[i][j].var()));
                pObj = (MO[i][j].sign()) ? Abc_ObjNot(pObj) : pObj;
                pObjNew = Abc_ObjNot(Abc_AigXor((Abc_Aig_t *)pNtkMiter->pManFunc, pObjTemp, pObj));
                Xnors.push_back(pObjNew);
            }
        }
    }
    pObj = Bmatch_NtkCreateAnd((Abc_Aig_t *)pNtkMiter->pManFunc, Xnors);
    Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, 0), pObj);

    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    // DC2
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

    // Io_Write(pNtkMiter, "miter_reduced.v", IO_FILE_VERILOG);

    return pNtkMiter;
}

void Bmatch_CnfReduce(std::vector<AutoBuffer<int> > &cnf) {
    for (auto iter = cnf.begin(); iter != cnf.end(); ) {
        auto &c1 = *iter;
        bool reduce = false;
        for (auto &c2 : cnf) {
            if (c1.size() != c2.size()) continue;

            int index = -1;
            for (int i = 0; i < c1.size(); ++i) {
                if (c1[i] != c2[i] && c1[i] == -c2[i]) {
                    if (true == reduce) {
                        reduce = false;
                        break;
                    } else {
                        index = i;
                        reduce = true;
                    }
                } else if (c1[i] != c2[i]) {
                    break;
                }
            }
            if (reduce) {
                for (int i = index; i < c2.size() - 1; ++i) {
                    c2[i] = c2[i + 1];
                }
                c2.resize(c2.size() - 1);
                break;
            }
        }
        if (reduce) {
            iter = cnf.erase(iter);
        } else {
            ++iter;
        }
    }
}

int Bmatch_PruneQbfSynthesizerByPossibleMI(CaDiCaL::Solver *pSatSyn, Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    auto &possibleMI = pMan->possibleMI;
    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1)));

    std::vector<AutoBuffer<int> > cnfs;
    for (int i = 0; i < pMan->ni; ++i) {
        for (int j = pMan->MapReduceMI[i].size() * 2; j < (1 << (nControlPi + 1)); ++j) {
        // for (int j = pMan->MapReduceMI[i].size() * 2; j < pMan->mi; ++j) {
            AutoBuffer<int> pLits(nControlPi + 1);
            Bmatch_EncodeControlSignal(i, j, nControlPi, 1, pLits);
            cnfs.emplace_back(pLits);
        }
    }

    Bmatch_CnfReduce(cnfs);

    for (auto &pLits : cnfs) {
        Bmatch_sat_solver_addclause(pSatSyn, pLits, pLits + pLits.size());
    }

    return 1;
}

InputMapping Bmatch_SolveQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    Abc_Ntk_t *pNtkMiter = Bmatch_NtkQbfMiter(pMan, pNtk1, pNtk2, MO);
    Aig_Man_t *pAig = Abc_NtkToDar(pNtkMiter, 0, 0);
    Gia_Man_t *pGia = Gia_ManFromAig(pAig);

    int nControlPi = (int)(std::ceil(std::log2(Abc_NtkPiNum(pNtk1) + 1))); // nPi + const (not include inv)
    int nPars = (nControlPi + 1) * Abc_NtkPiNum(pNtk2);

    Bmatch_Qbf_Man_t *pQbfMan = Bmatch_Gia_QbfAlloc(pGia, nPars, 0);
    CaDiCaL::Solver  *pSatSyn = pQbfMan->pSatSyn;

    Bmatch_PruneQbfSynthesizerByPossibleMI(pSatSyn, pMan, pNtk1, pNtk2, MO);
    // Bmatch_PruneSynthesizerByFunSupport(pSatSyn, pMan, pNtk1, pNtk2, MO);

    Vec_Int_t *vControl = Vec_IntAlloc(nPars);
    for (int i = 0; i < nPars; ++i)
        Vec_IntSetEntry(vControl, i, i);

    int result = Bmatch_Gia_QbfSolveValueInt(pQbfMan, pGia, vControl, nPars, 1024, 0, 100, 0, 0);
    Bmatch_Gia_QbfFree(pQbfMan);

    Gia_ManStop(pGia);
    Aig_ManStop(pAig);
    Abc_NtkDelete(pNtkMiter);

    int RetValue = 1;
    if (result == 1 || result == -1) {
        RetValue = 0;
    } else if (result == 0) {
        for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
            int decode = 0;
            // printf("ReduceMI: %d MI: %d\n", pMan->MapReduceMI[i].size(), pMan->mi);
            for (int j = nControlPi - 1; j >= 0; --j) {
                int value = (1 << j) * Vec_IntEntry(vControl, i * (nControlPi + 1) + j);
                decode = std::min(decode + value, (int)pMan->MapReduceMI[i].size());
            }
            decode = pMan->MapReduceMI[i][decode];
            decode = decode * 2 + Vec_IntEntry(vControl, i * (nControlPi + 1) + nControlPi);
            // printf("Decode: %d %d\n", decode, pMan->mi);
            for (int j = 0; j < pMan->mi; ++j) {
                if (j == decode) {
                    // printf("Resize: %d\n", resize);
                    MI[decode / 2].push_back(Literal(i, decode & 1));
                }
            }
        }
    }

    Vec_IntFree(vControl);

    return {RetValue, MI};
}

void Bmatch_FillPossibleMIbyStrSupp(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;
    auto &possibleMI = pMan->possibleMI;
    
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    possibleMI[n * mi + m * 2] = 1;
                    possibleMI[n * mi + m * 2 + 1] = 1;
                }
            }
        }
    }
    for (int yi = 0; yi < Abc_NtkPiNum(pNtk2); ++yi) {
        possibleMI[yi * mi + mi - 1] = 1;
        possibleMI[yi * mi + mi - 2] = 1;
    }
}

void Bmatch_ReducePossibleMIbyUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    auto &unateMat1 = pMan->unateMat1;
    auto &unateMat2 = pMan->unateMat2;
    auto &possibleMI = pMan->possibleMI;

    int mi = pMan->mi;
    int ni = pMan->ni;

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (int xi = 0; xi < mi / 2 - 1; ++xi) {
                int unateness1 = unateMat1[fi][xi];

                for (int yi = 0; yi < ni; ++yi) {
                    int unateness2 = unateMat2[gi][yi];

                    if (unateness1 == -1 || unateness2 == -1 || unateness1 == 3 || unateness2 == 3) continue;

                    if (unateness1 == unateness2) {
                        possibleMI[yi * mi + xi * 2 + (1 - g.sign())] = 0;
                    } else if (unateness1 != unateness2) {
                        possibleMI[yi * mi + xi * 2 + g.sign()] = 0;
                    }
                }
            }
        }
    }
}

ABC_NAMESPACE_IMPL_END