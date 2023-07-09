#include "bmatch.hpp"
#include "bmatchQbf.hpp"

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

/*=== aig/gia/giaQbf.c ================================================*/
extern int Gia_QbfSolveValue( Gia_Man_t * pGia, Vec_Int_t * vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose );

// bmatchQbf.cpp
extern Bmatch_Qbf_Man_t *Bmatch_Gia_QbfAlloc(Gia_Man_t *pGia, int nPars, int fVerbose);
extern void Bmatch_Gia_QbfFree( Bmatch_Qbf_Man_t * p );
extern int Bmatch_Gia_QbfSolveValue(Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);
extern int Bmatch_Gia_QbfSolveValueInt(Bmatch_Qbf_Man_t *p, Gia_Man_t *pGia, Vec_Int_t *vValues, int nPars, int nIterLimit, int nConfLimit, int nTimeOut, int fGlucose, int fVerbose);

Abc_Obj_t *Bmatch_NtkCreateMultiplexer3(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal);
Abc_Obj_t *Bmatch_NtkCreateConditionalEqual(Abc_Aig_t *pMan, Abc_Obj_t *pControl, Abc_Obj_t *p1, Abc_Obj_t *p2);
void Bmatch_NtkCreateSortingCircuit(Abc_Ntk_t *pNtkMiter, std::vector<Abc_Obj_t *> &pN, std::vector<Abc_Obj_t *> &pOut);
Abc_Ntk_t *Bmatch_NtkQbfMiter3(Bmatch_Man_t *pManBmatch, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Mat &possibleMI, Mat &possibleMO);
Mat Bmatch_CalculatePossibleMI(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
Mat Bmatch_CalculatePossibleMO(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
InputMapping Bmatch_SolveQbfInputSolver3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

Mat Bmatch_CalculatePossibleMI2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
Mat Bmatch_CalculatePossibleMO2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
InputMapping Bmatch_SolveQbfInputSolver4(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

#ifdef __cplusplus
}
#endif

InputMapping Bmatch_SolveQbfInputSolver4(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    Mat possibleMI(Abc_NtkPiNum(pNtk2), std::vector<int>(2 * (Abc_NtkPiNum(pNtk1) + 1), 1));
    Mat possibleMO(Abc_NtkPoNum(pNtk2), std::vector<int>(2 * Abc_NtkPoNum(pNtk1), 0));
    auto &oStrSupp1 = pMan->oStrSupp1;
    auto &oStrSupp2 = pMan->oStrSupp2;

    struct SortSupp {
        int id;
        std::set<int> supp;
    };
    std::vector<SortSupp> supp1, supp2;
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i)
        supp1.push_back({i, oStrSupp1[i]});
    for (int i = 0; i < Abc_NtkPoNum(pNtk2); ++i)
        supp2.push_back({i, oStrSupp2[i]});
    std::sort(supp1.begin(), supp1.end(), [](SortSupp &a, SortSupp &b) { return a.supp.size() < b.supp.size(); });
    std::sort(supp2.begin(), supp2.end(), [](SortSupp &a, SortSupp &b) { return a.supp.size() < b.supp.size(); });

    int poNum1 = Abc_NtkPoNum(pNtk1), poNum2 = Abc_NtkPoNum(pNtk2);
    int cur1 = 0, cur2 = 0, lsMatch = 0;
    while (cur2 < poNum2) {
        int suppNum2 = supp2[cur2].supp.size();
        while (cur1 < poNum1 && supp1[cur1].supp.size() <= suppNum2)
            ++cur1;

        int gi = supp2[cur2].id;

        std::set<int> allsupp1;
        for (int i = lsMatch; i < cur1; ++i) {
            int fi = supp1[i].id;
            possibleMO[gi][2 * fi + 0] = 1;
            possibleMO[gi][2 * fi + 1] = 1;
            allsupp1.insert(oStrSupp1[fi].begin(), oStrSupp1[fi].end());
        }

        for (int xi = 0; xi < Abc_NtkPiNum(pNtk1); ++xi) {
            if (allsupp1.count(xi)) continue;
            for (auto yi : oStrSupp2[gi]) {
                possibleMI[yi][xi * 2 + 0] = 0;
                possibleMI[yi][xi * 2 + 1] = 0;
            }
        }

        if (cur1 <= cur2 + 1) 
            lsMatch = cur1;
        ++cur2;
    }

    printf("possibleMI\n");
    for (int i = 0; i < possibleMI.size(); ++i) {
        for (int j = 0; j < possibleMI[i].size(); ++j) {
            printf("%d ", possibleMI[i][j]);
        }
        printf("\n");
    }
    printf("possibleMO\n");
    for (int i = 0; i < possibleMO.size(); ++i) {
        for (int j = 0; j < possibleMO[i].size(); ++j) {
            printf("%d ", possibleMO[i][j]);
        }
        printf("\n");
    }

    int limit_1 = 0;
    int limit_2 = 0;

    int try_freeze = 0, failed_freeze = 0;
    int best_score = 0;
    // int ideal_max_score = Abc_NtkPoNum(pNtk1) + Abc_NtkPoNum(pNtk2);
    int ideal_max_score = 2 * Abc_NtkPoNum(pNtk1);
    printf("ideal maximal score: %d\n", ideal_max_score);

    Mat possibleMIOld = possibleMI;
    Mat possibleMOOld = possibleMO;

    while (best_score < ideal_max_score) {
        while (limit_1 < Abc_NtkPoNum(pNtk1) - 1) {
            if (supp1[limit_1].supp.size() == supp1[limit_1 + 1].supp.size()) {
                ++limit_1;
            } else {
                break;
            }
        }
        while (limit_2 < Abc_NtkPoNum(pNtk2) - 1) {
            if (supp2[limit_2].supp.size() == supp2[limit_2 + 1].supp.size()) {
                ++limit_2;
            } else {
                break;
            }
        }

        for (int i = limit_1 + 1; i < Abc_NtkPoNum(pNtk1); ++i) {
            int fi = supp1[i].id;
            for (int gi = 0; gi < Abc_NtkPoNum(pNtk2); gi++) {
                possibleMO[gi][2 * fi + 0] = 0;
                possibleMO[gi][2 * fi + 1] = 0;
            }
        }
        for (int i = limit_2 + 1; i < Abc_NtkPoNum(pNtk2); ++i) {
            int gi = supp2[i].id;
            for (int fi = 0; fi < Abc_NtkPoNum(pNtk1); fi++) {
                possibleMO[gi][2 * fi + 0] = 0;
                possibleMO[gi][2 * fi + 1] = 0;
            }
        }

        int all_zero = 1;
        for (int i = 0; i < 2 * Abc_NtkPoNum(pNtk1); i++) {
            for (int j = 0; j < Abc_NtkPoNum(pNtk2); j++) {
                if (possibleMO[j][i] != 0) {
                    all_zero = 0;
                    break;
                }
            }
        }

        int should_reset = 1, should_go = 1;

        if (!all_zero) {
            int num_matchPo_1 = 0, num_matchPo_2 = 0;
            int nPo_1 = limit_1 + 1;
            int nPo_2 = limit_2 + 1;
            int best_nGroup = 0;
            Vec_Int_t *vPiValues = Vec_IntAlloc(0);

            Abc_Ntk_t *pNtkTemp;
            Abc_Ntk_t *pNtkMiter = Bmatch_NtkQbfMiter3(pMan, pNtk1, pNtk2, possibleMI, possibleMO);
            int nControlPi = Abc_NtkPiNum(pNtk2) * (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1))));
            int nControlPo = Abc_NtkPoNum(pNtk1) * Abc_NtkPoNum(pNtk2) * 2;
            int nPars = nControlPi + nControlPo;

            int nControlSortMain = Abc_NtkPoNum(pNtk1);
            int nControlSortSub  = Abc_NtkPoNum(pNtk2);

            int assert_index = nControlSortMain - nPo_1;

            while (assert_index + best_nGroup < nControlSortMain) {
                Abc_Obj_t *pOut = Abc_NtkPo(pNtkMiter, assert_index);
                Abc_Ntk_t *pNtkCone = Abc_NtkCreateCone(pNtkMiter, Abc_ObjFanin0(pOut), Abc_ObjName(pOut), 1);
                if (Abc_ObjFaninC0(pOut))
                    Abc_NtkPo(pNtkCone, 0)->fCompl0  ^= 1;
                pNtkCone = Abc_NtkDC2( pNtkTemp = pNtkCone, 0, 0, 1, 0, 0 );
                Abc_NtkDelete(pNtkTemp);

                Aig_Man_t * pAig = Abc_NtkToDar( pNtkCone, 0, 0 );
                Gia_Man_t * pGia = Gia_ManFromAig( pAig );
                Bmatch_Qbf_Man_t *pQbfMan = Bmatch_Gia_QbfAlloc(pGia, nPars, 0);

                abctime clkStart = Abc_Clock();
                int result = Bmatch_Gia_QbfSolveValueInt(pQbfMan, pGia, vPiValues, nPars, 1024, 0, 100, 0, 1);
                // int result = Gia_QbfSolveValue(pGia, vPiValues, nPars, 1024, 0, 100, 0, 1);
                ABC_PRT( "Time:", Abc_Clock() - clkStart );

                assert(nPars == Vec_IntSize(vPiValues));
                Bmatch_Gia_QbfFree(pQbfMan);
                Gia_ManStop( pGia );
                Aig_ManStop( pAig );
                Abc_NtkDelete(pNtkCone);
                if (result == 1) {
                    assert_index++;
                    continue;
                } else {
                    break;
                }
            }
            num_matchPo_1 = num_matchPo_2 = nControlSortMain - assert_index;
            if (best_nGroup < num_matchPo_1)
                best_nGroup = num_matchPo_1;
            best_score = num_matchPo_1 + num_matchPo_2;
            printf("best_score: %d\n", best_score);

            int nControl = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1))));
            for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
                int decode = 0;
                for (int j = nControl - 1; j >= 0; --j) {
                    int value = (1 << j) * Vec_IntEntry(vPiValues, i * nControl + j);
                    decode += value;
                    if (decode >= pMan->MapReduceMI[i].size())
                        decode -= value;
                }
                decode = pMan->MapReduceMI[i][decode];
                MI[decode / 2].push_back(Literal(i, decode & 1));
            }

            MO.clear();
            MO = vMatch(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
            for (int i = 0; i < Abc_NtkPoNum(pNtk1); i++) {
                for (int j = 0; j < 2 * Abc_NtkPoNum(pNtk2); j += 2) {
                    if (Vec_IntEntry(vPiValues, nControlPi + i * Abc_NtkPoNum(pNtk2) * 2 + j) == 1) {
                        MO[i].push_back(Literal(j / 2, false));
                    } else if (Vec_IntEntry(vPiValues, nControlPi + i * Abc_NtkPoNum(pNtk2) * 2 + j + 1) == 1) {
                        MO[i].push_back(Literal(j / 2, true));
                    }
                }
            }
        } else if (try_freeze) {
            failed_freeze = 1;
            should_reset = 1;
            should_go = 0;
            try_freeze = 0;
        }

        printf("%d %d\n", should_reset, should_go);

        if (should_go) {
            limit_1++;
            limit_2++;
            if (limit_1>=Abc_NtkPoNum(pNtk1) && limit_2>=Abc_NtkPoNum(pNtk2)) break;
            if (limit_1 >= Abc_NtkPoNum(pNtk1)) limit_1 = Abc_NtkPoNum(pNtk1)-1;
            if (limit_2 >= Abc_NtkPoNum(pNtk2)) limit_2 = Abc_NtkPoNum(pNtk2)-1;
        }

        possibleMO = possibleMIOld;
        // reset matrix if we didn't freeze it...
        if (should_reset) {
            possibleMI = possibleMIOld;
        }
    }

    return {1, MI};
}

InputMapping Bmatch_SolveQbfInputSolver3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    Mat possibleMO = Bmatch_CalculatePossibleMO(pMan, pNtk1, pNtk2, MO);
    Mat possibleMI = Bmatch_CalculatePossibleMI(pMan, pNtk1, pNtk2, MO);

    Abc_Ntk_t *pNtkTemp;
    Abc_Ntk_t *pNtkMiter = Bmatch_NtkQbfMiter3(pMan, pNtk1, pNtk2, possibleMI, possibleMO);

    int nControlPI = Abc_NtkPiNum(pNtk2) * (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1))));
    int nControlPO = Abc_NtkPoNum(pNtk1) * Abc_NtkPoNum(pNtk2) * 2;
    int nPars = nControlPI + nControlPO;
    Vec_Int_t *vPiValues = Vec_IntAlloc(0);

    int target = 0;
    for (auto &fi : MO)
        target += fi.size();

    int index = Abc_NtkPoNum(pNtk1) - target, result = -1;

    while (index < Abc_NtkPoNum(pNtk1)) {
        Abc_Obj_t *pOut = Abc_NtkPo(pNtkMiter, index);
        Abc_Ntk_t *pNtkCone = Abc_NtkCreateCone(pNtkMiter, Abc_ObjFanin0(pOut), Abc_ObjName(pOut), 1);
        if (Abc_ObjFaninC0(pOut))
            Abc_NtkPo(pNtkCone, 0)->fCompl0  ^= 1;
        pNtkCone = Abc_NtkDC2( pNtkTemp = pNtkCone, 0, 0, 1, 0, 0 );
        Abc_NtkDelete(pNtkTemp);

        Aig_Man_t * pAig = Abc_NtkToDar( pNtkCone, 0, 0 );
        Gia_Man_t * pGia = Gia_ManFromAig( pAig );
        Bmatch_Qbf_Man_t *pQbfMan = Bmatch_Gia_QbfAlloc(pGia, nPars, 0);

        abctime clkStart = Abc_Clock();
        result = Bmatch_Gia_QbfSolveValueInt(pQbfMan, pGia, vPiValues, nPars, 1024, 0, 100, 0, 1);
        // int result = Gia_QbfSolveValue(pGia, vPiValues, nPars, 1024, 0, 100, 0, 1);
        ABC_PRT( "Time:", Abc_Clock() - clkStart );

        assert(nPars == Vec_IntSize(vPiValues));
        Bmatch_Gia_QbfFree(pQbfMan);
        Gia_ManStop( pGia );
        Aig_ManStop( pAig );
        Abc_NtkDelete(pNtkCone);
        if (result == 1) {
            index++;
            continue;
        } else {
            break;
        }
    }

    if (Abc_NtkPoNum(pNtk1) - index != target) {
        printf("Target: %d MaxScore: %d\n", target,Abc_NtkPoNum(pNtk1) - index);
        return {0, MI};
    }

    int nControl = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1))));
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        int decode = 0;
        for (int j = nControl - 1; j >= 0; --j) {
            int value = (1 << j) * Vec_IntEntry(vPiValues, i * nControl + j);
            decode += value;
            if (decode >= pMan->MapReduceMI[i].size())
                decode -= value;
        }
        decode = pMan->MapReduceMI[i][decode];
        MI[decode / 2].push_back(Literal(i, decode & 1));
    }

    MO.clear();
    MO = vMatch(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); i++) {
        for (int j = 0; j < 2 * Abc_NtkPoNum(pNtk2); j += 2) {
            if (Vec_IntEntry(vPiValues, nControlPI + i*Abc_NtkPoNum(pNtk2)*2 + j) == 1) {
                MO[i].push_back(Literal(j / 2, false));
            } else if (Vec_IntEntry(vPiValues, nControlPI + i*Abc_NtkPoNum(pNtk2)*2 + j + 1) == 1) {
                MO[i].push_back(Literal(j / 2, true));
            }
        }
    }

    return {1, MI};
}

Abc_Ntk_t *Bmatch_NtkQbfMiter3(Bmatch_Man_t *pManBmatch, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, Mat &possibleMI, Mat &possibleMO) {
    char Buffer[1000];
    Abc_Ntk_t *pNtkMiter = Abc_NtkAlloc( ABC_NTK_STRASH, ABC_FUNC_AIG, 1 );
    sprintf(Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName);
    pNtkMiter->pName = Extra_UtilStrsav(Buffer);
    Abc_Aig_t *pMan = (Abc_Aig_t *)pNtkMiter->pManFunc;
    auto &MapReduceMI = pManBmatch->MapReduceMI;
    MapReduceMI = Mat(Abc_NtkPiNum(pNtk2), std::vector<int>());

    int i;
    int nControlPI = (int)(std::ceil(std::log2(2 * (Abc_NtkPiNum(pNtk1) + 1))));
    Abc_Obj_t *pObj, *pObjNew;

    // Const assign
    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // control PI
    std::vector<std::vector<Abc_Obj_t *> > controlPIs;
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        std::vector<Abc_Obj_t *> controlPI;
        for (int j = 0; j < nControlPI; ++j) {
            pObjNew = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "control_pi_%d_%d", i, j);
            Abc_ObjAssignName(pObjNew, Buffer, NULL);
            controlPI.push_back(pObjNew);
        }
        controlPIs.push_back(controlPI);
    }
    assert(Abc_NtkPiNum(pNtk2) == controlPIs.size());

    // control PO
    std::vector<std::vector<Abc_Obj_t *> > controlPOs;
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i) {
        std::vector<Abc_Obj_t *> controlPO;
        for (int j = 0; j < Abc_NtkPoNum(pNtk2); ++j) {
            pObjNew = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "control_po_%d_%d", i, j);
            Abc_ObjAssignName(pObjNew, Buffer, NULL);
            controlPO.push_back(pObjNew);
            pObjNew = Abc_NtkCreatePi(pNtkMiter);
            sprintf(Buffer, "control_po_%d_%d_inv", i, j);
            Abc_ObjAssignName(pObjNew, Buffer, NULL);
            controlPO.push_back(pObjNew);
        }
        controlPOs.push_back(controlPO);
    }

    // Ntk1 Input
    std::vector<Abc_Obj_t *> Ntk1Pis;
    Abc_NtkForEachPi(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), "_ntk1");
        pObj->pCopy = pObjNew;
        Ntk1Pis.push_back(pObjNew);
        Ntk1Pis.push_back(Abc_ObjNot(pObjNew));
    }
    Ntk1Pis.push_back(Abc_AigConst1(pNtkMiter));
    Ntk1Pis.push_back(Abc_ObjNot(Abc_AigConst1(pNtkMiter)));

    // Ntk2 Input
    Abc_NtkForEachPi(pNtk2, pObj, i) {
        auto &controlPI = controlPIs[i];
        std::vector<Abc_Obj_t *> reducePI;
        for (int j = 0; j < 2 * (Abc_NtkPiNum(pNtk1) + 1); ++j) {
            if (possibleMI[i][j]) {
                MapReduceMI[i].push_back(j);
                reducePI.push_back(Ntk1Pis[j]);
            }
        }
        pObj->pCopy = Bmatch_NtkCreateMultiplexer3(pMan, controlPI, reducePI);;
    }

    // Ntk1 Output
    Abc_NtkForEachPo(pNtk1, pObj, i) {
        pObjNew = Abc_NtkCreatePo(pNtkMiter);
        sprintf(Buffer, "cir1_%d", i);
        Abc_ObjAssignName(pObjNew, "miter_", Buffer);
    }

    // Ntk2 Output
    Abc_NtkForEachPo(pNtk2, pObj, i) {
        pObjNew = Abc_NtkCreatePo(pNtkMiter);
        sprintf(Buffer, "cir2_%d", i);
        Abc_ObjAssignName(pObjNew, "miter_", Buffer);
    }
    assert(Abc_NtkPoNum(pNtkMiter) == Abc_NtkPoNum(pNtk1) + Abc_NtkPoNum(pNtk2));

    // construct the network
    Bmatch_NtkMiterAddOne(pNtk1, pNtkMiter);
    Bmatch_NtkMiterAddOne(pNtk2, pNtkMiter);

    // equivalence checking
    std::vector<Abc_Obj_t *> ANDs(1, Abc_AigConst1(pNtkMiter));
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i)  {
        for (int j = 0; j < 2 * Abc_NtkPoNum(pNtk2); ++j) {
            Abc_Obj_t * Ntk1_PO = Abc_ObjChild0Copy( Abc_NtkPo(pNtk1, i)   );
            Abc_Obj_t * Ntk2_PO = Abc_ObjChild0Copy( Abc_NtkPo(pNtk2, j / 2) );
            Abc_Obj_t * control = controlPOs[i][j];
            
            if (possibleMO[i][j] == 0) {
                pObjNew = Abc_AigXor( pMan, control, Abc_AigConst1(pNtkMiter));
            } else {
                pObjNew = Bmatch_NtkCreateConditionalEqual(pMan, control, Ntk1_PO, (j & 1) ? Abc_ObjNot(Ntk2_PO) : Ntk2_PO);
            }
            ANDs.push_back(pObjNew);
        }
    }
    Abc_Obj_t *pUni = Bmatch_NtkCreateAnd(pMan, ANDs);

    // ensure at least one output matching be used
    std::vector<Abc_Obj_t *> POs_1;
    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i) {
        std::vector<Abc_Obj_t *> ORs;
        for (int j = 0; j < 2 * Abc_NtkPoNum(pNtk2); ++j) {
            if (possibleMO[i][j] == 0) continue;

            ORs.push_back(controlPOs[i][j]);
        }
        Abc_Obj_t *pEnsure = (ORs.empty()) ? Abc_ObjNot(Abc_AigConst1(pNtkMiter)) : Bmatch_NtkCreateOr(pMan, ORs);
        POs_1.push_back(Abc_AigAnd(pMan, pEnsure, pUni));
    }
    std::vector<Abc_Obj_t *> POs_2;
    for (int j = 0; j < 2 * Abc_NtkPoNum(pNtk2); j += 2) {
        std::vector<Abc_Obj_t *> ORs;
        for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i) {
            if (possibleMI[i][j])
                ORs.push_back(controlPOs[i][j + 0]);
            if (possibleMI[i][j + 1])
                ORs.push_back(controlPOs[i][j + 1]);
        }
        Abc_Obj_t *pEnsure = (ORs.empty()) ? Abc_ObjNot(Abc_AigConst1(pNtkMiter)) : Bmatch_NtkCreateOr(pMan, ORs);
        POs_2.push_back(Abc_AigAnd(pMan, pEnsure, pUni));
    }
    std::vector<Abc_Obj_t *> pOut_1, pOut_2;
    Bmatch_NtkCreateSortingCircuit(pNtkMiter, POs_1, pOut_1);
    Bmatch_NtkCreateSortingCircuit(pNtkMiter, POs_2, pOut_2);
    assert(pOut_1.size() == Abc_NtkPoNum(pNtk1));
    assert(pOut_2.size() == Abc_NtkPoNum(pNtk2));

    for (int i = 0; i < pOut_1.size(); ++i) {
        Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, i), pOut_1[i]);
    }
    for (int i = 0; i < pOut_2.size(); ++i) {
        Abc_ObjAddFanin(Abc_NtkPo(pNtkMiter, pOut_1.size() + i), pOut_2[i]);
    }

    // Cleanup
    Abc_AigCleanup((Abc_Aig_t*)pNtkMiter->pManFunc);

    if (!Abc_NtkCheck(pNtkMiter))
        Abc_Print(1, "Network check has failed.\n");

    // Io_Write(pNtkMiter, "miter.v", IO_FILE_VERILOG);

    return pNtkMiter;
}

void Bmatch_FillPossibleMIByStrSupp(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, Mat &possibleMI) {
    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    possibleMI[n][m * 2] = 1;
                    possibleMI[n][m * 2 + 1] = 1;
                }
            }
        }
    }
    for (int n = 0; n < Abc_NtkPiNum(pNtk2); ++n) {
        possibleMI[n][mi - 1] = 1;
        possibleMI[n][mi - 2] = 1;
    }
}

void Bmatch_ReducePossibleMIbyUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, Mat &possibleMI) {
    auto &unateMat1 = pMan->unateMat1;
    auto &unateMat2 = pMan->unateMat2;

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
                        possibleMI[yi][xi * 2 + (1 - g.sign())] = 0;
                    } else if (unateness1 != unateness2) {
                        possibleMI[yi][xi * 2 + g.sign()] = 0;
                    }
                }
            }
        }
    }
}

Mat Bmatch_CalculatePossibleMI(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    Mat possibleMI(Abc_NtkPiNum(pNtk2), std::vector<int>(2 * (Abc_NtkPiNum(pNtk1) + 1), 0));

    Bmatch_FillPossibleMIByStrSupp(pMan, pNtk1, pNtk2, MO, possibleMI);
    Bmatch_ReducePossibleMIbyUnate(pMan, pNtk1, pNtk2, MO, possibleMI);

    return possibleMI;
}

Mat Bmatch_CalculatePossibleMO(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO) {
    Mat possibleMO(Abc_NtkPoNum(pNtk2), std::vector<int>(2 * (Abc_NtkPoNum(pNtk1)), 0));

    for (int gi = 0; gi < Abc_NtkPoNum(pNtk2); ++gi) {
        for (auto fi : MO[gi]) {
            possibleMO[gi][fi.var() * 2 + fi.sign()] = 1;
            // possibleMO[gi][fi.var() * 2 + 1] = 1;
        }
    }

    return possibleMO;
}

Mat Bmatch_CalculatePossibleMI2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    Mat possibleMI(Abc_NtkPiNum(pNtk2), std::vector<int>(2 * (Abc_NtkPiNum(pNtk1) + 1), 1));

    return possibleMI;
}

Mat Bmatch_CalculatePossibleMO2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    Mat possibleMO(Abc_NtkPoNum(pNtk2), std::vector<int>(2 * (Abc_NtkPoNum(pNtk1)), 0));

    // support
    auto &group = pMan->Groups;
    for (auto &g : group) {
        for (int i = 0; i < g.first.size(); ++i) {
            for (int j = 0; j < g.second.size(); ++j) {
                possibleMO[g.second[j]][g.first[i] * 2 + 0] = 1;
                possibleMO[g.second[j]][g.first[i] * 2 + 1] = 1;
            }
        }
    }

    // unateness
    Mat &unateMat1 = pMan->unateMat1;
    Mat &unateMat2 = pMan->unateMat2;
    AutoBuffer<int> binate1(Abc_NtkPoNum(pNtk1));
    AutoBuffer<int> unate1(Abc_NtkPoNum(pNtk1));
    AutoBuffer<int> binate2(Abc_NtkPoNum(pNtk2));
    AutoBuffer<int> unate2(Abc_NtkPoNum(pNtk2));

    Bmatch_GetUnateCount(unateMat1, binate1, unate1);
    Bmatch_GetUnateCount(unateMat2, binate2, unate2);

    for (int i = 0; i < Abc_NtkPoNum(pNtk1); ++i) {
        int nSupp1 = binate1[i] + unate1[i];
        int nEquivUnate1 = nSupp1 + binate1[i];
        for (int j = 0; j < Abc_NtkPoNum(pNtk2); ++j) {
            int nSupp2 = binate2[j] + unate2[j];
            int nEquivUnate2 = nSupp2 + binate2[j];
            if (nSupp2 < nSupp1 || nEquivUnate2 < nEquivUnate1) {
                possibleMO[j][i * 2 + 0] = 0;
                possibleMO[j][i * 2 + 1] = 0;
            }
        }
    }

    return possibleMO;
}

void Bmatch_NtkCreateSortingCircuit(Abc_Ntk_t *pNtkMiter, std::vector<Abc_Obj_t *> &pN, std::vector<Abc_Obj_t *> &pOut) {
    pOut.clear();

    if (pN.size() == 0) {
        pOut.push_back(Abc_AigConst1(pNtkMiter));
        return;
    }
        
    if (pN.size() == 1) {
        pOut.push_back(pN[0]);
        return;
    }
    
    for (int i = 0; i < pN.size(); i++) {
        pOut.push_back(pN[i]);
    }
    
    for (int i = pOut.size() - 1; i >= 1; --i) {
        for (int j = 0; j < i; ++j) {
            Abc_Obj_t * pObjLarge = Abc_AigAnd( (Abc_Aig_t *)pNtkMiter->pManFunc, pOut[j], pOut[j + 1] );
            Abc_Obj_t * pObjSmall = Abc_AigOr ( (Abc_Aig_t *)pNtkMiter->pManFunc, pOut[j], pOut[j + 1] );
        
            pOut[j + 0] = pObjLarge;
            pOut[j + 1] = pObjSmall;
        }
    }
    return;
}

Abc_Obj_t *Bmatch_NtkCreateConditionalEqual(Abc_Aig_t *pMan, Abc_Obj_t *pControl, Abc_Obj_t *p1, Abc_Obj_t *p2) {
    Abc_Obj_t *pObjNew = Abc_ObjNot(Abc_AigXor(pMan, p1, p2));
    pObjNew = Abc_AigOr(pMan, Abc_ObjNot(pControl), pObjNew);

    return pObjNew;
}

Abc_Obj_t *Bmatch_NtkCreateMultiplexer3(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal) {
    std::vector<Abc_Obj_t *> prevSignal = pSignal;
    std::vector<Abc_Obj_t *> currSignal;
    
    for (int i = 0; i < pControl.size(); ++i) {
        Abc_Obj_t *pC = pControl[i];
        for (int j = 0; j < prevSignal.size(); j += 2) {
            Abc_Obj_t *p0 = prevSignal[j];
            Abc_Obj_t *p1 = ((j + 1) == prevSignal.size()) ? prevSignal[j] : prevSignal[j + 1];
            currSignal.push_back(Abc_AigMux(pMan, pC, p1, p0));
        }
        prevSignal = currSignal;
        currSignal.clear();
    }

    return prevSignal[0];
}

ABC_NAMESPACE_IMPL_END