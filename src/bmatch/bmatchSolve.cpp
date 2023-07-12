#include "bmatch.hpp"
#include <vector>
#include <map>

// #include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option, char *output);
void Bmatch_WriteOutput(char *output, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
int Bmatch_CalculateScore(vMatch &MI, vMatch &MO);
void Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan);
std::vector<std::pair<int, int> > Bmatch_OneToOneCheck(Bmatch_Man_t *pMan);

//not used
void Bmatch_OutputLearnCase6(Bmatch_Man_t *pMan, int n, int m);

#ifdef __cplusplus
}
#endif

// x1
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}};

// x2
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {}, {}};
//#define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}, {}, {}, {}};
// #define OUTPUT_MAPPING vMatch MO_test = {{}, {Literal(1, false)}, {}, {}};
// case 0
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}, {}, {Literal(1, false), Literal(2, true)}};

// case 1
// #define OUTPUT_MAPPING vMatch MO = {{Literal()}};

// case 10
// #define OUTPUT_MAPPING vMatch MO = {{Literal(1, true)}, {Literal(0, true)}};

// case 14
// #define OUTPUT_MAPPING vMatch MO = {{Literal(5)}, {Literal(3)}, {Literal(6)}, {Literal(0)}, {Literal(2)}, {Literal(1)}, {Literal(4)}};

// case 15
// #define OUTPUT_MAPPING vMatch MO = {{Literal(8)}, {Literal(0)}, {Literal(2)}, {Literal(4)}, {Literal(5)}, {Literal(3)}, {Literal(1)}, {Literal(9)}, {Literal(6)}, {Literal(7)}};

// case 16
// #define OUTPUT_MAPPING vMatch MO = {{Literal(4, true)}, {Literal(6, true)}, {Literal(7, true)}, {Literal(1)}, {Literal(2)}, {Literal(5)}, {Literal(0)}, {Literal(3)}};

// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, true)}, {Literal(1, true)}, {Literal(2, true)}, {Literal(3, true)}};
// #define OUTPUT_MAPPING vMatch MO = {{Literal(3)}, {Literal(2)}, {Literal(1)}, {Literal(0)}};
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0)}, {Literal(1)}, {Literal(2)}, {Literal(3)}};

// OUTPUT_MAPPING

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option, char *output) {
    int maxIter = 100000, iter = 0, tried = 0, best = 0;
    EcResult result;
    vMatch MI, MO;
    int OutputSolverMode = 3;
    int inputSolverMode = 5;
    bool optimal = false;
    bool ret = true;

    // if (option & VERBOSE_MASK) Abc_NtkPrintIo(stdout, pNtk1, 0);
    // if (option & VERBOSE_MASK) Abc_NtkPrintIo(stdout, pNtk2, 0);
    if (option & VERBOSE_MASK) Bmatch_PrintBusInfo(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintInputSupport(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintOutputSupport(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintSymm(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintUnate(pMan, pNtk1, pNtk2);
    // if (option & VERBOSE_MASK) Bmatch_PrintEqual(pMan, pNtk1, pNtk2);



    //np3
    if (Abc_NtkPiNum(pNtk1) * Abc_NtkPiNum(pNtk2) < 50) {
        if (Abc_NtkPoNum(pNtk1) * Abc_NtkPoNum(pNtk2) < 50)
            inputSolverMode = 1;
    }
    Bmatch_SolveOutputGroup(pMan);
    // if (option & VERBOSE_MASK) Bmatch_PrintOutputGroup(pNtk1, pNtk2, pMan->Groups);

    int idealMax = Abc_NtkPoNum(pNtk1) + Abc_NtkPoNum(pNtk2);
    printf("\n-------------------------------------------------\n");
    printf("Optimal: %d\n", idealMax);
    abctime clkTotal = Abc_Clock();

    // Exhausting matching using QBF to solve all possible MI/MO in once
    if (Bmatch_ExhaustingMatching(pMan, pNtk1, pNtk2, MI, MO)) {
        result = Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO, 1, 0);
        printf("Exhausting matching verify: %s\n", result.status == 3 ? "PASS" : "FAIL");
        if (result.status == 3) {
            Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO);
            optimal = (best = Bmatch_CalculateScore(MI, MO)) == idealMax;
            Bmatch_WriteOutput(output, pNtk1, pNtk2, MI, MO);
        }
    } else {
        printf("Exhausting matching fail...\nDo Normal matching process...\n");
    }

    // Partition grouping
    Bmatch_PPCheck(pMan, pNtk1, pNtk2);

    //preprocess
    if (inputSolverMode == 1) Bmatch_InitControllableInputOutputMiter(pMan, pNtk1, pNtk2);
    Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);
    Bmatch_InitOutputSolver(pMan, pNtk1, pNtk2);
    ret &= Bmatch_PruneOutputSolverByUnate(pMan, pNtk1, pNtk2);

    if (option & VERBOSE_MASK) Bmatch_PrintOutputGroup(pNtk1, pNtk2, pMan->Groups);
    printf("Optimal: %d\n", 2*Abc_NtkPoNum(pNtk2));


    vMatch MO_try(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
    
    if(OutputSolverMode == 0 || OutputSolverMode == 1){
        MO = Bmatch_SolveOutput(pMan, pNtk1, pNtk2, NULL, NULL, 0);
        Bmatch_OutputLearn(pMan, false, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
    }
    // vMatch_Group MO;


    //for one to one partition case
    std::vector<std::pair<int, int> >OneToOneMap;
    int counter = 0;
    bool neg = true;
    int iterNeg = 0;
    if(OutputSolverMode == 3){
        OneToOneMap = Bmatch_OneToOneCheck(pMan);
        if(OneToOneMap.size() == 0){
            printf("not one to one map\n");
            OutputSolverMode = 0;
        }
        else printf("one to one map\n");
       
    }


    while (!optimal && ret) {
        //find new pair of output matching
        EcResult result;
        if(OutputSolverMode == 0 ){
            MO = Bmatch_SolveOutput(pMan, pNtk1, pNtk2, NULL, NULL, 0);
            // std::cout<<pMan->MO.size()<<" "<<pMan->MO.back().size()<<std::endl;
            // Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
        }
        else if(OutputSolverMode == 1){// assume there is optimal
            while(true){
                MO = Bmatch_SolveOutput(pMan, pNtk1, pNtk2, NULL, NULL, 0);
                std::cout<<pMan->MO.size()<<" "<<pMan->MO.back().size()<<std::endl;
                Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO);
                int score = 0;
                for (int i = 0; i < MO.size(); ++i) {
                    score += (int)!MO[i].empty() + MO[i].size();
                }
                if(score == 2*Abc_NtkPoNum(pNtk2)){
                    break;
                }
                // Match temp;
                if(MO.size() == 0) return;
                Bmatch_OutputLearn(pMan, true, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
            }
           
           
        }
        else if(OutputSolverMode == 3){
            MO_try[OneToOneMap[counter].first].push_back(Literal(OneToOneMap[counter].second, neg));
            // std::cout<<OneToOneMap[counter].first<<" "<<OneToOneMap[counter].second<<std::endl;
        }


        // MO_new = MO;
        //Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
        // Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_test);
        // if (MO_new.size() == 0) { break;} //all path traced

        
       
        
        //input solve
        if ((inputSolverMode == 5) && (OutputSolverMode == 0)) {
            auto Mapping = Bmatch_SolveQbfInputSolver3(pMan, pNtk1, pNtk2, MO);
            result.status = (Mapping.status == 0) ? NON_EQUIVALENT : EQUIVALENT;
            MI = Mapping.MI;
        } else if ((inputSolverMode == 5) && (OutputSolverMode == 3)) {
            auto Mapping = Bmatch_SolveQbfInputSolver3(pMan, pNtk1, pNtk2, MO_try);
            result.status = (Mapping.status == 0) ? NON_EQUIVALENT : EQUIVALENT;
            MI = Mapping.MI;
        } else if (inputSolverMode == 4) {
            Bmatch_InitQbfInputSolver(pMan, pNtk1, pNtk2);
            Bmatch_FillPossibleMIbyStrSupp(pMan, pNtk1, pNtk2, MO);
            Bmatch_ReducePossibleMIbyUnate(pMan, pNtk1, pNtk2, MO);
            auto Mapping = Bmatch_SolveQbfInputSolver(pMan, pNtk1, pNtk2, MO);
            result.status = (Mapping.status == 0) ? NON_EQUIVALENT : EQUIVALENT;
            MI = Mapping.MI;
        } else {
            Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);
            Bmatch_PruneInputSolverByStrSupport(pMan, MO);
            Bmatch_PruneInputSolverByFuncSupport(pMan, MO);
            Bmatch_PruneInputSolverBySymmetryProperty(pMan, MO);
            Bmatch_PruneInputSolverByUnate(pMan, MO);
            Bmatch_PruneInputSolverByBusOrdered(pMan, pNtk1, pNtk2);
            Bmatch_PruneInputSolverByBusExactMap(pMan, pNtk1, pNtk2);
            Bmatch_ApplyInputSolverRowConstraint(pMan, pNtk1, pNtk2);

            if (inputSolverMode == 2) Bmatch_InitControllableInputMiter(pMan, pNtk1, pNtk2, MO);

            if (inputSolverMode == 3) {
                auto Mapping = Bmatch_SolveInputQbf(pMan, pNtk1, pNtk2, MO);
                result.status = (Mapping.status == 0) ? NON_EQUIVALENT : EQUIVALENT;
                MI = Mapping.MI;
            } else {
                while (result.status != EQUIVALENT && iter++ < maxIter) {
                    Bmatch_PruneInputSolverByCounterPart(pMan, pNtk1, pNtk2, result.model, MI, MO);
                    Bmatch_PruneInputSolverBySymmetry(pMan, MI);
                    auto Mapping = Bmatch_SolveInput(pMan, pNtk1, pNtk2, NULL, NULL, 0);
                    if (Mapping.status == 0) break;
                    MI = Mapping.MI;

                    // result = Bmatch_NtkEcGia(pNtk1, pNtk2, MI, MO); // trash
                    result = (inputSolverMode == 1) ? Bmatch_NtkControllableInputOutputEcFraig(pMan, pNtk1, pNtk2, MI, MO)
                           : (inputSolverMode == 2) ? Bmatch_NtkControllableInputEcFraig(pMan, pNtk1, pNtk2, MI)
                           :                          Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO, 1, 0);
                }
            }
        }

        // Abc_PrintTime(1, "Current time", Abc_Clock() - clkTotal);

        if (result.status == EQUIVALENT) {
            printf("Find matching at iteration %d!!!\n", iter);
            if(OutputSolverMode == 0) {
                Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO);
                Bmatch_OutputLearn(pMan, true, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
            }
            else if(OutputSolverMode == 3) {
                Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_try);
                neg = true;
                counter ++;
                iterNeg = 0;
            }

            int score = 0;
            if(OutputSolverMode == 0){
               score = Bmatch_CalculateScore(MI, MO);
            } else if(OutputSolverMode == 3){
                score = counter*2;
            }
            optimal = score == 2*Abc_NtkPoNum(pNtk2);
            if (score > best) {
                Abc_PrintTime(1, "Current time", Abc_Clock() - clkTotal);
                printf("Optimal: %d Current: %d\n", 2*Abc_NtkPoNum(pNtk2), best = score);
                if(OutputSolverMode == 0) Bmatch_WriteOutput(output, pNtk1, pNtk2, MI, MO);
                else if(OutputSolverMode == 3) Bmatch_WriteOutput(output, pNtk1, pNtk2, MI, MO_try);
            }

        } else {
            if (option & VERBOSE_DETAIL_MASK) {
                if (iter - 1 == maxIter) printf("Reach maximum iteration (%d)!\n", maxIter);
                else printf("Input Solver UNSAT Mapping is infeasible using %d iterations\n", iter);
            }
            if(OutputSolverMode == 0) Bmatch_OutputLearn(pMan, false, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
            else if(OutputSolverMode == 3) {
                MO_try[OneToOneMap[counter].first].pop_back();
                neg = false;
                iterNeg++;
            }
            if(iterNeg > 1) {
                printf("partition map fail");
                OutputSolverMode = 0;
            }
            // Bmatch_OutputLearnCase6(pMan, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
            

        }
        iter = 0;
        ret = 1;
        tried++;
    }
    printf("Output Solver tried %d times\n", tried);
    printf("Best score: %d\n", best);
    Abc_PrintTime(1, "Total time", Abc_Clock() - clkTotal);
}

void Bmatch_WriteOutput(char *output, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO) {
    FILE *f = fopen(output, "w");
    if (!f) return;

    for (int xi = 0; xi < Abc_NtkPiNum(pNtk1); ++xi) {
        if (MI[xi].empty()) continue;
        fprintf(f, "INGROUP\n");
        fprintf(f, "1 + %s\n", Abc_ObjName(Abc_NtkPi(pNtk1, xi)));
        for (auto &yi : MI[xi]) {
            fprintf(f, "2 %c %s\n", (yi.sign()) ? '-' : '+', Abc_ObjName(Abc_NtkPi(pNtk2, yi.var())));
        }
        fprintf(f, "END\n");
    }

    for (int fi = 0; fi < Abc_NtkPoNum(pNtk1); ++fi) {
        if (MO[fi].empty()) continue;
        fprintf(f, "OUTGROUP\n");
        fprintf(f, "1 + %s\n", Abc_ObjName(Abc_NtkPo(pNtk1, fi)));
        for (auto &gi : MO[fi]) {
            fprintf(f, "2 %c %s\n", (gi.sign()) ? '-' : '+', Abc_ObjName(Abc_NtkPo(pNtk2, gi.var())));
        }
        fprintf(f, "END\n");
    }

    if (!MI.back().empty()) {
        fprintf(f, "CONSTGROUP\n");
        for (auto &yi : MI.back()) {
            fprintf(f, "%c %s\n", (yi.sign()) ? '+' : '-', Abc_ObjName(Abc_NtkPi(pNtk2, yi.var())));
        }
        fprintf(f, "END\n");
    }

    fclose(f);
}

std::vector<std::pair<int, int> > Bmatch_OneToOneCheck(Bmatch_Man_t *pMan){
    int c1, c2 = 0;
    std::vector<std::pair<int, int> > mapping;
    // std::cout<<c1<<" "<<c2<<std::endl;
    while((c1<pMan->oPartition1.size()) && (c2<pMan->oPartition2.size())){
        if((pMan->oPartition1[c1].size() > 1) || (pMan->oPartition2[c2].size() > 1)) {
            std::cout<<"par1 size:"<<pMan->oPartition1[c1].size()<<" par2 size:"<<pMan->oPartition2[c2].size()<<std::endl;
            mapping.clear();
            return mapping;
        }
        if(pMan->oPartition1[c1].size() == 0 ){
            c1++;
        }
        else if(pMan->oPartition2[c2].size() == 0 ){
            c2++;
        }
        else if((pMan->oPartition1[c1].size() == 1) && (pMan->oPartition2[c2].size() == 1)){
            // std::cout<<pMan->oPartition1[c1][0]<<" "<<pMan->oPartition2[c2][0]<<std::endl;
            mapping.emplace_back(std::make_pair(pMan->oPartition1[c1][0], pMan->oPartition2[c2][0]));
            c1++;
            c2++;
                        
        }
        
        // if(c1 != c2) break; 
    }
    return mapping;
}

void Bmatch_OutputLearnCase6(Bmatch_Man_t *pMan, int n, int m){
    vMatch_Group &MO = pMan->MO;
    AutoBuffer<int> pLits(n*m+1, 0);
    for(int i = 0; i<MO[0].size(); i++){
        pLits[i] = Bmatch_toLitCond(MO[0][i], 1);
    }
    Bmatch_sat_solver_addclause(pMan->pOutputSolver, pLits, pLits + MO[0].size());
    pMan->MO.pop_back();
}

int Bmatch_CalculateScore(vMatch &MI, vMatch &MO) {
    int score = 0;
    for (int i = 0; i < MO.size(); ++i) {
        score += (int)!MO[i].empty() + MO[i].size();
    }
    return score;
}


void Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan) {
    vGroup groups;
    auto &supp1 = pMan->vSuppInfo1;
    auto &supp2 = pMan->vSuppInfo2;
    int suppFunc1 = 0, suppFunc2 = 0;
    int n1 = supp1.size() - 1;
    int n2 = supp2.size() - 1;

    groups.emplace_back(std::vector<int>(), std::vector<int>());
    for (int i = n1, j = n2; i >= 0 || j >= 0; --i, --j) {
        auto &group = groups.back();
        if (i >= 0) group.first.emplace_back(std::get<PO>(supp1[i]));
        if (j >= 0) group.second.emplace_back(std::get<PO>(supp2[j]));

        // suppFunc1 = std::get<SUPPFUNC>(supp1[i]);
        suppFunc1 = std::max(suppFunc1, std::get<SUPPFUNC>(supp1[i]));
        suppFunc2 = (j - 1 >= 0) ? std::get<SUPPFUNC>(supp2[j - 1]) : std::get<SUPPFUNC>(supp2[0]);

        if (suppFunc1 > suppFunc2 && (i > 0 || j > 0)) {
            suppFunc1 = 0;
            groups.emplace_back(std::vector<int>(), std::vector<int>());
        }
    }

    pMan->Groups = groups;
    // return groups;
}

ABC_NAMESPACE_IMPL_END
