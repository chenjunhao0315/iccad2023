#include "bmatch.hpp"
#include <vector>
#include <algorithm>

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

int Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_ApplyInputSolverRowConstraint(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverBySymmetryProperty(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverByUnate(Bmatch_Man_t *pMan, vMatch &MO);
int Bmatch_PruneInputSolverBySymmetry(Bmatch_Man_t *pMan, vMatch &MI);
int Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI, vMatch& MO);
InputMapping Bmatch_HeuristicSolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
InputMapping Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);

vMatch Bmatch_SolveOutput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);
void Bmatch_InitOutputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
int Bmatch_PruneOutputSolverByUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan);
void Bmatch_OutputLearn(Bmatch_Man_t *pMan, bool status, int n, int m);
bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m);
void Bmatch_New_Or(Bmatch_Man_t *pMan, int n, int m);

#ifdef __cplusplus
}
#endif

// x1
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}};

// x2
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {}, {}};
// #define OUTPUT_MAPPING vMatch MO_test = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}};
// #define OUTPUT_MAPPING vMatch MO_test = {{}, {Literal(1, false)}, {}, {}};
// case 0
// #define OUTPUT_MAPPING vMatch MO = {{Literal(0, false)}, {}, {Literal(1, false), Literal(2, true)}};

// case 10
// #define OUTPUT_MAPPING vMatch MO = {{Literal(1, true)}, {Literal(0, true)}};

// case 14
#define OUTPUT_MAPPING vMatch MO = {{Literal(5)}, {Literal(3)}, {Literal(6)}, {Literal(0)}, {Literal(2)}, {Literal(1)}, {Literal(4)}};

// case 15
// #define OUTPUT_MAPPING vMatch MO = {{Literal(8)}, {Literal(0)}, {Literal(2)}, {Literal(4)}, {Literal(5)}, {Literal(3)}, {Literal(1)}, {Literal(9)}, {Literal(6)}, {Literal(7)}};

// case 16
// #define OUTPUT_MAPPING vMatch MO = {{Literal(4, true)}, {Literal(6, true)}, {Literal(7, true)}, {Literal(1)}, {Literal(2)}, {Literal(5)}, {Literal(0)}, {Literal(3)}};

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    int maxIter = 5000, iter = 0, tried = 0, best = 0;
    int ret = 1;
    EcResult result;

    Abc_NtkPrintIo(stdout, pNtk1, 0);
    Abc_NtkPrintIo(stdout, pNtk2, 0);
    if (option & VERBOSE_MASK) Bmatch_PrintBusInfo(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintInputSupport(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintOutputSupport(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintSymm(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintUnate(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintEqual(pMan, pNtk1, pNtk2);

    //preprocess
    Bmatch_InitControllableInputOutputMiter(pMan, pNtk1, pNtk2);
    Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);
    Bmatch_SolveOutputGroup(pMan);
    Bmatch_InitOutputSolver(pMan, pNtk1, pNtk2);
    ret &= Bmatch_PruneOutputSolverByUnate(pMan, pNtk1, pNtk2);

    if (option & VERBOSE_MASK) Bmatch_PrintOutputGroup(pNtk1, pNtk2, pMan->Groups);

    abctime clkTotal = Abc_Clock();

    vMatch MI, MO_new;
    // vMatch_Group MO;
    bool optimal = false;
    while (!optimal && ret) {
        //find new pair of output matching
        EcResult result;
        MO_new = Bmatch_SolveOutput(pMan, pNtk1, pNtk2, NULL, NULL, 0);

        // MO_new = MO_test;
        // Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
        // Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_test);
        if (MO_new.size() == 0) { break;}
        
        //input solve
        Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);
        ret &= Bmatch_PruneInputSolverByStrSupport(pMan, MO_new);
        ret &= Bmatch_PruneInputSolverByFuncSupport(pMan, MO_new);
        ret &= Bmatch_PruneInputSolverBySymmetryProperty(pMan, MO_new);
        ret &= Bmatch_PruneInputSolverByUnate(pMan, MO_new);
        ret &= Bmatch_ApplyInputSolverRowConstraint(pMan, pNtk1, pNtk2);

        // int controllableMiter = 1;
        
        // if (controllableMiter) Bmatch_InitControllableInputMiter(pMan, pNtk1, pNtk2, MO_new);

        while (ret && result.status != EQUIVALENT && iter++ < maxIter) {
            ret &= Bmatch_PruneInputSolverByCounterPart(pMan, pNtk1, pNtk2, result.model, MI, MO_new);
            ret &= Bmatch_PruneInputSolverBySymmetry(pMan, MI);
            if (!ret) break;
            auto Mapping = Bmatch_SolveInput(pMan, pNtk1, pNtk2, NULL, NULL, 0);
            if (Mapping.status == 0) break;
            MI = Mapping.MI;

            // result = (controllableMiter) ? Bmatch_NtkControllableInputEcFraig(pMan, pNtk1, pNtk2, MI) : Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO_new, 1, 0);
            result = Bmatch_NtkControllableInputOutputEcFraig(pMan, pNtk1, pNtk2, MI, MO_new);
        }

        Abc_PrintTime(1, "Current time", Abc_Clock() - clkTotal);

        if (result.status == EQUIVALENT) {
            printf("Find matching at iteration %d!!!\n", iter);
            Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
            Bmatch_OutputLearn(pMan, true, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));

            int score = 0;
            for (int i = 0; i < MO_new.size(); ++i) {
                score += (int)!MO_new[i].empty() + MO_new[i].size();
            }
            optimal = score == 2*Abc_NtkPoNum(pNtk2);
            if (score > best) {
                printf("Optimal: %d Current: %d\n", 2*Abc_NtkPoNum(pNtk2), best = score);
            }
        } else {
            if (option & VERBOSE_DETAIL_MASK) {
                if (iter - 1 == maxIter) printf("Reach maximum iteration (%d)!\n", maxIter);
                else printf("Input Solver UNSAT Mapping is infeasible using %d iterations\n", iter);
            }
            Bmatch_OutputLearn(pMan, false, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
        }
        iter = 0;
        ret = 1;
        tried++;
    }
    printf("Output Solver tried %d times\n", tried);
    printf("Best score: %d\n", best);
    Abc_PrintTime(1, "Total time", Abc_Clock() - clkTotal);
}

void Bmatch_OutputLearn(Bmatch_Man_t *pMan, bool status, int n, int m){
    std::vector<int> &LearnedLevel = pMan->LearnedLevel;
    auto *pSolver = pMan->pOutputSolver;   
    std::vector<int> &LearnedAssumption = pMan->LearnedAssumption;
    vMatch_Group &MO = pMan->MO;

    if(status){
        //add new assumption base on i level(c)
        auto &MoBack=MO.back();
        LearnedLevel.emplace_back(MoBack.size());
        for(auto &match:MoBack){
            // sat_solver_push(pSolver, Bmatch_toLit(match));
            LearnedAssumption.emplace_back(Bmatch_toLit(match));
        }

        //add clause to constraint at lease one pair
        Bmatch_New_Or(pMan, n, m);

    }
    else{
        //add new assumptin base on i level(~c)
        auto &MoBack=MO.back();
        if(LearnedLevel.size() == 0) LearnedLevel.emplace_back(MoBack.size());
        else LearnedLevel.back() += MoBack.size();

        for(auto &match:MoBack){
            LearnedAssumption.emplace_back(Bmatch_toLitCond(match, 1));
        }
        MO.pop_back();

    }
    pMan->pOutputSolver = pSolver;
    pMan->LearnedLevel = LearnedLevel;
    pMan->LearnedAssumption = LearnedAssumption;
}

bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m){
    std::vector<int> LearnedLevel = pMan->LearnedLevel;
    auto *pSolver = pMan->pOutputSolver;
    std::vector<int> LearnedAssumption = pMan->LearnedAssumption;
    vMatch_Group &MO = pMan->MO;
    std::cout<<"start backtrack\n";
    for(int i = 0; i < LearnedLevel.back(); i++){
        LearnedAssumption.pop_back();
    }
    
    LearnedLevel.pop_back();
    if (LearnedLevel.size() == 0) LearnedLevel.emplace_back(1);
    else LearnedLevel.back() += 1;
    if (MO.size() != 0){
        LearnedAssumption.emplace_back(Bmatch_toLitCond(MO.back()[0], 1));
        MO.pop_back();
        Bmatch_New_Or(pMan, n, m);
        
    }
    else{
        printf("all possible path traced\n");
        return false;
    }

    pMan->pOutputSolver = pSolver;
    pMan->LearnedLevel = LearnedLevel;
    pMan->LearnedAssumption = LearnedAssumption;
    return true;
}

void Bmatch_New_Or(Bmatch_Man_t *pMan, int n, int m){
    std::vector<int> LearnedAssumption = pMan->LearnedAssumption;
    auto *pSolver = pMan->pOutputSolver;
    
    AutoBuffer<int> pLits(n*m+1);
    int cont = 0;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            if(std::find(LearnedAssumption.begin(), LearnedAssumption.end(), Bmatch_toLit(i*m+j)) == LearnedAssumption.end() & \
                std::find(LearnedAssumption.begin(), LearnedAssumption.end(), Bmatch_toLitCond(i*m+j, 1)) == LearnedAssumption.end()){
                pLits[cont] = Bmatch_toLit(i*m+j);
                cont++;
            }
        }
    }
    
    Bmatch_sat_solver_addvar(pSolver);
    pLits[m*n-LearnedAssumption.size()] = Bmatch_toLit(Bmatch_sat_solver_nvars(pSolver));
    pMan->ClauseControl.emplace_back(Bmatch_sat_solver_nvars(pSolver));

    int LitSize = cont + 1;
    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + LitSize);
    // ABC_FREE(pLits);
}


vMatch Bmatch_SolveOutput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose){
    vMatch MO_new(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
    auto *pSolver = pMan->pOutputSolver;
    // pSolver->set("verbose", 3);
    std::vector<int> learnedAssumption = pMan->LearnedAssumption;
    std::vector<int> ClauseControl = pMan->ClauseControl;
    //pSolver->verbosity = 2;
    int n = Abc_NtkPoNum(pNtk2);
    int m = 2 * (Abc_NtkPoNum(pNtk1));
    int LitSize = learnedAssumption.size()+ClauseControl.size()+1;
    // int LitSize = ClauseControl.size()+learnedAssumption.size();
    vMatch_Group &MO = pMan->MO;

    //add assumption
    AutoBuffer<int> pLit(learnedAssumption.size()+ClauseControl.size()+pMan->LearnedLevel.size()+1);
    for(int i =0; i<learnedAssumption.size(); i++){
        pLit[i] = learnedAssumption[i];
        // std::cout<<"assump"<<learnedAssumption[i]<<" ";
    }
    // std::cout<<std::endl;
    //clause control
    for(int i = 0; i<ClauseControl.size(); i++){
        if (i == ClauseControl.size()-1) pLit[learnedAssumption.size()+i] = Bmatch_toLitCond(ClauseControl[i], 1);
        else pLit[i+learnedAssumption.size()] = Bmatch_toLit(ClauseControl[i]);
        // std::cout<<"clause"<<ClauseControl[i]<<" ";
    }
    // std::cout<<std::endl;
    
    // //projective
    if (pMan->AllowProjection) pLit[learnedAssumption.size()+ClauseControl.size()] = Bmatch_toLit(pMan->Projective);
    else pLit[learnedAssumption.size()+ClauseControl.size()] = Bmatch_toLitCond(pMan->Projective, 1);
    // std::cout<<pLit[learnedAssumption.size()+ClauseControl.size()]<<" "<<pMan->AllowProjection<<std::endl;
    
    if (fVerbose) std::cout<<"output solve start\n";
    int status = Bmatch_sat_solver_solve(pSolver, pLit, pLit+LitSize, 0, 0, 0, 0);
    while(status == 20){
        if (fVerbose) std::cout<<"output match failed"<<std::endl;
        
        if (status == 20 && !pMan->AllowProjection){
            if (fVerbose) std::cout<<"projection on"<<std::endl;
            //projection on
            pMan->AllowProjection = true;
            pLit[LitSize-1] = Bmatch_toLit(pMan->Projective);
            status = Bmatch_sat_solver_solve(pSolver, pLit, pLit+LitSize, 0, 0, 0, 0);
        }
        // std::cout<<status<<std::endl;

        if (status == 20){
            if (fVerbose) std::cout<<"projection off"<<std::endl;
            bool endloop = Bmatch_OutputBacktrack(pMan, n, m);
            learnedAssumption = pMan->LearnedAssumption;
            ClauseControl = pMan->ClauseControl;

            vMatch temp;
            if (!endloop) return temp;
            //projection off
            pMan->AllowProjection = false;
            
            for(int i =0; i<learnedAssumption.size(); i++){
                pLit[i] = learnedAssumption[i];
            }
            // std::cout<<"clause control:";
            for(int i = 0; i<ClauseControl.size(); i++){
                // std::cout<<ClauseControl[i]<<" ";
                if (i == ClauseControl.size()-1) pLit[learnedAssumption.size()+i] = Bmatch_toLitCond(ClauseControl[i], 1);
                else pLit[i+learnedAssumption.size()] = Bmatch_toLit(ClauseControl[i]);
            }
            // std::cout<<std::endl;

            LitSize = ClauseControl.size()+learnedAssumption.size()+1;
            pLit[LitSize-1] = Bmatch_toLitCond(pMan->Projective, 1);

            // std::cout<<"assumptions ";
            // for(int i = 0;i<LitSize; i++){
            //     std::cout<<pLit[i]<<" ";
            // }
            // std::cout<<std::endl;
            status = Bmatch_sat_solver_solve(pSolver, pLit, pLit+LitSize, 0, 0, 0, 0);
        }
    }

    if (fVerbose) std::cout<<"new match found"<<std::endl;
    std::vector<int> new_match;
    if (fVerbose) { printf("       "); for (int j = 0; j < m; ++j) { printf("%c%-3s", ((j & 1) != 0) ? '~' : ' ', Abc_ObjName(Abc_NtkPo(pNtk1, j / 2))); }; printf("\n");  }
    for (int i = 0; i < n; ++i) {
        if (fVerbose) printf("%3s: ", Abc_ObjName(Abc_NtkPo(pNtk2, i)));
        for (int j = 0; j < m; ++j) {
            if (fVerbose) printf(" %3d", Bmatch_sat_solver_var_value(pSolver, i * m + j));
            if (Bmatch_sat_solver_var_value(pSolver, i * m + j)) {
                MO_new[j / 2].emplace_back(i, (int)((j & 1) != 0));
                bool add = true;
                for(int k=0; k<MO.size(); k++){
                    if(std::find(MO[k].begin(), MO[k].end(), i*m+j) != MO[k].end()){
                        add = false;
                        break;
                    }
                }
                
                if (add) new_match.emplace_back(i*m+j);
            }
        }
        if (fVerbose) printf("\n");
    }
    if (fVerbose) printf("\n");
    MO.emplace_back(new_match);
    if(fVerbose){
        for(auto &level:MO){
            printf("MO");
            for(int j = 0;j<level.size(); j++){
                printf(" %3d", level[j]);
            }
            printf("\n");
        }
    }
    
    pMan->MO = MO;
    // ABC_FREE(pLit);
    return MO_new;
}

int Bmatch_PruneOutputSolverByUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int ret = 1;
    auto *pSolver = pMan->pOutputSolver;
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
                int Lit = Bmatch_toLitCond(i * (2 * Abc_NtkPoNum(pNtk2)) + j * 2, 1);
                // printf("(%d, %d) ", i, j * 2);
                Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                Lit = Bmatch_toLitCond(i * (2 * Abc_NtkPoNum(pNtk2)) + j * 2 + 1, 1);
                // printf("(%d, %d) ", i, j * 2 + 1);
                Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
            }
        }
    }

    return ret;
}

void Bmatch_InitOutputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2){
    auto *pSolver = pMan->pOutputSolver;
    if (pSolver) Bmatch_sat_solver_delete(pSolver);
    pSolver = Bmatch_sat_solver_new();
    //pSolver->verbosity = 2;

    int n = Abc_NtkPoNum(pNtk2);
    int m = 2 * (Abc_NtkPoNum(pNtk1));
    AutoBuffer<int> pLits(m*n);
    
    Bmatch_sat_solver_setnvars(pSolver, n * m + 1);
    //at least one pair
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            pLits[i*m+j] = Bmatch_toLit(i*m+j);
        }
    }
    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + m*n);

    //(c+d)<=1 ==> (~c+~d)
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            pLits[0] = Bmatch_toLitCond(i*m+j, 1);
            for(int k = j+1; k<m; k++){
                if (k != j){
                    pLits[1] = Bmatch_toLitCond(i*m+k, 1);
                    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
                }
            }
        }
    }

    //output group
    vGroup vGroup = pMan->Groups;

    for(int l = 0; l< vGroup.size(); l++){
        auto Group_ntk1 = vGroup[l].first;
        for(int k = l+1; k<vGroup.size(); k++){
            auto Group_ntk2 = vGroup[k].second;
            for(int i = 0;i<Group_ntk1.size(); i++){
                for(int j = 0;j<Group_ntk2.size(); j++){
                    //for output_ntk1_i and output_ntk2_j not in same group
                    //add clause ~cij ~dij
                    pLits[0] = Bmatch_toLitCond(Group_ntk1[i]*2 + m*Group_ntk2[j],1);
                    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
                    pLits[0] = Bmatch_toLitCond(Group_ntk1[i]*2 + m*Group_ntk2[j] + 1,1);
                    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
                }
            }
            
        }
    }
    
    // allow projection
    pMan->AllowProjection = false;
    pMan->Projective = n*m;
    pLits[0] = Bmatch_toLit(n*m);//allow projection
    for(int i = 0; i<n; i++){
        for(int j = i+1; j<n; j++){
            for(int k = 0; k<m/2; k++){
                pLits[1] = Bmatch_toLitCond(i*m+k*2, 1);
                pLits[2] = Bmatch_toLitCond(j*m+k*2+1, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits+3);
                // std::cout<<i*m+k*2<<" "<<j*m+k*2+1<<std::endl;

                pLits[1] = Bmatch_toLitCond(i*m+k*2+1, 1);
                pLits[2] = Bmatch_toLitCond(j*m+k*2, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits+3);

                pLits[1] = Bmatch_toLitCond(i*m+k*2, 1);
                pLits[2] = Bmatch_toLitCond(j*m+k*2, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits+3);

                pLits[1] = Bmatch_toLitCond(i*m+k*2+1, 1);
                pLits[2] = Bmatch_toLitCond(j*m+k*2+1, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits+3);
            }
        }
    }

    //equal

    // ABC_FREE(pLits);
    pMan->pOutputSolver = pSolver;

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
