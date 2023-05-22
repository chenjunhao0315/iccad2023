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
bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m, int verbose);
void Bmatch_New_Or(Bmatch_Man_t *pMan, int n, int m, int verbose);

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
        Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
        // Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_test);
        if (MO_new.size() == 0) { break;}
        
        //input solve
        Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);
        ret &= Bmatch_PruneInputSolverByStrSupport(pMan, MO_new);
        ret &= Bmatch_PruneInputSolverByFuncSupport(pMan, MO_new);
        ret &= Bmatch_PruneInputSolverBySymmetryProperty(pMan, MO_new);
        ret &= Bmatch_PruneInputSolverByUnate(pMan, MO_new);
        ret &= Bmatch_ApplyInputSolverRowConstraint(pMan, pNtk1, pNtk2);

        while (ret && result.status != EQUIVALENT && iter++ < maxIter) {
            ret &= Bmatch_PruneInputSolverByCounterPart(pMan, pNtk1, pNtk2, result.model, MI, MO_new);
            ret &= Bmatch_PruneInputSolverBySymmetry(pMan, MI);
            if (!ret) break;
            auto Mapping = Bmatch_SolveInput(pMan, pNtk1, pNtk2, NULL, NULL, 0);
            if (Mapping.status == 0) break;
            MI = Mapping.MI;

            assert(MI.size() == Abc_NtkPiNum(pNtk1) + 1);
            result = Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO_new, 1, 0);
        }

        Abc_PrintTime(1, "Current time", Abc_Clock() - clkTotal);

        if (result.status == EQUIVALENT) {
            if (option & VERBOSE_DETAIL_MASK) printf("Find matching at iteration %d!!!\n", iter);
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
        Bmatch_New_Or(pMan, n, m, 0);

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

bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m, int verbose){
    std::vector<int> LearnedLevel = pMan->LearnedLevel;
    auto *pSolver = pMan->pOutputSolver;
    std::vector<int> LearnedAssumption = pMan->LearnedAssumption;
    vMatch_Group &MO = pMan->MO;
    std::cout<<"start backtrack"<<std::endl;
    std::cout<<"learn level:"<<LearnedLevel.size()<<std::endl;

    if (verbose){

        for(int i =0; i<LearnedLevel.size(); i++){
            std::cout<<LearnedLevel[i]<<" ";
        }
        std::cout<<"learn assumption:"<<LearnedAssumption.size()<<std::endl;
        for(int i =0; i<LearnedAssumption.size(); i++){
            std::cout<<LearnedAssumption[i]<<" ";
        }
        std::cout<<std::endl;
        std::cout<<"pop ";
        for(int i = 0; i < LearnedLevel.back(); i++){
            std::cout<<LearnedAssumption.back()<<" ";
            LearnedAssumption.pop_back();
        }
        std::cout<<std::endl;

        std::cout<<"new assumption:"<<LearnedAssumption.size()<<std::endl;
        for(int i =0; i<LearnedAssumption.size(); i++){
            std::cout<<LearnedAssumption[i]<<" ";
        }
        std::cout<<std::endl;
    }

    
    LearnedLevel.pop_back();
    if (LearnedLevel.size() == 0) LearnedLevel.emplace_back(1);
    else LearnedLevel.back()++;

    if (verbose){

        std::cout<<"new learn level:"<<LearnedLevel.size()<<std::endl;
        for(int i =0; i<LearnedLevel.size(); i++){
            std::cout<<LearnedLevel[i]<<" ";
        }
        std::cout<<std::endl;
    }
    if (MO.size() != 0){
        LearnedAssumption.emplace_back(Bmatch_toLitCond(MO.back()[0], 1));
        MO.pop_back();
        pMan->LearnedLevel = LearnedLevel;
        pMan->LearnedAssumption = LearnedAssumption;
        Bmatch_New_Or(pMan, n, m, 0);
        
        
    }
    else{
        printf("all possible path traced\n");
        return false;
    }

    int test = 0;
    for(int i = 0; i<LearnedLevel.size(); i++){
        test+= LearnedLevel[i];
    }
    assert(test == LearnedAssumption.size());

    pMan->pOutputSolver = pSolver;
    pMan->LearnedLevel = LearnedLevel;
    pMan->LearnedAssumption = LearnedAssumption;
    return true;
}

void Bmatch_New_Or(Bmatch_Man_t *pMan, int n, int m, int verbose){
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
                if (verbose) std::cout<<i*m+j<<" ";
            }
        }
    }
    std::cout<<std::endl;
    
    Bmatch_sat_solver_addvar(pSolver);
    pLits[m*n - LearnedAssumption.size()] = Bmatch_toLit(Bmatch_sat_solver_nvars(pSolver)-1);
    if (verbose){
        for(int i = 0;i<n*m; i++){
            std::cout<<pLits[i]<<" ";

        }
        std::cout<<std::endl;
    }
    pMan->ClauseControl.emplace_back(Bmatch_sat_solver_nvars(pSolver)-1);

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

        if (status == 20){
            if (fVerbose) std::cout<<"projection off"<<std::endl;
            bool endloop = Bmatch_OutputBacktrack(pMan, n, m, 0);
            vMatch temp;
            if (!endloop) return temp;

            learnedAssumption = pMan->LearnedAssumption;
            ClauseControl = pMan->ClauseControl;
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

            if (fVerbose){
                std::cout<<"backtrack debug"<<" clauseControl:"<<ClauseControl.size()<<std::endl;
                std::cout<<pLit[LitSize-1]<<" "<<pLit[LitSize-2]<<" "<<pLit[LitSize-3]<<std::endl;
                std::cout<<ClauseControl[ClauseControl.size()-1]<<" "<<ClauseControl[ClauseControl.size()-2];
            }
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
            if (fVerbose) printf(" %3d", Bmatch_sat_solver_var_value(pSolver, i * m + j) > 0);
            if (Bmatch_sat_solver_var_value(pSolver, i * m + j) > 0) {
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

InputMapping Bmatch_HeuristicSolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());

}

InputMapping Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    auto *pSolver = pMan->pInputSolver;
    // pSolver->set("verbose", 3);

    int status = Bmatch_sat_solver_simplify(pSolver, 0);
    if (status == 20) return {0, MI};
    status = Bmatch_sat_solver_solve(pSolver, bLits, eLits, 0, 0, 0, 0);
    if (status == 20 || status == 0) return {0, MI};

    int n = pMan->ni;
    int m = pMan->mi;

    if (fVerbose) { printf("       "); for (int j = 0; j < m - 2; ++j) { printf("%c%-3s", ((j & 1) != 0) ? '~' : ' ', Abc_ObjName(Abc_NtkPi(pNtk1, j / 2))); } printf(" 1   0\n"); }
    for (int i = 0; i < n; ++i) {
        if (fVerbose) printf("%3s: ", Abc_ObjName(Abc_NtkPi(pNtk2, i)));
        for (int j = 0; j < m; ++j) {
            if (fVerbose) printf(" %3d", Bmatch_sat_solver_var_value(pSolver, i * m + j) > 0);
            if (Bmatch_sat_solver_var_value(pSolver, i * m + j) > 0) {
                MI[j / 2].emplace_back(i, (int)((j & 1) != 0));
            }
        }
        if (fVerbose) printf("\n");
    }
    if (fVerbose) printf("\n");

    return {1, MI};
}

int Bmatch_PruneInputSolverByUnate(Bmatch_Man_t *pMan, vMatch &MO) {
    auto *pSolver = pMan->pInputSolver;

    Mat &unateMat1 = pMan->unateMat1;
    Mat &unateMat2 = pMan->unateMat2;
    auto &impossibleMI = pMan->impossibleMI;

    int ret = 1;

    int mi = pMan->mi;
    int ni = pMan->ni;

    int Lit;

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (int xi = 0; xi < mi / 2 - 1; ++xi) {
                int unateness1 = unateMat1[fi][xi];

                for (int yi = 0; yi < ni; ++yi) {
                    int unateness2 = unateMat2[gi][yi];

                    if (unateness1 == -1 || unateness2 == -1 || unateness1 == 3 || unateness2 == 3) continue;

                    if (unateness1 == unateness2) {
                        Lit = Bmatch_toLitCond(yi * mi + xi * 2 + (1 - g.sign()), 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        // printf("(%d, %d) ", yi, xi * 2 + 1);
                    } else if (unateness1 != unateness2) {
                        Lit = Bmatch_toLitCond(yi * mi + xi * 2 + g.sign(), 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        // printf("(%d, %d) ", yi, xi * 2);
                    }
                }
            }
            // printf("\n");
        }
    }

    return ret;
}

int Bmatch_PruneInputSolverBySymmetryProperty(Bmatch_Man_t *pMan, vMatch &MO) {
    if (pMan->vSymm1.empty() || pMan->vSymm2.empty()) return 1;
    auto *pSolver = pMan->pInputSolver;

    vSymm &symm1 = pMan->vSymm1;
    vSymm &symm2 = pMan->vSymm2;

    int ret = 1;
    int mi = pMan->mi;

    AutoBuffer<int> pLits(pMan->mi * pMan->ni);

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &g : MO[fi]) {
            int gi = g.var();

            for (auto &vSymm1 : symm1[fi]) {      // vector of fi input symmetry group
                for (auto &vSymm2 : symm2[gi]) {  // vector of gi input symmetry group
                    for (auto &xi : vSymm1) {
                        for (auto &yi : vSymm2) {
                            for (int j = 0; j < 2; ++j) {
                                int start = 0;
                                pLits[start++] = Bmatch_toLitCond(yi * mi + xi * 2 + j, (1 - j));
                                pLits[start++] = Bmatch_toLitCond(yi * mi + xi * 2 + (1 - j), j);
                                // printf("%c(%d, %d) ", (1 - j) ? '!' : ' ', yi, xi * 2 + j);
                                // printf("%c(%d, %d) ", (j) ? '!' : ' ', yi, xi * 2 + (1 - j));
                                for (auto &xir : vSymm1) {
                                    if (xi == xir) continue;
                                    for (auto &yia : vSymm2) {
                                        if (yi == yia) continue;
                                        pLits[start++] = Bmatch_toLit(yia * mi + xir * 2);
                                        pLits[start++] = Bmatch_toLit(yia * mi + xir * 2 + 1);
                                        // printf("(%d, %d) ", yia, xir * 2);
                                        // printf("(%d, %d) ", yia, xir * 2 + 1);
                                    }
                                }
                                // printf("\n");
                                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + start);
                            }
                        }
                    }
                }
            }
            break;
        }
    }

    // ABC_FREE(pLits);

    return ret;
}

int Bmatch_PruneInputSolverBySymmetry(Bmatch_Man_t *pMan, vMatch &MI) {
    if (pMan->vSymmPair1.empty() || pMan->vSymmPair2.empty()) return 1;
    if (MI.empty()) return 1;
    int ret = 1;

    auto *pSolver = pMan->pInputSolver;

    vSymmPair &symm2 = pMan->vSymmPair2;

    int mi = pMan->mi;

    int pLits[2];
    for (int i = 0; i < MI.size(); ++i) {
        for (auto &p : MI[i]) {
            for (auto &s : symm2) { // loop symmetry pairs
                if ((p.var() == s.first || p.var() == s.second)) {
                    int symm_port = (p.var() == s.first) ? s.second : s.first;
                    int find = 0;
                    for (int j = 0; j < MI.size(); ++j) {
                        for (auto &q : MI[j]) {
                            if (q.var() == symm_port) { // find the corresponding pair
                                find = 1;
                                pLits[0] = Bmatch_toLitCond(p.var() * mi + j * 2 + q.sign(), 1);
                                pLits[1] = Bmatch_toLitCond(q.var() * mi + i * 2 + p.sign(), 1);
                                // printf("(%d, %d) ", p.var(), j * 2 + q.sign());
                                // printf("(%d, %d)\n", q.var(), i * 2 + p.sign());
                                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
                                break;
                            }
                        }
                        if (find) break;
                    }
                }
            }
        }
    }

    return ret;
}

int Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI, vMatch& MO) {
    if (!pModel1) return 1;

    int ret = 1;
    auto *pSolver = pMan->pInputSolver;
    AutoBuffer<int> pModel2(Abc_NtkPiNum(pNtk2));
    // get the input pattern of pNtk2
    for (int i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
        for (auto &p : MI[i]) {
            pModel2[p.var()] = (p.sign()) ? !pModel1[i] : pModel1[i];
        }
    }
    for (auto &p : MI[Abc_NtkPiNum(pNtk1)]) {
        pModel2[p.var()] = (p.sign()) ? 0 : 1;
    }

    int mi = pMan->mi;
    int k = 0;
    AutoBuffer<int> pLits((Abc_NtkPiNum(pNtk1) + 2) * Abc_NtkPiNum(pNtk2));

    auto &sRedund1 = pMan->sRedund1;
    auto &sRedund2 = pMan->sRedund2;
    auto &impossibleMI = pMan->impossibleMI;

    // counter example
    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        if (sRedund2[i]) continue;
        for (int j = 0; j < Abc_NtkPiNum(pNtk1); ++j) {
            if (sRedund1[j]) continue;
            else if (impossibleMI[i * mi + j * 2] || impossibleMI[i * mi + j * 2 + 1]) continue;
            pLits[k++] = Bmatch_toLit(i * mi + j * 2 + (pModel2[i] == pModel1[j]));
            // printf("(%d %d) ", i, j * 2 + (pModel2[i] == pModel1[j]));
        }
        pLits[k++] = Bmatch_toLit((i * mi) + mi - 1 - (pModel2[i] == 0)); // don't know how this work???
    }
    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + k);

    // discard current
    // k = 0;
    // for (int i = 0; i < Abc_NtkPiNum(pNtk1) + 1; ++i) {
    //     for (auto &p : MI[i]) {
    //         pLits[k++] = Bmatch_toLitCond(p.var() * mi + i * 2 + p.sign(), 1);
    //     }
    // }
    // Bmatch_sat_solver_addclause(pSolver, pLits, pLits + k);

    ABC_FREE(pModel1);
    // ABC_FREE(pModel2);
    // ABC_FREE(pLits);

    return ret;
}

int Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oFuncSupp1;
    auto &oSupp2 = pMan->oFuncSupp2;
    auto &impossibleMI = pMan->impossibleMI;

    AutoBuffer<int> pLits(ni * 2);

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];       // functional support of fi
            auto &cond2 = oSupp2[gi.var()]; // functional support of gi

            for (auto &m : cond1) {         // every support(fi)
                int i = 0;
                for (auto &n : cond2) {     // every support(gi)
                    pLits[i++] = Bmatch_toLit(n * mi + m * 2);
                    pLits[i++] = Bmatch_toLit(n * mi + m * 2 + 1);
                    if (cond1.size() == cond2.size()) { // if functional support is the same, then the input of gi cannot be const
                        int Lit = Bmatch_toLitCond(n * mi + mi - 1, 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        impossibleMI[n * mi + mi - 1] = 1;
                        Lit = Bmatch_toLitCond(n * mi + mi - 2, 1);
                        Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                        impossibleMI[n * mi + mi - 2] = 1;
                    }
                }
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + i);
            }
        }
    }

    // ABC_FREE(pLits);

    return ret;
}

int Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;
    auto &impossibleMI = pMan->impossibleMI;

    AutoBuffer<int> valid(ni * mi, 0);
    
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    valid[n * mi + m * 2] = 1;
                    valid[n * mi + m * 2 + 1] = 1;
                }
            }
        }
    }

    for (int n = 0; n < ni; ++n) {
        for (int m = 0; m < mi - 2; ++m) {
            if (valid[n * mi + m] == 0) {
                int Lit = Bmatch_toLitCond(n * mi + m, 1);
                Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                impossibleMI[n * mi + m] = 1;
            }
        }
    }

    // ABC_FREE(valid);

    return ret;
}

int Bmatch_ApplyInputSolverRowConstraint(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;

    int start = 0;
    int n = pMan->ni;
    int m = pMan->mi;
    auto &sRedund1 = pMan->sRedund1;
    auto &sRedund2 = pMan->sRedund2;
    auto &impossibleMI = pMan->impossibleMI;

    AutoBuffer<int> pLits(m);

    for (int i = 0; i < n; ++i) {
        // if input of cir2 is redundant
        if (sRedund2[i]) {
            pLits[0] = Bmatch_toLit(i * m + m - 1); // set it to constant
            Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
            continue;
        }
        // V(aij v bij)
        start = 0;
        for (int j = 0; j < m; ++j) {
            if (j < m - 2 && sRedund1[j / 2]) continue; // do not consider cir1's redundant input
            else if (impossibleMI[i * m + j]) continue; // impossible MI
            pLits[start++] = Bmatch_toLit(i * m + j);
        }
        Bmatch_sat_solver_addclause(pSolver, pLits, pLits + start);

        // VC~(row, 2)
        for (int j = 0; j < m - 1; ++j) {
            if (j < m - 2 && sRedund1[j / 2]) continue; // do not consider cir1's redundant input
            else if (impossibleMI[i * m + j]) continue; // impossible MI
            pLits[0] = Bmatch_toLitCond(i * m + j, 1);
            for (int k = j + 1; k < m; ++k) {
                pLits[1] = Bmatch_toLitCond(i * m + k, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
            }
        }
    }

    return ret;
}

int Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int ret = 1;
    auto *pSolver = pMan->pInputSolver;
    if (pSolver) Bmatch_sat_solver_delete(pSolver);
    pSolver = Bmatch_sat_solver_new();
    // pSolver->fVerbose = 1;
    // pSolver->fPrintClause = 1;
    // pSolver->verbosity = 2;

    int n = pMan->ni = Abc_NtkPiNum(pNtk2);
    int m = pMan->mi = 2 * (Abc_NtkPiNum(pNtk1) + 1);

    Bmatch_sat_solver_setnvars(pSolver, n * m);

    pMan->heuristicStage = 0;
    pMan->impossibleMI.resize(n * m);
    pMan->impossibleMI.fill(0);
    pMan->pInputSolver = pSolver;

    return ret;
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
