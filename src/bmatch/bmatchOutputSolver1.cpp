#include "bmatch.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_OutputLearn(Bmatch_Man_t *pMan, bool status, int n, int m);
bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m, int verbose);
void Bmatch_New_Or(Bmatch_Man_t *pMan, int n, int m, int verbose);
int Bmatch_PruneOutputSolverByUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_InitOutputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
vMatch Bmatch_SolveOutput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);
#ifdef __cplusplus
}
#endif


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
            // std::cout<<Bmatch_toLit(match)<<std::endl;
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
        // LearnedAssumption.emplace_back(Bmatch_toLitCond(MoBack[0], 1));
        for(auto &match:MoBack){
                // std::cout<<"check"<<Bmatch_toLitCond(match, 1)<<std::endl;
                LearnedAssumption.emplace_back(Bmatch_toLitCond(match, 1));
        }
        // if(MoBack.size() == 1) LearnedAssumption.emplace_back(Bmatch_toLitCond(MoBack[0], 1));
        // else{
        //     for(auto &match:MoBack){
        //         // std::cout<<Bmatch_toLitCond(match, 1)<<std::endl;
        //         LearnedAssumption.emplace_back(Bmatch_toLitCond(match, 1));
        //     }
        // }
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

    for(int i = 0; i < LearnedLevel.back(); i++){
            // std::cout<<LearnedAssumption.back()<<" ";
            LearnedAssumption.pop_back();
        }

    
    LearnedLevel.pop_back();
    if (LearnedLevel.size() == 0) LearnedLevel.emplace_back(1);
    else LearnedLevel.back()++;

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
    
    AutoBuffer<int> pLits(n*m+1, 0);
    int cont = 0;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            if(std::find(LearnedAssumption.begin(), LearnedAssumption.end(), Bmatch_toLit(i*m+j)) == LearnedAssumption.end() && \
                std::find(LearnedAssumption.begin(), LearnedAssumption.end(), Bmatch_toLitCond(i*m+j, 1)) == LearnedAssumption.end()){
                pLits[cont] = Bmatch_toLit(i*m+j);
                cont++;
                if (verbose) std::cout<<Bmatch_toLit(i*m+j)<<" ";
            }
        }
    }
    if (verbose) std::cout<<std::endl;
    
    Bmatch_sat_solver_addvar(pSolver);
    pLits[m*n - LearnedAssumption.size()] = Bmatch_toLit(Bmatch_sat_solver_nvars(pSolver)-1);
    if(verbose){
        for(auto &i:LearnedAssumption){
            std::cout<<i<<" ";
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
            learnedAssumption = pMan->LearnedAssumption;
            ClauseControl = pMan->ClauseControl;

            vMatch temp;
            if (!endloop) return temp;
            //projection off
            pMan->AllowProjection = false;
            
            for(int i =0; i<learnedAssumption.size(); i++){
                pLit[i] = learnedAssumption[i];
            }

            for(int i = 0; i<ClauseControl.size(); i++){
                if (i == ClauseControl.size()-1) pLit[learnedAssumption.size()+i] = Bmatch_toLitCond(ClauseControl[i], 1);
                else pLit[i+learnedAssumption.size()] = Bmatch_toLit(ClauseControl[i]);
            }

            LitSize = ClauseControl.size()+learnedAssumption.size()+1;
            pLit[LitSize-1] = Bmatch_toLitCond(pMan->Projective, 1);

            if (fVerbose){
                std::cout<<"backtrack debug"<<" clauseControl:"<<ClauseControl.size()<<std::endl;
                std::cout<<pLit[LitSize-1]<<" "<<pLit[LitSize-2]<<" "<<pLit[LitSize-3]<<std::endl;
                std::cout<<ClauseControl[ClauseControl.size()-1]<<" "<<ClauseControl[ClauseControl.size()-2];
            }
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
                int Lit = Bmatch_toLitCond(j * (2 * Abc_NtkPoNum(pNtk2)) + i * 2, 1);
                // printf("(%d, %d) ", j, i * 2);
                Bmatch_sat_solver_addclause(pSolver, &Lit, &Lit + 1);
                Lit = Bmatch_toLitCond(j * (2 * Abc_NtkPoNum(pNtk2)) + i * 2 + 1, 1);
                // printf("(%d, %d) ", j, i * 2 + 1);
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
                pLits[1] = Bmatch_toLitCond(i*m+k, 1);
                Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
                // std::cout<<pLits[0]<<pLits[1]<<std::endl;
            }
        }
    }
    // std::cout<<std::endl;

    // pLits[0] = Bmatch_toLit(1);
    // Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
    // int status = Bmatch_sat_solver_solve(pSolver, NULL, NULL, 0, 0, 0, 0);
    // std::cout<<status<<std::endl;

    //output group
    vGroup vGroup = pMan->Groups;

    for(int l = 0; l< vGroup.size(); l++){
        auto Group_ntk1 = vGroup[l].first;
        for(int k = 0; k<vGroup.size(); k++){
            if (k == l) continue;
            auto Group_ntk2 = vGroup[k].second;
            for(int i = 0;i<Group_ntk1.size(); i++){
                for(int j = 0;j<Group_ntk2.size(); j++){
                    //for output_ntk1_i and output_ntk2_j not in same group
                    //add clause ~cij ~dij
                    pLits[0] = Bmatch_toLitCond(Group_ntk1[i]*2 + m*Group_ntk2[j],1);
                    // std::cout<<Group_ntk1[i]*2 + m*Group_ntk2[j]<<" ";
                    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
                    pLits[0] = Bmatch_toLitCond(Group_ntk1[i]*2 + m*Group_ntk2[j] + 1,1);
                    // std::cout<<Group_ntk1[i]*2 + m*Group_ntk2[j]+1<<" ";
                    Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 1);
                }
            }
            
        }
    }
    // std::cout<<std::endl;

    //test case6
    // int c1, c2 = 0;
    // std::cout<<c1<<" "<<c2<<std::endl;
    // while(c1<pMan->oPartition1.size() && c2<pMan->oPartition2.size()){
    //     if(pMan->oPartition1[c1].size() == 0 ){
    //         c1++;
    //     }

    //     if(pMan->oPartition2[c2].size() == 0 ){
    //         c2++;
    //     }
        
    //     if((pMan->oPartition1[c1].size() != 0) && (pMan->oPartition2[c2].size() != 0)){
    //         std::cout<<pMan->oPartition1[c1][0]<<" "<<pMan->oPartition2[c2][0]<<std::endl;
    //         pLits[0] = Bmatch_toLit(pMan->oPartition1[c1][0]*2 + pMan->oPartition2[c2][0]*m);
    //         pLits[1] = Bmatch_toLit(pMan->oPartition1[c1][0]*2 + pMan->oPartition2[c2][0]*m+1);
    //         Bmatch_sat_solver_addclause(pSolver, pLits, pLits + 2);
    //         c1++;
    //         c2++;
                        
    //     }
    //     // if(c1 != c2) break; 
    // }
    

    
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