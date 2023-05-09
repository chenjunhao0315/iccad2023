#include "bmatch.hpp"
#include <vector>

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

void Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO);
void Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO);
void Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI);
vMatch Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);

vMatch Bmatch_SolveOutput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);
void Bmatch_InitOutputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan);
void Bmatch_OutputLearn(Bmatch_Man_t *pMan, bool status, int n, int m);
bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m);
void Bmatch_New_Or(Bmatch_Man_t *pMan, int n, int m);

#ifdef __cplusplus
}
#endif

// vMatch MI = {{Literal(1, false), Literal(3, false)}, {Literal(2, true)}, {Literal(0, true)}, {}};
// vMatch MO = {{Literal(0, false)}};

// vMatch MI = {{Literal(2, false)}, {Literal(3, false)}, {Literal(0, false)}, {Literal(1, false)}, {Literal(4, false)}, {}};
// vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {}, {}};

// vMatch MI = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}, {Literal(4, false)}, {}};
// vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}};

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    int maxIter = 50, iter = 0;
    

    if (option & VERBOSE_MASK) Bmatch_PrintBusInfo(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintInputSense(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintOutputSupport(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintSymm(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintEqual(pMan, pNtk1, pNtk2);

    //preprocess
    Bmatch_InitInputSolver(pMan, pNtk1, pNtk2);
    Bmatch_SolveOutputGroup(pMan);
    Bmatch_InitOutputSolver(pMan, pNtk1, pNtk2);

    if (option & VERBOSE_MASK) Bmatch_PrintOutputGroup(pNtk1, pNtk2, pMan->Groups);
    


    vMatch MI, MO_new;
    // vMatch_Group MO;
    bool optimal = false;
    while(!optimal){
        //find new pair of output matching
        EcResult result;
        MO_new = Bmatch_SolveOutput(pMan, pNtk1, pNtk2, NULL, NULL, 1);
        Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
        if (MO_new.size() == 0) return;
        
        //input solve
        Bmatch_PruneInputSolverByStrSupport(pMan, MO_new);
        Bmatch_PruneInputSolverByFuncSupport(pMan, MO_new);
        while (result.status != EQUIVALENT && iter++ < maxIter) {
            Bmatch_PruneInputSolverByCounterPart(pMan, pNtk1, pNtk2, result.model, MI);
            MI = Bmatch_SolveInput(pMan, pNtk1, pNtk2, NULL, NULL, 0);

            assert(MI.size() == Abc_NtkPiNum(pNtk1) + 1);
            // assert(MO.size() == Abc_NtkPoNum(pNtk1));
            result = Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO_new, 0);
            print("Miter Status:", (result.status == EQUIVALENT) ? "EQUIVALENT" : "NON-EQUIVALENT");
        }

        if (result.status == EQUIVALENT) {
            printf("Find matching at iteration %d!!!\n", iter);
            Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO_new);
            Bmatch_OutputLearn(pMan, true, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
        } else {
            printf("Reach maximum iteration!\n");
            Bmatch_OutputLearn(pMan, false, Abc_NtkPoNum(pNtk2), 2*Abc_NtkPoNum(pNtk1));
            
        }
        iter = 0;
        optimal = result.status == EQUIVALENT;
    }
}

void Bmatch_OutputLearn(Bmatch_Man_t *pMan, bool status, int n, int m){

    std::vector<int> &LearnedLevel = pMan->LearnedLevel;
    sat_solver *pSolver = pMan->pOutputSolver;  
    std::vector<int> &LearnedAssumption = pMan->LearnedAssumption;
    vMatch_Group &MO = pMan->MO;

    if(status){
        //add new assumption base on i level(c)
        auto &MoBack=MO.back();
        LearnedLevel.emplace_back(MoBack.size()+1);
        for(auto &match:MoBack){
            // sat_solver_push(pSolver, toLit(match));
            
            LearnedAssumption.emplace_back(toLit(match));
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
            LearnedAssumption.emplace_back(toLitCond(match, 1));
        }
        MO.pop_back();

    }
    pMan->pOutputSolver = pSolver;
    pMan->LearnedLevel = LearnedLevel;
    pMan->LearnedAssumption = LearnedAssumption;
}

bool Bmatch_OutputBacktrack(Bmatch_Man_t *pMan, int n, int m){
    std::vector<int> LearnedLevel = pMan->LearnedLevel;
    sat_solver *pSolver = pMan->pOutputSolver;
    std::vector<int> LearnedAssumption = pMan->LearnedAssumption;
    vMatch_Group &MO = pMan->MO;

    for(int i = 0; i < LearnedLevel.back(); i++){
        LearnedAssumption.pop_back();
    }
    
    LearnedLevel.pop_back();
    if (LearnedLevel.size() == 0) LearnedLevel.emplace_back(1);
    else LearnedLevel.back() += 1;
    if (MO.size() != 0){
        LearnedAssumption.emplace_back(toLitCond(MO.back()[0], 1));
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
    sat_solver *pSolver = pMan->pOutputSolver;

    int *pLits = ABC_ALLOC(int, n*m);
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            if(std::find(LearnedAssumption.begin(), LearnedAssumption.end(), i*m+j) == LearnedAssumption.end()){
                pLits[i*m+j] = toLit(i*m+j);
            }
        }
    }
    sat_solver_addvar(pSolver);
    pLits[m*n-LearnedAssumption.size()] = toLit(sat_solver_nvars(pSolver));
    pMan->ClauseControl.emplace_back(sat_solver_nvars(pSolver));
    sat_solver_addclause(pSolver, pLits, pLits + m*n - LearnedAssumption.size() + 1);
}

vMatch Bmatch_SolveOutput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose){
    
    vMatch MO_new(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
    sat_solver *pSolver = pMan->pOutputSolver;
    std::vector<int> learnedAssumption = pMan->LearnedAssumption;
    std::vector<int> ClauseControl = pMan->ClauseControl;
    pSolver->verbosity = 2;
    int n = Abc_NtkPoNum(pNtk2);
    int m = 2 * (Abc_NtkPoNum(pNtk1));
    int *model = pSolver->model;
    int LitSize = learnedAssumption.size()+ClauseControl.size()+1;
    vMatch_Group &MO = pMan->MO;


    //add assumption
    int *pLit = ABC_ALLOC(int, learnedAssumption.size()+ClauseControl.size()+1);
    for(int i =0; i<learnedAssumption.size(); i++){
        pLit[i] = learnedAssumption[i];
        // std::cout<<"assump"<<learnedAssumption[i]<<std::endl;
    }
    //clause control
    for(int i = 0; i<ClauseControl.size(); i++){
        if (i == ClauseControl.size()-1) pLit[learnedAssumption.size()+i] = toLitCond(ClauseControl.back(), 1);
        else pLit[i+learnedAssumption.size()] = toLit(ClauseControl[i]);
        // std::cout<<"clause"<<ClauseControl[i];
    }
    
    // //projective
    if (pMan->AllowProjection) pLit[learnedAssumption.size()+ClauseControl.size()] = toLit(pMan->Projective);
    else pLit[learnedAssumption.size()+ClauseControl.size()] = toLitCond(pMan->Projective, 1);


    int status = sat_solver_solve(pSolver, pLit, pLit+LitSize, 0, 0, 0, 0);
    while(status == l_False){
        std::cout<<"output match failed"<<std::endl;
        if (status == l_False){
            //projection on
            pMan->AllowProjection = true;
            pLit[LitSize-1] = toLit(pMan->Projective);
            std::cout<<pMan->Projective<<std::endl;
            // std::cout<<LitSize<<" "<<learnedAssumption.size()<<" "<<ClauseControl.size()<<" "<<pMan->AllowProjection<<std::endl;
            status = sat_solver_solve(pSolver, pLit, pLit+LitSize, 0, 0, 0, 0);
        }
        std::cout<<"projection on"<<std::endl;
        if (status == l_False){
            bool endloop = Bmatch_OutputBacktrack(pMan, n, m);
            vMatch temp;
            if (!endloop) return temp;
            std::cout<<"backtrack end"<<std::endl;
            //projection off
            pMan->AllowProjection = false;
            // std::cout<<LitSize<<" "<<learnedAssumption.size()<<" "<<ClauseControl.size()<<" "<<pMan->AllowProjection<<std::endl;
            for(int i =0; i<learnedAssumption.size(); i++){
                pLit[i] = learnedAssumption[i];
            }
            for(int i = 0; i<ClauseControl.size(); i++){
                if (i == ClauseControl.size()-1) pLit[learnedAssumption.size()+i] = toLitCond(ClauseControl.back(), 1);
                else pLit[i+learnedAssumption.size()] = toLit(ClauseControl[i]);
            }
            pLit[LitSize-1] = toLitCond(pMan->Projective, 1);
            status = sat_solver_solve(pSolver, pLit, pLit+LitSize, 0, 0, 0, 0);
        }
    }

    std::vector<int> new_match;
    if (fVerbose) { printf("       "); for (int j = 0; j < m; ++j) { printf("%c%-3s", ((j & 1) != 0) ? '~' : ' ', Abc_ObjName(Abc_NtkPo(pNtk1, j / 2))); }; printf("\n");  }
    for (int i = 0; i < n; ++i) {
        if (fVerbose) printf("%3s: ", Abc_ObjName(Abc_NtkPo(pNtk2, i)));
        for (int j = 0; j < m; ++j) {
            if (fVerbose) printf(" %3d", model[i * m + j] == l_True);
            if (model[i * m + j] == l_True) {
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
                printf(" %3d", level[j], " ");
            }
            printf("\n");
        }
    }

    
    pMan->MO = MO;
    ABC_FREE(pLit);
    return MO_new;

}
void Bmatch_InitOutputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2){
    sat_solver *pSolver = pMan->pOutputSolver;
    if (pSolver) sat_solver_delete(pSolver);
    pSolver = sat_solver_new();
    pSolver->verbosity = 2;

    int n = Abc_NtkPoNum(pNtk2);
    int m = 2 * (Abc_NtkPoNum(pNtk1));
    int *pLits = ABC_ALLOC(int, m*n);
    
    sat_solver_setnvars(pSolver, n * m + 1);
    //at least one pair
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            pLits[i*m+j] = toLit(i*m+j);
        }
    }
    sat_solver_addclause(pSolver, pLits, pLits + m*n);

    //(c+d)<=1 ==> (~c+~d)
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            pLits[0] = toLitCond(i*m+j, 1);
            for(int k = j+1; k<m; k++){
                if (k != j){
                    pLits[1] = toLitCond(i*m+k, 1);
                    sat_solver_addclause(pSolver, pLits, pLits + 2);
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
                    pLits[0] = toLitCond(Group_ntk1[i]*2 + m*Group_ntk2[j],1);
                    sat_solver_addclause(pSolver, pLits, pLits + 1);
                    pLits[0] = toLitCond(Group_ntk1[i]*2 + m*Group_ntk2[j] + 1,1);
                    sat_solver_addclause(pSolver, pLits, pLits + 1);
                }
            }
            
        }
    }
    
    // allow projection
    pMan->AllowProjection = false;
    pMan->Projective = n*m+1;
    pLits[0] = toLit(n*m+1);//allow projection
    for(int i = 0; i<n; i++){
        for(int j = i+1; i<n; i++){
            for(int k =0; k<m/2; k++){
                pLits[1] = toLitCond(i*m+k*2, 1);
                pLits[2] = toLitCond(j*m+k*2+1, 1);
                sat_solver_addclause(pSolver, pLits, pLits+3);

                pLits[1] = toLitCond(i*m+k*2+1, 1);
                pLits[2] = toLitCond(j*m+k*2, 1);
                sat_solver_addclause(pSolver, pLits, pLits+3);

                pLits[1] = toLitCond(i*m+k*2, 1);
                pLits[2] = toLitCond(j*m+k*2, 1);
                sat_solver_addclause(pSolver, pLits, pLits+3);

                pLits[1] = toLitCond(i*m+k*2+1, 1);
                pLits[2] = toLitCond(j*m+k*2+1, 1);
                sat_solver_addclause(pSolver, pLits, pLits+3);
            }
        }
    }

    //equal


    ABC_FREE(pLits);
    pMan->pOutputSolver = pSolver;

}

vMatch Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose) {
    vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    sat_solver *pSolver = pMan->pInputSolver;

    int status = sat_solver_solve(pSolver, bLits, eLits, 0, 0, 0, 0);
    print(status == l_True ? "InputSolver: SAT" : "InputSolver: UNSAT");
    int *model = pSolver->model;
    if (status == l_False) return MI;

    int n = pMan->ni;
    int m = pMan->mi;

    if (fVerbose) { printf("       "); for (int j = 0; j < m - 2; ++j) { printf("%c%-3s", ((j & 1) != 0) ? '~' : ' ', Abc_ObjName(Abc_NtkPi(pNtk1, j / 2))); } printf(" 1   0\n"); }
    for (int i = 0; i < n; ++i) {
        if (fVerbose) printf("%3s: ", Abc_ObjName(Abc_NtkPi(pNtk2, i)));
        for (int j = 0; j < m; ++j) {
            if (fVerbose) printf(" %3d", model[i * m + j] == l_True);
            if (model[i * m + j] == l_True) {
                MI[j / 2].emplace_back(i, (int)((j & 1) != 0));
            }
        }
        if (fVerbose) printf("\n");
    }

    return MI;
}

void Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI) {
    if (!pModel1) return;

    sat_solver *pSolver = pMan->pInputSolver;
    int *pModel2 = ABC_ALLOC(int, Abc_NtkPiNum(pNtk2));
    // get the input pattern of pNtk2
    for (int i = 0; i < Abc_NtkPiNum(pNtk1); ++i) {
        for (auto &p : MI[i]) {
            pModel2[p.var()] = (p.sign()) ? !pModel1[i] : pModel1[i];
        }
    }

    int mi = pMan->mi;
    int k = 0;
    int *pLits = ABC_ALLOC(int, Abc_NtkPiNum(pNtk1) * Abc_NtkPiNum(pNtk2));

    for (int i = 0; i < Abc_NtkPiNum(pNtk2); ++i) {
        for (int j = 0; j < Abc_NtkPiNum(pNtk1); ++j) {
            pLits[k++] = toLit(i * mi + j * 2 + (pModel2[i] == pModel1[j]));
        }
    }

    sat_solver_addclause(pSolver, pLits, pLits + k);

    ABC_FREE(pModel1);
    ABC_FREE(pModel2);
    ABC_FREE(pLits);
}

void Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    sat_solver *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oFuncSupp1;
    auto &oSupp2 = pMan->oFuncSupp2;

    int *pLits = ABC_ALLOC(int, ni * 2);

    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &m : cond1) {
                int i = 0;
                for (auto &n : cond2) {
                    pLits[i++] = toLit(n * mi + m * 2);
                    pLits[i++] = toLit(n * mi + m * 2 + 1);
                }
                sat_solver_addclause(pSolver, pLits, pLits + i);
            }
        }
    }

    ABC_FREE(pLits);
}

void Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO) {
    sat_solver *pSolver = pMan->pInputSolver;

    int ni = pMan->ni;
    int mi = pMan->mi;

    auto &oSupp1 = pMan->oStrSupp1;
    auto &oSupp2 = pMan->oStrSupp2;

    int *valid = ABC_CALLOC(int, ni * mi);
    
    for (int fi = 0; fi < MO.size(); ++fi) {
        for (auto &gi : MO[fi]) {
            auto &cond1 = oSupp1[fi];
            auto &cond2 = oSupp2[gi.var()];

            for (auto &n : cond2) {
                for (auto &m : cond1) {
                    valid[n * mi + m * 2] = 1;
                    valid[n * mi + m * 2 + 1] = 1;
                    
                }
                valid[n * mi + mi - 2] = 1;
                valid[n * mi + mi - 1] = 1;
            }
        }
    }

    for (int n = 0; n < ni; ++n) {
        for (int m = 0; m < mi; ++m) {
            if (valid[n * mi + m] == 0) {
                int Lit = toLitCond(n * mi + m, 1);
                sat_solver_addclause(pSolver, &Lit, &Lit + 1);
            }
        }
    }

    ABC_FREE(valid);
}

void Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    sat_solver *pSolver = pMan->pInputSolver;
    if (pSolver) sat_solver_delete(pSolver);
    pSolver = sat_solver_new();
    // pSolver->fVerbose = 1;
    // pSolver->fPrintClause = 1;
    // pSolver->verbosity = 2;

    int n = pMan->ni = Abc_NtkPiNum(pNtk2);
    int m = pMan->mi = 2 * (Abc_NtkPiNum(pNtk1) + 1);

    int *pLits = ABC_ALLOC(int, m);

    sat_solver_setnvars(pSolver, n * m);
    for (int i = 0; i < n; ++i) {
        // V(aij v bij)
        for (int j = 0; j < m; ++j) {
            pLits[j] = toLit(i * m + j);
            // std::cout<<pLits[j]<<" ";
        }
        sat_solver_addclause(pSolver, pLits, pLits + m);
        // std::cout<<std::endl;

        // VC~(row, 2)
        for (int j = 0; j < m - 1; ++j) {
            pLits[0] = toLitCond(i * m + j, 1);
            for (int k = j + 1; k < m; ++k) {
                pLits[1] = toLitCond(i * m + k, 1);
                sat_solver_addclause(pSolver, pLits, pLits + 2);
            }
        }
    }

    ABC_FREE(pLits);
    pMan->pInputSolver = pSolver;
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