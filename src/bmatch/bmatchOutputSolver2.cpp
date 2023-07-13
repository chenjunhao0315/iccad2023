#include "bmatch.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif
vMatch Bmatch_SolveOutput2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int fVerbose);
void Bmatch_PossibleOutputCalculate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int fVerbose);
void Bmatch_OutputLearn2(Bmatch_Man_t *pMan, bool status, int fVerbose);

#ifdef __cplusplus
}
#endif


void MatchTree::TreeConstruct(std::pair<std::vector<int>, std::vector<int> > Group, MatchNode *Node, int Level, int fVerbose){
    if(Level==Group.first.size()){
        if(fVerbose) std::cout<<"return"<<std::endl;
        return;
    }
    for(auto &g:Group.second){
        MatchNode gPNode = MatchNode(Group.first[Level], Literal(g, false));
        MatchNode gNNode = MatchNode(Group.first[Level], Literal(g, true));
        //sort g list
        std::vector<int> gNew;
        int temp;
        for(auto &gsort:Group.second){
            if(g == gsort) temp = gsort;
            else gNew.emplace_back(gsort);
        }
        gNew.emplace_back(temp);
        std::pair<std::vector<int>, std::vector<int> > GroupNew = std::make_pair(Group.first, gNew);
        if (fVerbose){
            // std::cout<<Level<<" "<<Group.first[Level]<<" "<<g<<"{";
            // for(auto &gsort:gNew){
            //     std::cout<<gsort<<" ";
            // }
            // std::cout<<"}"<<std::endl;
            // std::cout<<Level<<" "<<Group.first[Level]<<"->"<<g<<std::endl;
        }
        if(fVerbose) std::cout<<Level<<" "<<gPNode.fPort<<"->"<<gPNode.gPort.Var/2<<std::endl;
        TreeConstruct(GroupNew, &gPNode, Level+1, fVerbose);
        if(fVerbose) std::cout<<Level<<" "<<gNNode.fPort<<"->-"<<gNNode.gPort.Var/2<<std::endl;
        TreeConstruct(GroupNew, &gNNode, Level+1, fVerbose);
        Node->Child.emplace_back(gPNode);
        Node->Child.emplace_back(gNNode);
    }
    //f map to null
    // MatchNode gNode = MatchNode(Group.first[Level], Literal());
    // if(fVerbose) std::cout<<Level<<" "<<gNode.fPort<<"->"<<"-0"<<std::endl;
    // TreeConstruct(Group, &gNode, Level+1, fVerbose);
    // Node->Child.emplace_back(gNode);

}
void MatchTree::ResetTree(){
    for(auto &node:Root.Child){
        ResetTreeRecur(&node);
        node.ResetVisited();
    }
};
void MatchTree::TreePrint(){
    std::cout<<"root{"<<std::endl;
    for(auto &node:Root.Child){
        // std::cout<<node.fPort<<"->"<<node.gPort.Var<<std::endl;
        TreePrintRecur(&node);
    }
    std::cout<<"}"<<std::endl;
}
void MatchTree::TreePrintRecur(MatchNode *Node){
    std::cout<<Node->fPort<<"->"<<Node->gPort.Var<<std::endl;
    if(Node->Child.size() == 0) return;
    std::cout<<"{";
    for(auto &node:Node->Child){
        TreePrintRecur(&node);
    }
    std::cout<<"}"<<std::endl;
}


void MatchTree::ResetTreeRecur(MatchNode *Node){
    for(auto &n:Node->Child){
        if(n.Child.size() == 0){
            n.ResetVisited();
        }
        else {
            for(auto &nn:n.Child){
                ResetTreeRecur(&nn);
            }
            n.ResetVisited();
        }
    }
}
void MatchTree::FindNewMatchRecur(MatchNode *Node, MatchList *Result){
    if(Node->Child.size() == 0){
        Result->emplace_back(std::make_pair(Node->fPort, Node->gPort));
        // std::cout<<Node->fPort<<" "<<Node->gPort.Var/2<<std::endl;
        Node->SetVisited();
        return;
    }
    for(auto &node:Node->Child){
        if(!node.Visited) {
            int Current = Result->size();   
            FindNewMatchRecur(&node, Result);
            if(Result->size() == Current) node.SetVisited();
            else {
                if(Node->fPort != -1) Result->emplace_back(std::make_pair(Node->fPort, Node->gPort));
                // if(Node->fPort != -1) std::cout<<"check:"<<Node->fPort<<" "<<Node->gPort.Var/2<<std::endl;
                if(Node->Child.back().Visited) Node->SetVisited();
                return;
            }
        }
    }
}




vMatch Bmatch_SolveOutput2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int fVerbose){
    vMatch MO(Abc_NtkPoNum(pNtk1), std::vector<Literal>());
    // std::cout<<pMan->CurrentTree+1<<std::endl;
    for(int i = 0; i<pMan->CurrentTree+1; i++){
        MatchList Result;
        pMan->PossibleOutput[i].FindNewMatchRecur(&pMan->PossibleOutput[i].Root, &Result);
        if(fVerbose){
            for(auto &pair:Result){
                
                std::cout<<"pair:"<<pair.first<<pair.second.Var/2<<std::endl;
            }
        }
        for(auto &pair:Result){
            if(!pair.second.isUndef())MO[pair.first].emplace_back(pair.second);
        }
    }
    // Bmatch_PrintMatching(pNtk1, pNtk2, MO, MO); 
    return MO;
}

void Bmatch_PossibleOutputCalculate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int fVerbose){
    vGroup &group = pMan->Groups;
    PossList PossibleOutput;
    
    for(auto &Group:group){
        MatchTree Tree;
        //maybe add some pruning here 
        Tree.TreeConstruct(Group, &Tree.Root, 0, 0);
        PossibleOutput.emplace_back(Tree);
        if(fVerbose) std::cout<<"tree"<<PossibleOutput.size()<<" end"<<std::endl;
    }
    pMan->PossibleOutput = PossibleOutput;
    // if(fVerbose){
    //     for(auto &Tree:pMan->PossibleOutput){
    //             std::cout<<"tree"<<std::endl;
    //             Tree.TreePrint();
    //             std::cout<<std::endl; 
    //         }
    // }
}

void Bmatch_OutputLearn2(Bmatch_Man_t *pMan, bool status, int fVerbose){
    if(status) {
        pMan->CurrentTree++;
    }
    else{
        if(pMan->PossibleOutput[pMan->CurrentTree].Root.Visited){
            pMan->PossibleOutput[pMan->CurrentTree].ResetTree();
            pMan->CurrentTree--;
        }
    }
}
