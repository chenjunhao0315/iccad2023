#include "bmatch.hpp"
#include <map>
ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C"
{
#endif
    void Bmatch_PPBusPrune(Bmatch_Man_t *pMan, std::map<int, std::vector<int>> PossibleListI, std::map<int, std::vector<int>> PossibleListO, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
    void Bmatch_PPCheck(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);

    void Bmatch_DegreePrune(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
    void Bmatch_MintermPrune(Bmatch_Man_t *pMan, std::map<int, std::vector<int>> PossibleListI, std::map<int, std::vector<int>> PossibleListO, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
    void Bmatch_SignPrune(Bmatch_Man_t *pMan, std::map<int, std::vector<int>> PossibleListI, std::map<int, std::vector<int>> PossibleListO, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
    void Bmatch_MintermWCofactor(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int fVerbose);
#ifdef __cplusplus
}
#endif

void Bmatch_PPCheck(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2)
{
    // printf("degree prune\n");
    Bmatch_DegreePrune(pMan, pNtk1, pNtk2);

    std::map<int, std::vector<int>> PossibleListI, PossibleListO;
    Bmatch_PPBusPrune(pMan, PossibleListI, PossibleListO, pNtk1, pNtk2);
    Bmatch_PrintPartition(pMan, pNtk1, pNtk2);
    // Bmatch_MintermPrune(pMan, PossibleListI, PossibleListO, pNtk1, pNtk2);
    // Bmatch_SignPrune(pMan, PossibleListI, PossibleListO, pNtk1, pNtk2);
    // vMatch test;
    // return test;
}

void Bmatch_PPBusPrune(Bmatch_Man_t *pMan, std::map<int, std::vector<int>> PossibleListI, std::map<int, std::vector<int>> PossibleListO, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2)
{

#define BUS_TABLE(Table, Bus)                              \
    do                                                     \
    {                                                      \
        for (auto &i : Bus)                                \
        {                                                  \
            for (auto &j : i)                              \
            {                                              \
                Table.insert(std::make_pair(j, i.size())); \
            }                                              \
        }                                                  \
    } while (0)

    std::map<int, int> BusSizeTableI1, BusSizeTableI2, BusSizeTableO1, BusSizeTableO2;
    BUS_TABLE(BusSizeTableI1, pMan->BI1);
    BUS_TABLE(BusSizeTableI2, pMan->BI2);
    BUS_TABLE(BusSizeTableO1, pMan->BO1);
    BUS_TABLE(BusSizeTableO2, pMan->BO2);

    // possible list init
    std::vector<int> temp;
    for (int i = 0; i < pMan->iFuncSupp1.size(); i++)
    {
        PossibleListI[i] = temp;
    }
    for (int i = 0; i < pMan->oFuncSupp1.size(); i++)
    {
        PossibleListO[i] = temp;
    }

    // input pruning
    for (int i = 0; i < std::max(pMan->iPartition1.size(), pMan->iPartition2.size()); i++)
    {
        for (auto &j : pMan->iPartition1[i])
        {
            for (auto &k : pMan->iPartition2[i])
            {
                bool bus1 = true;
                bool bus2 = true;
                if (BusSizeTableI1.find(j) == BusSizeTableI1.end())
                    bus1 = false;
                if (BusSizeTableI2.find(k) == BusSizeTableI2.end())
                    bus2 = false;

                if (!bus1 && !bus2)
                {
                    PossibleListI[j].emplace_back(k);
                }
                else if (bus1 && bus2 && (BusSizeTableI1[j] <= BusSizeTableI2[k]))
                {
                    PossibleListI[j].emplace_back(k);
                }
            }
        }
    }
    // check
    // for (int i = 0; i < PossibleListI.size(); i++)
    // {
    //     std::cout << Abc_ObjName(Abc_NtkPi(pNtk1, i)) << ':';
    //     for (auto &j : PossibleListI[i])
    //     {
    //         std::cout << Abc_ObjName(Abc_NtkPi(pNtk2, j)) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // output pruning

    BusSizeTableI1.clear();
    BusSizeTableI2.clear();
    BusSizeTableO1.clear();
    BusSizeTableO2.clear();
}

void Bmatch_DegreePrune(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2)
{

#define DEGREE_PARTITION(FuncSupp, Partition)                     \
    do                                                            \
    {                                                             \
        std::map<int, int> buffer;                                \
        std::vector<int> temp;                                    \
        for (int i = 0, n = FuncSupp.size(); i < n; i++)          \
        {                                                         \
            buffer.insert(std::make_pair(i, FuncSupp[i].size())); \
        }                                                         \
        for (int i = 0, n = buffer.size(); i < n; i++)            \
        {                                                         \
            while (buffer[i] + 1 > Partition.size())              \
            {                                                     \
                Partition.emplace_back(temp);                     \
            }                                                     \
            Partition[buffer[i]].emplace_back(i);                 \
        }                                                         \
        buffer.clear();                                           \
    } while (0)

    DEGREE_PARTITION(pMan->iFuncSupp1, pMan->iPartition1);
    DEGREE_PARTITION(pMan->iFuncSupp2, pMan->iPartition2);
    DEGREE_PARTITION(pMan->oFuncSupp1, pMan->oPartition1);
    DEGREE_PARTITION(pMan->oFuncSupp2, pMan->oPartition2);

    // printf("degree end\n");
    // Bmatch_PrintPartition(pMan, pNtk1, pNtk2);
}

void Bmatch_SignPrune(Bmatch_Man_t *pMan, std::map<int, std::vector<int>> PossibleListI, std::map<int, std::vector<int>> PossibleListO, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2)
{

    std::map<int, std::vector<int>> signI1, signI2, signO1, signO2;

#define SIGN_INIT(FuncSupp, Partition, sign, pNtk)                                                          \
    do                                                                                                \
    {                                                                                                 \
        for (int i = 0; i < FuncSupp.size(); i++)                                                     \
        {                                                                                             \
            std::vector<int> temp;                                                                    \
            for (auto &j : FuncSupp[i])                                                         \
            {                                                                                         \
                for (int k = 0; k < Partition.size(); k++)                                            \
                {                                                                                     \
                    if (std::find(Partition[k].begin(), Partition[k].end(), j) != Partition[k].end()) \
                    {   std::cout<<Abc_ObjName(Abc_NtkPi(pNtk, i))<<" "<<Abc_ObjName(Abc_NtkPo(pNtk, j))<<" "<<k<<std::endl;                                                                              \
                        temp.emplace_back(k);                                                         \
                    }                                                                                 \
                }                                                                                     \
            }                                                                                         \
            sign[i] = temp;                                                                           \
        }                                                                                             \
    } while (0);                                                                                      \
    
    SIGN_INIT(pMan->iFuncSupp1, pMan->oPartition1, signI1,pNtk1);
    SIGN_INIT(pMan->iFuncSupp2, pMan->oPartition2, signI2, pNtk2);
    // SIGN_INIT(pMan->oFuncSupp1, pMan->iPartition1, signO1);
    // SIGN_INIT(pMan->oFuncSupp2, pMan->iPartition2, signO2);

    for(int i = 0; i<pMan->iPartition1.size(); i++){
        std::cout<<"{";
        for(auto &j:pMan->iPartition1[i]){
            std::cout<<"{";
            for(auto &k:signI1[j]){
                std::cout<<k<<" ";
            }
            std::cout<<"}";
        }
        std::cout<<"}"<<std::endl;
    }
    for(int i = 0; i<pMan->iPartition2.size(); i++){
        std::cout<<"{";
        for(auto &j:pMan->iPartition2[i]){
            std::cout<<"{";
            for(auto &k:signI2[j]){
                std::cout<<k<<" ";
            }
            std::cout<<"}";
        }
        std::cout<<"}"<<std::endl;
    }
    
}

void Bmatch_MintermPrune(Bmatch_Man_t *pMan, std::map<int, std::vector<int>> PossibleListI, std::map<int, std::vector<int>> PossibleListO, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2)
{
    DdManager *bdd1 = pMan->bdd1;
    DdManager *bdd2 = pMan->bdd2;

    int i;
    Abc_Obj_t *pObj;
    std::map<int, double> MintermTable1, MintermTable2;

    Abc_NtkForEachPo(pNtk1, pObj, i)
    {
        double count = Cudd_CountMinterm(bdd1, (DdNode *)Abc_ObjGlobalBdd(pObj), Abc_NtkPiNum(pNtk1));
        MintermTable1.insert(std::make_pair(i, count));
    }

    Abc_NtkForEachPo(pNtk2, pObj, i)
    {
        double count = Cudd_CountMinterm(bdd2, (DdNode *)Abc_ObjGlobalBdd(pObj), Abc_NtkPiNum(pNtk2));
        MintermTable2.insert(std::make_pair(i, count));
    }
    // for(auto &i)
    // std::cout<<"o1"<<std::endl;
    // for(auto &j:pMan->oPartition1){
    //     std::cout<<"{";
    //     for(auto &k:j){
    //         std::cout<<MintermTable1[k]<<" ";
    //     }
    //     std::cout<<"}"<<std::endl;
    // }
    // std::cout<<"o2"<<std::endl;
    // for(auto &j:pMan->oPartition2){
    //     std::cout<<"{";
    //     for(auto &k:j){
    //         std::cout<<MintermTable2[k]<<" ";
    //     }
    //     std::cout<<"}"<<std::endl;
    // }
}

//case3 test
void Bmatch_MintermWCofactor(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int fVerbose){
    DdManager *bdd1 = pMan->bdd1;
    DdManager *bdd2 = pMan->bdd2;

//     Abc_NtkForEachPo(pNtk1, pObj, i)
//     {
//         if()
//         double count = Cudd_CountMinterm(bdd1, (DdNode *)Abc_ObjGlobalBdd(pObj), Abc_NtkPiNum(pNtk1));
//     }
}