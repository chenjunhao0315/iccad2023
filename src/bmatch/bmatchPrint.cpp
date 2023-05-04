#include "bmatch.hpp"

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_ObjPrint(Abc_Obj_t *pObj);
void Bmatch_NtkPrint(Abc_Ntk_t *pNtk);
void Bmatch_NtkPrintIO(Abc_Ntk_t *pNtk);
void Bmatch_PrintOutputGroup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vGroup &group);
void Bmatch_PrintMatching(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch& MO);
void Bmatch_PrintBusInfo(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_PrintInputSense(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_PrintOutputSupport(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
void Bmatch_PrintSymm(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);

#ifdef __cplusplus
}
#endif

void Bmatch_NtkPrint(Abc_Ntk_t *pNtk) {
    Abc_Print(1, "%s: ", Abc_NtkName(pNtk));
    Abc_Print(1, "i/o: %d/%d ", Abc_NtkPiNum(pNtk), Abc_NtkPoNum(pNtk));
    Abc_Print(1, "nodes: %d ", Abc_NtkNodeNum(pNtk));
    Abc_Print(1, "\n");
}

void Bmatch_NtkPrintIO(Abc_Ntk_t *pNtk) {
    int i;
    Abc_Obj_t *pObj;

    Abc_NtkForEachPi(pNtk, pObj, i) {
        Abc_Print(1, "PI %6s: ", Abc_ObjName(pObj));
        Abc_ObjPrint(stdout, pObj);
    }
    Abc_NtkForEachPo(pNtk, pObj, i) {
        Abc_Print(1, "PO %6s: ", Abc_ObjName(pObj));
        Abc_ObjPrint(stdout, pObj);
    }
}

void Bmatch_ObjPrint(Abc_Obj_t *pObj) {
    Abc_Print(1, "%6s: ", Abc_ObjName(pObj));
    Abc_ObjPrint(stdout, pObj);
}

void Bmatch_PrintOutputGroup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vGroup &group) {
    Abc_Print(1, "Guess Output Grouping:\n");
    for (auto &g : group) {
        Abc_Print(1, "  ([");
        for (int i = 0; i < g.first.size(); ++i) {
            Abc_Print(1, "%s", Abc_ObjName(Abc_NtkPo(pNtk1, g.first[i])));
            if (i != g.first.size() - 1) Abc_Print(1, ", ");
        }
        Abc_Print(1, "], [");
        for (int i = 0; i < g.second.size(); ++i) {
            Abc_Print(1, "%s", Abc_ObjName(Abc_NtkPo(pNtk2, g.second[i])));
            if (i != g.second.size() - 1) Abc_Print(1, ", ");
        }
        Abc_Print(1, "])\n");
    }
}

void Bmatch_PrintMatching(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch& MO) {
    int score = 0;

    Abc_Print(1, "NP3 matching:\n");
    Abc_Print(1, "  MI:");
    for (int i = 0; i < MI.size(); ++i) {
        for (auto &p : MI[i]) {
            if (i != MI.size() - 1) Abc_Print(1, " (%c%s, %s)",p.sign() ? '~' : '\0', Abc_ObjName(Abc_NtkPi(pNtk1, i)), Abc_ObjName(Abc_NtkPi(pNtk2, p.var())));
            else Abc_Print(1, " (%s, %s)", p.sign() ? "CONST0" : "CONST1", Abc_ObjName(Abc_NtkPi(pNtk2, p.var())));
        }
    }
    Abc_Print(1, "\n  MO:");
    for (int i = 0; i < MO.size(); ++i) {
        score += (int)!MO[i].empty() + MO[i].size();
        for (auto &p : MO[i]) {
            Abc_Print(1, " (%c%s, %s)",p.sign() ? '~' : '\0', Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, p.var())));
        }
    }
    Abc_Print(1, "\n  Score: %d\n", score);
}

void Bmatch_PrintBusInfo(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    #define BUS_PRINT(TYPENAME, TYPE, BIO, pNtk)                                    \
    do {                                                                            \
        Abc_Print(1, "  ");                                                         \
        Abc_Print(1, #TYPENAME);                                                    \
        Abc_Print(1, " bus:");                                                      \
        for (int i = 0; i < BIO.size(); ++i) {                                      \
            Abc_Print(1, " (");                                                     \
            for (int j = 0; j < BIO[i].size(); ++j) {                               \
                Abc_Print(1, "%s", Abc_ObjName(Abc_Ntk##TYPE(pNtk, BIO[i][j])));    \
                if (j != BIO[i].size() - 1) Abc_Print(1, ", ");                     \
            }                                                                       \
            Abc_Print(1, ")");                                                      \
        }                                                                           \
        Abc_Print(1, "\n");                                                         \
    } while (0)

    if (pMan->BI1.empty() && pMan->BI2.empty() && pMan->BO1.empty() && pMan->BO2.empty())
        Abc_Print(1, "No Bus Information!\n");
    if (!pMan->BI1.empty() || !pMan->BO1.empty())
        Abc_Print(1, "Cir1:\n");
    if (!pMan->BI1.empty()) BUS_PRINT(input , Pi, pMan->BI1, pNtk1);
    if (!pMan->BO1.empty()) BUS_PRINT(output, Po, pMan->BO1, pNtk1);
    if (!pMan->BI2.empty() || !pMan->BO2.empty())
        Abc_Print(1, "Cir2:\n");
    if (!pMan->BI2.empty()) BUS_PRINT(input , Pi, pMan->BI2, pNtk2);
    if (!pMan->BO2.empty()) BUS_PRINT(output, Po, pMan->BO2, pNtk2);

    #undef BUS_PRINT
}

void Bmatch_PrintInputSense(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int i;
    Abc_Obj_t *pObj;
    
    #define PRINT_SENSE(SI, pNtk)                                     \
    do {                                                              \
        Abc_NtkForEachPi(pNtk, pObj, i) {                             \
            Abc_Print(1, "    %s:", Abc_ObjName(pObj));                 \
            for (auto &p : SI[i]) {                                   \
                Abc_Print(1, " %s", Abc_ObjName(Abc_NtkPo(pNtk, p))); \
            }                                                         \
            Abc_Print(1, "\n");                                       \
        }                                                             \
    } while (0)

    Abc_Print(1, "Input Support\n");
    Abc_Print(1, "  Cir1:\n");
    PRINT_SENSE(pMan->iFuncSupp1, pNtk1);
    Abc_Print(1, "  Cir2:\n");
    PRINT_SENSE(pMan->iFuncSupp2, pNtk2);

    #undef PRINT_SENSE
}

void Bmatch_PrintOutputSupport(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    int i;
    Abc_Obj_t *pObj;
    
    #define PRINT_SENSE(FO, SO, RO, pNtk)                             \
    do {                                                              \
        Abc_NtkForEachPo(pNtk, pObj, i) {                             \
            Abc_Print(1, "    %s: functional (", Abc_ObjName(pObj));  \
            for (auto &p : FO[i]) {                                   \
                Abc_Print(1, " %s", Abc_ObjName(Abc_NtkPi(pNtk, p))); \
            }                                                         \
            Abc_Print(1, " ) structural (");                          \
            for (auto &p : SO[i]) {                                   \
                Abc_Print(1, " %s", Abc_ObjName(Abc_NtkPi(pNtk, p))); \
            }                                                         \
            Abc_Print(1, " ) redundant (");                           \
            for (auto &p : RO[i]) {                                   \
                Abc_Print(1, " %s", Abc_ObjName(Abc_NtkPi(pNtk, p))); \
            }                                                         \
            Abc_Print(1, " )\n");                                     \
        }                                                             \
    } while (0)

    Abc_Print(1, "Output Support\n");
    Abc_Print(1, "  Cir1:\n");
    PRINT_SENSE(pMan->oFuncSupp1, pMan->oStrSupp1, pMan->oRedundSupp1, pNtk1);
    Abc_Print(1, "  Cir2:\n");
    PRINT_SENSE(pMan->oFuncSupp2, pMan->oStrSupp2, pMan->oRedundSupp2, pNtk2);

    #undef PRINT_SENSE
}

void Bmatch_PrintSymm(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    #define PRINT_SYMM(vSymm, pNtk)                                       \
    do {                                                                  \
        for (int i = 0; i < Abc_NtkPoNum(pNtk); ++i) {                    \
            Abc_Print(1, "    %s:", Abc_ObjName(Abc_NtkPo(pNtk, i)));     \
            for (auto &g : vSymm[i]) {                                    \
                Abc_Print(1, " (");                                       \
                for (auto &p : g) {                                       \
                    Abc_Print(1, " %s", Abc_ObjName(Abc_NtkPi(pNtk, p))); \
                }                                                         \
                Abc_Print(1, " )");                                       \
            }                                                             \
            Abc_Print(1, "\n");                                           \
        }                                                                 \
    } while (0)

    Abc_Print(1, "Symmetry information\n");
    Abc_Print(1, "  Cir1:\n");
    PRINT_SYMM(pMan->vSymm1, pNtk1);
    Abc_Print(1, "  Cir2:\n");
    PRINT_SYMM(pMan->vSymm2, pNtk2);

    #undef PRINT_SYMM
}

ABC_NAMESPACE_IMPL_END