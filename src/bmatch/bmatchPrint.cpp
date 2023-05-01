#include "bmatch.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_ObjPrint(Abc_Obj_t *pObj);
void Bmatch_NtkPrint(Abc_Ntk_t *pNtk);
void Bmatch_NtkPrintIO(Abc_Ntk_t *pNtk);
void Bmatch_PrintOutputGroup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vGroup &group);
void Bmatch_PrintMatching(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch& MO);

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
    Abc_Print(1, "Output Grouping:\n");
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
            Abc_Print(1, " (%s, %c%s)", Abc_ObjName(Abc_NtkPi(pNtk1, i)), p.sign() ? '~' : '\0', Abc_ObjName(Abc_NtkPi(pNtk2, p.var())));
        }
    }
    Abc_Print(1, "\n  MO:");
    for (int i = 0; i < MO.size(); ++i) {
        score += (int)!MO[i].empty() + MO[i].size();
        for (auto &p : MO[i]) {
            Abc_Print(1, " (%s, %c%s)", Abc_ObjName(Abc_NtkPo(pNtk1, i)), p.sign() ? '~' : '\0', Abc_ObjName(Abc_NtkPo(pNtk2, p.var())));
        }
    }
    Abc_Print(1, "\n  Score: %d\n", score);
}

ABC_NAMESPACE_IMPL_END