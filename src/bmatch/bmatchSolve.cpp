#include "bmatch.hpp"

#include <tuple>

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

vGroup Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan);
std::tuple<vMatch, vMatch> Bmatch_SolveInputOutputMatch(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);

#ifdef __cplusplus
}
#endif

void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option) {
    int Status;

    if (option & VERBOSE_MASK) Bmatch_PrintBusInfo(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintInputSense(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintOutputSupport(pMan, pNtk1, pNtk2);
    if (option & VERBOSE_MASK) Bmatch_PrintSymm(pMan, pNtk1, pNtk2);

    // testing flow
    auto groups = Bmatch_SolveOutputGroup(pMan);

    if (option & VERBOSE_MASK) Bmatch_PrintOutputGroup(pNtk1, pNtk2, groups);

    vMatch MI, MO;
    std::tie(MI, MO) = Bmatch_SolveInputOutputMatch(pMan, pNtk1, pNtk2);

    if (option & VERBOSE_MASK) Bmatch_PrintMatching(pNtk1, pNtk2, MI, MO);

    assert(MI.size() == Abc_NtkPiNum(pNtk1) + 1);
    assert(MO.size() == Abc_NtkPoNum(pNtk1));
    Status = Bmatch_NtkEcFraig(pNtk1, pNtk2, MI, MO, 0);

    print("Status:", (Status == EQUIVALENT) ? "EQUIVALENT" : "NON-EQUIVALENT");
}

std::tuple<vMatch, vMatch> Bmatch_SolveInputOutputMatch(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2) {
    vMatch MI = {{Literal(1, false), Literal(3, false)}, {Literal(2, true)}, {Literal(0, true)}, {}};
    vMatch MO = {{Literal(0, false)}};

    // vMatch MI = {{Literal(2, false)}, {Literal(3, false)}, {Literal(0, false)}, {Literal(1, false)}, {Literal(4, false)}, {}};
    // vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {}, {}};

    // vMatch MI = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}, {Literal(4, false)}, {}};
    // vMatch MO = {{Literal(0, false)}, {Literal(1, false)}, {Literal(2, false)}, {Literal(3, false)}};

    // vMatch MI(Abc_NtkPiNum(pNtk1) + 1, std::vector<Literal>());
    // vMatch MO(Abc_NtkPoNum(pNtk1), std::vector<Literal>());

    return std::make_tuple(MI, MO);
}

vGroup Bmatch_SolveOutputGroup(Bmatch_Man_t *pMan) {
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

    return groups;
}

ABC_NAMESPACE_IMPL_END