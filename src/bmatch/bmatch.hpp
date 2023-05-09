#ifndef __BMATCH_H__
#define __BMATCH_H__

#include <string>
#include <vector>
#include <set>
#include <tuple>

#include "base/abc/abc.h"
#include "sat/bsat/satSolver.h"

#define VERBOSE_MASK        (1 << 0)
#define VERBOSE_DETAIL_MASK (1 << 1)
#define RESYNTH_MASK        (1 << 2)

#ifdef __cplusplus
extern "C" {
#endif

typedef std::vector<std::pair<std::vector<int>, std::vector<int> > > vGroup;
typedef std::vector<std::vector<std::string> > vsBus;
typedef std::vector<std::vector<int> > vBus;
typedef std::vector<std::vector<std::set<int> > > vSymm;
typedef std::vector<std::pair<int, int> > vSymmPair;
typedef std::vector<std::set<int> > vSupp;
typedef std::vector<std::vector<int>> vEqual;

enum { PO = 0, SUPPFUNC = 1, STRFUNC = 2};
typedef std::vector<std::tuple<int, int, int> > vSuppInfo;

struct Literal {
    int Var;

    Literal() : Var(-1) {}
    Literal(int Var, bool Sign = false) : Var(Var * 2 + (int)Sign) {}
    
    bool     sign()      const { return Var & 1; }
    int      var()       const { return Var >> 1; }
    bool     isUndef()   const { return Var == -1; }
    operator int()       const { return Var; }
};

typedef std::vector<std::vector<Literal> > vMatch;
typedef std::vector<std::vector<int>>  vMatch_Group;

class Bmatch_Man_t {
public:
    Bmatch_Man_t() {}

public:
    std::string cir1;
    std::string cir2;

    // bus information
    vsBus sBIO1;
    vsBus sBIO2;
    vBus  BI1, BO1;
    vBus  BI2, BO2;

    // functional support information
    vSupp iFuncSupp1, oFuncSupp1;
    vSupp iFuncSupp2, oFuncSupp2;

    //group information
    vGroup Groups;

    //equal group
    vEqual oEqual1;
    vEqual oEqual2;

    // structrual support information
    vSupp oStrSupp1;
    vSupp oStrSupp2;

    // redundant support information
    vSupp oRedundSupp1;
    vSupp oRedundSupp2;
    std::set<int> sRedund1;
    std::set<int> sRedund2;

    // symmetry group
    vSymm vSymm1;
    vSymm vSymm2;
    vSymmPair vSymmPair1;
    vSymmPair vSymmPair2;

    // functional support information
    // tuple<PO, suppFunc, strFunc>
    vSuppInfo vSuppInfo1;
    vSuppInfo vSuppInfo2;

    // input solver, output solver
    sat_solver *pInputSolver;
    sat_solver *pOutputSolver;
    int ni, mi;
    int no, mo;

    //learned clause level
    std::vector<int> LearnedLevel;
    std::vector<int> LearnedAssumption;
    int Projective;
    bool AllowProjection;
    std::vector<int> ClauseControl;
    vMatch_Group MO;
};

enum {
    MITER_FAIL,
    RESOURCE_LIMIT,
    PROVE_ERROR,
    EQUIVALENT,
    NON_EQUIVALENT,
    UNDEF
};

struct EcResult {
    EcResult() : status(UNDEF), model(NULL) {}
    EcResult(int status, int *model) : status(status), model(model) {}

    int status;
    int *model;
};

struct InputMapping {
    InputMapping() : status(0) {}
    InputMapping(int status, vMatch MI) : status(status), MI(MI) {}

    int status;
    vMatch MI;
};

// bmatchMan.cpp
extern Bmatch_Man_t* Bmatch_ManStart();
extern void Bmatch_ManStop(Bmatch_Man_t* p);

// bmatchParser.cpp
extern void Bmatch_ParseInput(Bmatch_Man_t *pMan, char *filename);
extern int Bmatch_ReadNtk(Bmatch_Man_t *pMan, Abc_Ntk_t **ppNtk1, Abc_Ntk_t **ppNtk2);

// bmatchFunc.cpp
extern void Bmatch_CalCir1Redund(Abc_Ntk_t *pNtk1, vSupp &oStrSupp, std::set<int> &sRedund);
extern void Bmatch_CalCir2RedundWithGivenMapping(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, std::set<int> &sRedund);
extern void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

// bmatchSolve.cpp
extern void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

// bmatchEc.cpp
extern EcResult Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int fVerbose);

// bmatchMiter.cpp
extern Abc_Ntk_t* Bmatch_NtkMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);

// bmatchSetup.cpp
extern void Bmatch_NtkSetup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);
extern Abc_Ntk_t* Bmatch_NtkResynth(Abc_Ntk_t *pNtk);

// bmatchPrint.cpp
extern void Bmatch_ObjPrint(Abc_Obj_t *pObj);
extern void Bmatch_NtkPrint(Abc_Ntk_t *pNtk);
extern void Bmatch_NtkPrintIO(Abc_Ntk_t *pNtk);
extern void Bmatch_PrintOutputGroup(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vGroup &group);
extern void Bmatch_PrintMatching(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch& MO);
extern void Bmatch_PrintBusInfo(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern void Bmatch_PrintInputSupport(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern void Bmatch_PrintOutputSupport(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern void Bmatch_PrintSymm(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern void Bmatch_PrintEqual(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);

#ifdef __cplusplus
}
#endif

#endif // __BMATCH_H__