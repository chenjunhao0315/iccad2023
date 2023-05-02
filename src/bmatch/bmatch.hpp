#ifndef __BMATCH_H__
#define __BMATCH_H__

#include <string>
#include <vector>
#include <set>

#include "base/abc/abc.h"

#define VERBOSE_MASK        (1 << 0)
#define VERBOSE_DETAIL_MASK (1 << 1)
#define RESYNTH_MASK        (1 << 2)

#ifdef __cplusplus
extern "C" {
#endif

typedef std::vector<std::pair<int, int> > vSupp;
typedef std::vector<std::pair<std::vector<int>, std::vector<int> > > vGroup;
typedef std::vector<std::vector<std::string> > vsBus;
typedef std::vector<std::vector<int> > vBus;

struct Literal {
    int Var;

    Literal() : Var(-4) {}
    Literal(int Var, int Sign = false) : Var(Var * 2 + (int)Sign) {}
    
    bool     sign()      const { return Var & 1; }
    int      var()       const { return Var >> 1; }
    bool     isConst()   const { return Var < 0; }
    bool     isUndef()   const { return Var == -4; }
};

const Literal CONST1 = Literal(-1, 0);
const Literal CONST0 = Literal(-1, -1);

typedef std::vector<std::vector<Literal> > vMatch;

class Bmatch_Man_t {
public:
    Bmatch_Man_t() {}

public:
    std::string cir1;
    std::string cir2;

    vsBus sBIO1;
    vsBus sBIO2;
    vBus  BI1;
    vBus  BO1;
    vBus  BI2;
    vBus  BO2;

    // pair<PO, suppFunc>
    vSupp suppFunc1;
    vSupp suppFunc2;
};

enum {
    MITER_FAIL,
    RESOURCE_LIMIT,
    PROVE_ERROR,
    EQUIVALENT,
    NON_EQUIVALENT
};

// bmatchMan.cpp
extern Bmatch_Man_t* Bmatch_ManStart();
extern void Bmatch_ManStop(Bmatch_Man_t* p);

// bmatchParser.cpp
extern void Bmatch_ParseInput(Bmatch_Man_t *pMan, char *filename);
extern int Bmatch_ReadNtk(Bmatch_Man_t *pMan, Abc_Ntk_t **ppNtk1, Abc_Ntk_t **ppNtk2);

// bmatchFunc.cpp
extern void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

// bmatchSolve.cpp
extern void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

// bmatchEc.cpp
extern int Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int fVerbose);

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

#ifdef __cplusplus
}
#endif

#endif // __BMATCH_H__