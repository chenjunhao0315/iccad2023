#ifndef __BMATCH_H__
#define __BMATCH_H__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <tuple>
#include <array>
#include <cmath>
#include <algorithm>

#include "base/abc/abc.h"
#include "AutoBuffer.hpp"
#include "bmatchCadicalSat.hpp"

#define VERBOSE_MASK        (1 << 0)
#define VERBOSE_DETAIL_MASK (1 << 1)
#define RESYNTH_MASK        (1 << 2)

#ifdef __cplusplus
extern "C" {
#endif

#define S_ONE (~0)
#define S_ZERO (0)

typedef std::vector<std::pair<std::vector<int>, std::vector<int> > > vGroup;
typedef std::vector<std::vector<std::string> > vsBus;
typedef std::vector<std::vector<int> > vBus;
typedef std::vector<std::vector<std::set<int> > > vSymm;
typedef std::vector<std::pair<int, int> > vSymmPair;
typedef std::vector<std::set<int> > vSupp;
typedef std::vector<std::vector<int> > Mat;
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

struct Prob {
    float data;
    int id;
};

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
    AutoBuffer<int> sRedund1;
    AutoBuffer<int> sRedund2;

    // symmetry group
    vSymm vSymm1;
    vSymm vSymm2;
    vSymmPair vSymmPair1;
    vSymmPair vSymmPair2;

    // functional support information
    // tuple<PO, suppFunc, strFunc>
    vSuppInfo vSuppInfo1;
    vSuppInfo vSuppInfo2;

    // unateness
    Mat unateMat1;
    Mat unateMat2;

    // qbf input solver
    int nControlPi;
    AutoBuffer<int> possibleMI;
    Mat MapReduceMI;

    // input solver, output solver
    CaDiCaL::Solver *pInputSolver;
    CaDiCaL::Solver *pOutputSolver;
    int ni, mi;
    int no, mo;
    AutoBuffer<int> impossibleMI;

    // controllable miter
    CaDiCaL::Solver *pMiterSolver;
    AutoBuffer<int> inputControl;
    int controlOffset;
    CaDiCaL::Solver *pMiterSolverNew;

    //learned clause level
    std::vector<int> LearnedLevel;
    std::vector<int> LearnedAssumption;
    int Projective;
    bool AllowProjection;
    std::vector<int> ClauseControl;
    vMatch_Group MO;

    AutoBuffer<Prob> prob1;
    AutoBuffer<Prob> prob2;
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

// bmatchUnate.cpp
extern void Bmatch_RandomSimUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int iter);
extern void Bmatch_SatUnate(Abc_Ntk_t *pNtk, Mat &unateMat, int Po, int iter);
extern void Bmatch_GetUnateCount(Mat &unateMat, int *Binate, int *Unate);

// bmatchFunc.cpp
extern void Bmatch_CalCir1Redund(Abc_Ntk_t *pNtk1, vSupp &oStrSupp, std::set<int> &sRedund);
extern void Bmatch_CalCir2RedundWithGivenMapping(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, std::set<int> &sRedund);
extern void Bmatch_Preprocess(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option);

// bmatchSolve.cpp
extern void Bmatch_SolveNP3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int option, char *output);

// bmatchInputSolver.cpp
extern void Bmatch_InitControllableInputMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern void Bmatch_InitControllableInputOutputMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern int Bmatch_InitInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern int Bmatch_ApplyInputSolverRowConstraint(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern int Bmatch_PruneInputSolverByBusOrdered(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern int Bmatch_PruneInputSolverByBusExactMap(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern int Bmatch_PruneInputSolverByStrSupport(Bmatch_Man_t *pMan, vMatch &MO);
extern int Bmatch_PruneInputSolverByFuncSupport(Bmatch_Man_t *pMan, vMatch &MO);
extern int Bmatch_PruneInputSolverBySymmetryProperty(Bmatch_Man_t *pMan, vMatch &MO);
extern int Bmatch_PruneInputSolverByUnate(Bmatch_Man_t *pMan, vMatch &MO);
extern int Bmatch_PruneInputSolverBySymmetry(Bmatch_Man_t *pMan, vMatch &MI);
extern int Bmatch_PruneInputSolverByCounterPart(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *pModel1, vMatch& MI, vMatch& MO);
extern InputMapping Bmatch_HeuristicSolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern InputMapping Bmatch_SolveInput(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int *bLits, int *eLits, int fVerbose);
extern InputMapping Bmatch_SolveInputQbf(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

// bmatchInputSolver2.cpp
extern void Bmatch_InitQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern void Bmatch_FillPossibleMIbyStrSupp(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern void Bmatch_ReducePossibleMIbyUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern Abc_Ntk_t *Bmatch_NtkQbfMiter(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern InputMapping Bmatch_SolveQbfInputSolver(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern InputMapping Bmatch_SolveQbfInputSolver2(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern InputMapping Bmatch_SolveQbfInputSolver3(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);
extern InputMapping Bmatch_SolveQbfInputSolver4(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO);

// bmatchEc.cpp
extern EcResult Bmatch_NtkEcFraig(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO, int cadicalSat, int fVerbose);
extern EcResult Bmatch_NtkControllableInputEcFraig(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI);
extern EcResult Bmatch_NtkControllableInputOutputEcFraig(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
extern EcResult Bmatch_NtkEcGia(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
extern CaDiCaL::Solver *Bmatch_ControllableInputOutputSat(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, int &controlPiOffset);
extern CaDiCaL::Solver *Bmatch_ControllableInputSat(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, int &controlPiOffset);

// bmatchMiter.cpp
extern Abc_Ntk_t *Bmatch_NtkMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
extern Abc_Ntk_t *Bmatch_NtkControllableInputMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MO, int inv = 0);
extern Abc_Ntk_t *Bmatch_NtkControllableInputOutputMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern Abc_Ntk_t *Bmatch_NtkGiaMiter(Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2, vMatch &MI, vMatch &MO);
extern void Bmatch_NtkMiterAddOne(Abc_Ntk_t *pNtk, Abc_Ntk_t *pNtkMiter);
extern Abc_Obj_t *Bmatch_NtkCreateOr(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal);
extern Abc_Obj_t *Bmatch_NtkCreateAnd(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pSignal);
extern Abc_Obj_t *Bmatch_NtkCreateParallelCase(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal);
extern Abc_Obj_t *Bmatch_NtkCreateMultiplexer(Abc_Aig_t *pMan, std::vector<Abc_Obj_t *> &pControl, std::vector<Abc_Obj_t *> &pSignal, int dontApplyNot);
extern void Bmatch_NtkBuildWithCone(Abc_Ntk_t *pNtk, Abc_Ntk_t *pNtkNew, std::vector<int> &outs);

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
extern void Bmatch_PrintUnate(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);
extern void Bmatch_PrintEqual(Bmatch_Man_t *pMan, Abc_Ntk_t *pNtk1, Abc_Ntk_t *pNtk2);

#ifdef __cplusplus
}
#endif

#endif // __BMATCH_H__
