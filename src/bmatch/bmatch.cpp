#include "bmatch.hpp"
#include "base/main/mainInt.h"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_Init(Abc_Frame_t *pAbc);
void Bmatch_End(Abc_Frame_t *pAbc);

static int Abc_CommandBmatch(Abc_Frame_t *pAbc, int argc, char *argv[]);
static int Abc_CommandBmatchGroup(Abc_Frame_t *pAbc, int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

void Bmatch_Init(Abc_Frame_t *pAbc) {
    Cmd_CommandAdd(pAbc, "z Bmatch", "bmatch", Abc_CommandBmatch, 0);
    Cmd_CommandAdd(pAbc, "z BmatchGroup", "bmatchgroup", Abc_CommandBmatchGroup, 0);
}

static int Abc_CommandBmatch(Abc_Frame_t *pAbc, int argc, char **argv) {
    int c = 0;
    int option = 0;
    int nArgc;
    char **nArgv;

    int pDelete1 = 0, pDelete2 = 0;
    Abc_Ntk_t *pNtk1, *pNtk2;
    Bmatch_Man_t *pMan = Bmatch_ManStart();
    
    Extra_UtilGetoptReset();
    while ((c = Extra_UtilGetopt(argc, argv, "vdrh")) != EOF) {
        switch (c) {
            case 'v': option ^= VERBOSE_MASK; break;
            case 'd': option ^= (VERBOSE_MASK | VERBOSE_DETAIL_MASK); break;
            case 'r': option ^= RESYNTH_MASK; break;
            case 'h': goto usage;
            default : goto usage;
        }
    }

    nArgc = argc - globalUtilOptind;
    nArgv = argv + globalUtilOptind;

    // reserve for reading group information

    if (nArgc != 3) {
        Abc_Print(-1, "Invalid command!\n");
        goto usage;
    }

    pMan->cir1 = std::string(nArgv[0]);
    pMan->cir2 = std::string(nArgv[1]);
    Bmatch_ReadNtk(pMan, &pNtk1, &pNtk2);
    Bmatch_Preprocess(pMan, pNtk1, pNtk2, option);
    Bmatch_SolveNP3(pMan, pNtk1, pNtk2, option, nArgv[2]);

    Bmatch_ManStop(pMan);

    return 0;
usage:
    Abc_Print(-2, "usage: bmatch <cir1> <cir2>\n");
    Abc_Print(-2, "\t-v         : verbosity [default = %d]\n", option & VERBOSE_MASK);
    Abc_Print(-2, "\t-v         : detail verbosity [default = %d]\n", option & VERBOSE_DETAIL_MASK);
    Abc_Print(-2, "\t-r         : resynthesize [default = %d]\n", option & RESYNTH_MASK);
    Abc_Print(-2, "\t-h         : print the command usage\n");

    return 1;
}

static int Abc_CommandBmatchGroup(Abc_Frame_t *pAbc, int argc, char *argv[]) {
    int c = 0;
    int option = 0;
    int nArgc;
    char **nArgv;

    Abc_Ntk_t *pNtk1, *pNtk2;
    Bmatch_Man_t *pMan = Bmatch_ManStart();
    
    Extra_UtilGetoptReset();
    while ((c = Extra_UtilGetopt(argc, argv, "vdrh")) != EOF) {
        switch (c) {
            case 'v': option ^= VERBOSE_MASK; break;
            case 'd': option ^= (VERBOSE_MASK | VERBOSE_DETAIL_MASK); break;
            case 'r': option ^= RESYNTH_MASK; break;
            case 'h': goto usage;
            default : goto usage;
        }
    }

    nArgc = argc - globalUtilOptind;
    nArgv = argv + globalUtilOptind;

    // reserve for reading group information
    if (nArgc != 2) {
        Abc_Print(-1, "Invalid command!\n");
        goto usage;
    }

    Bmatch_ParseInput(pMan, nArgv[0]);
    Bmatch_ReadNtk(pMan, &pNtk1, &pNtk2);
    Bmatch_Preprocess(pMan, pNtk1, pNtk2, option);

    Bmatch_SolveNP3(pMan, pNtk1, pNtk2, option, nArgv[1]);

    Bmatch_ManStop(pMan);

    return 0;
usage:
    Abc_Print(-2, "usage: bmatchgroup <def file>\n");
    Abc_Print(-2, "\t-v         : verbosity [default = %d]\n", option & VERBOSE_MASK);
    Abc_Print(-2, "\t-v         : detail verbosity [default = %d]\n", option & VERBOSE_DETAIL_MASK);
    Abc_Print(-2, "\t-r         : resynthesize [default = %d]\n", option & RESYNTH_MASK);
    Abc_Print(-2, "\t-h         : print the command usage\n");

    return 1;
}

ABC_NAMESPACE_IMPL_END