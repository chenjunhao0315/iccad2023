#include "bmatch.hpp"
#include "base/main/mainInt.h"

#include "print.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

void Bmatch_ParseInput(Bmatch_Man_t *pMan, char *filename);
int Bmatch_ReadNtk(Bmatch_Man_t *pMan, Abc_Ntk_t **ppNtk1, Abc_Ntk_t **ppNtk2);

#ifdef __cplusplus
}
#endif

void Bmatch_ParseInput(Bmatch_Man_t *pMan, char *filename) {
    FILE *pFile;
    char buffer[100];
    int i, j, nBuses, nPorts;

    pFile = fopen(filename, "r");
    if (!pFile) {
        Abc_Print(-1, "Cannot open the file \"%s\"", filename);
    }

    fscanf(pFile, "%s\n", buffer);
    pMan->cir1 = std::string(buffer);

    fscanf(pFile, "%d\n", &nBuses);
    for (i = 0; i < nBuses; ++i) {
        std::vector<std::string> bus;
        fscanf(pFile, " %d", &nPorts);
        for (j = 0; j < nPorts; ++j) {
            fscanf(pFile, " %s", buffer);
            bus.emplace_back(buffer);
        }
        pMan->sBIO1.emplace_back(std::move(bus));
    }

    fscanf(pFile, "%s\n", buffer);
    pMan->cir2 = std::string(buffer);

    fscanf(pFile, "%d\n", &nBuses);
    for (i = 0; i < nBuses; ++i) {
        std::vector<std::string> bus;
        fscanf(pFile, " %d", &nPorts);
        for (j = 0; j < nPorts; ++j) {
            fscanf(pFile, " %s", buffer);
            bus.emplace_back(buffer);
        }
        pMan->sBIO2.emplace_back(std::move(bus));
    }

    fclose(pFile);
}

int Bmatch_ReadNtk(Bmatch_Man_t *pMan, Abc_Ntk_t **ppNtk1, Abc_Ntk_t **ppNtk2) {
    Abc_Ntk_t *pNtk1, *pNtk2, *pNtkTemp;
    
    pNtk1 = Io_Read((char*)pMan->cir1.c_str(), Io_ReadFileType((char*)pMan->cir1.c_str()), 1, 0);
    if (!pNtk1) return 0;
    pNtk2 = Io_Read((char*)pMan->cir2.c_str(), Io_ReadFileType((char*)pMan->cir2.c_str()), 1, 0);
    if (!pNtk2) {
        Abc_NtkDelete(pNtk1);
        return 0;
    }

    if (!Abc_NtkIsStrash(pNtk1)) {
        pNtkTemp = Abc_NtkStrash(pNtk1, 0, 1, 0);
        Abc_NtkDelete(pNtk1);
        pNtk1 = pNtkTemp;
    }
    if (!Abc_NtkIsStrash(pNtk2)) {
        pNtkTemp = Abc_NtkStrash(pNtk2, 0, 1, 0);
        Abc_NtkDelete(pNtk2);
        pNtk2 = pNtkTemp;
    }

    *ppNtk1 = pNtk1;
    *ppNtk2 = pNtk2;

    return 1;
}

ABC_NAMESPACE_IMPL_END