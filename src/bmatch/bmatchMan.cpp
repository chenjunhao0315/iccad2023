#include "bmatch.hpp"

ABC_NAMESPACE_IMPL_START

#ifdef __cplusplus
extern "C" {
#endif

Bmatch_Man_t* Bmatch_ManStart();
void Bmatch_ManStop(Bmatch_Man_t* p);

#ifdef __cplusplus
}
#endif

Bmatch_Man_t* Bmatch_ManStart() {
    Bmatch_Man_t *p = new Bmatch_Man_t;
    p->pInputSolver = NULL;
    p->pOutputSolver = NULL;

    return p;
}

void Bmatch_ManStop(Bmatch_Man_t* p) {
    if (p->pInputSolver) Bmatch_sat_solver_delete(p->pInputSolver);
    if (p->pOutputSolver) Bmatch_sat_solver_delete(p->pOutputSolver);
    delete p;
}

ABC_NAMESPACE_IMPL_END
