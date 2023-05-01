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

    return p;
}

void Bmatch_ManStop(Bmatch_Man_t* p) {
    delete p;
}

ABC_NAMESPACE_IMPL_END