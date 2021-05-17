/**CFile****************************************************************

  FileName    [giaMinLut.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Scalable AIG package.]

  Synopsis    [Collapsing AIG.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - June 20, 2005.]

  Revision    [$Id: giaMinLut.c,v 1.00 2005/06/20 00:00:00 alanmi Exp $]

***********************************************************************/

#include "gia.h"
#include "giaAig.h"
#include "base/main/mainInt.h"
#include "opt/sfm/sfm.h"

#ifdef ABC_USE_CUDD
#include "bdd/extrab/extraBdd.h"
#include "bdd/dsd/dsd.h"
#endif

ABC_NAMESPACE_IMPL_START

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

extern Abc_Ntk_t * Abc_NtkFromAigPhase( Aig_Man_t * pMan );

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Vec_WrdReadText( char * pFileName, Vec_Wrd_t ** pvSimI, Vec_Wrd_t ** pvSimO, int nIns, int nOuts )
{
    int i, nSize, iLine, nLines, nWords;
    char pLine[1000];
    Vec_Wrd_t * vSimI, * vSimO;
    FILE * pFile = fopen( pFileName, "rb" );
    if ( pFile == NULL )
    {
        printf( "Cannot open file \"%s\" for reading.\n", pFileName );
        return;
    }
    fseek( pFile, 0, SEEK_END );  
    nSize = ftell( pFile );  
    if ( nSize % (nIns + nOuts + 1) > 0 )
    {
        printf( "Cannot read file with simulation data that is not aligned at 8 bytes (remainder = %d).\n", nSize % (nIns + nOuts + 1) );
        fclose( pFile );
        return;
    }
    rewind( pFile );
    nLines = nSize / (nIns + nOuts + 1);
    nWords = (nLines + 63)/64;
    vSimI  = Vec_WrdStart( nIns *nWords );
    vSimO  = Vec_WrdStart( nOuts*nWords );
    for ( iLine = 0; fgets( pLine, 1000, pFile ); iLine++ )
    {
        for ( i = 0; i < nIns; i++ )
            if ( pLine[nIns-1-i] == '1' )
                Abc_TtXorBit( Vec_WrdArray(vSimI) + i*nWords, iLine );
            else assert( pLine[nIns-1-i] == '0' );
        for ( i = 0; i < nOuts; i++ )
            if ( pLine[nIns+nOuts-1-i] == '1' )
                Abc_TtXorBit( Vec_WrdArray(vSimO) + i*nWords, iLine );
            else assert( pLine[nIns+nOuts-1-i] == '0' );
    }
    fclose( pFile );
    *pvSimI = vSimI;
    *pvSimO = vSimO;
    printf( "Read %d words of simulation data for %d inputs and %d outputs (padded %d zero-patterns).\n", nWords, nIns, nOuts, nWords*64-nLines );
}
void Gia_ManSimInfoTransform( char * pFileName, char * pFileOut1, char * pFileOut2, int nIns, int nOuts )
{
    Vec_Wrd_t * vSimI, * vSimO;
    Vec_WrdReadText( pFileName, &vSimI, &vSimO, nIns, nOuts );
    Vec_WrdDumpBin( pFileOut1 ? pFileOut1 : Extra_FileNameGenericAppend(pFileName, ".simi"), vSimI, 1 );
    Vec_WrdDumpBin( pFileOut2 ? pFileOut2 : Extra_FileNameGenericAppend(pFileName, ".simo"), vSimO, 1 );
    Vec_WrdFree( vSimI );
    Vec_WrdFree( vSimO );
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
Vec_Wrd_t * Vec_WrdZoneExtract( int ZoneSize, Vec_Wrd_t * p, int iWord, int nWords )
{
    int z, nZones = Vec_WrdSize(p)/ZoneSize;
    int w, Limit = Abc_MinInt( nWords, ZoneSize-iWord );
    Vec_Wrd_t * pNew = Vec_WrdStart( nZones*nWords );
    for ( z = 0; z < nZones; z++ )
        for ( w = 0; w < Limit; w++ )
            Vec_WrdWriteEntry( pNew, z*nWords + w, Vec_WrdEntry(p, z*ZoneSize + iWord + w) );
    return pNew;
}
void Vec_WrdZoneInsert( Vec_Wrd_t * pNew, int ZoneSize, Vec_Wrd_t * p, int iWord, int nWords )
{
    int z, nZones = Vec_WrdSize(pNew)/ZoneSize;
    int w, Limit = Abc_MinInt( nWords, ZoneSize-iWord );
    for ( z = 0; z < nZones; z++ )
        for ( w = 0; w < Limit; w++ )
            Vec_WrdWriteEntry( pNew, z*ZoneSize + iWord + w, Vec_WrdEntry(p, z*nWords + w) );
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Gia_ManSimInfoPrintOne( Gia_Man_t * p, Vec_Wrd_t * vSimsIn, Vec_Wrd_t * vSimsOut, int nWords, int nPats )
{
    int Id, i, k;
    for ( k = 0; k < nPats; k++ )
    {
        Gia_ManForEachCiId( p, Id, i )
    //        printf( "%d", Vec_WrdEntry(p->vSims, p->nSimWords*Id) & 1 );
            printf( "%d", (int)(Vec_WrdEntry(vSimsIn,  nWords*i) >> k) & 1 );
        printf( " " );
        Gia_ManForEachCoId( p, Id, i )
    //        printf( "%d", Vec_WrdEntry(p->vSims, p->nSimWords*Id) & 1 );
            printf( "%d", (int)(Vec_WrdEntry(vSimsOut, nWords*i) >> k) & 1 );
        printf( "\n" );
    }
}
Vec_Wrd_t * Gia_ManSimInfoTryOne( Gia_Man_t * p, Vec_Wrd_t * vSimI, int fPrint )
{
    extern Vec_Wrd_t * Gia_ManSimulateWordsOut( Gia_Man_t * p, Vec_Wrd_t * vSimsIn );
    Vec_Wrd_t * vSimsOut = Gia_ManSimulateWordsOut( p, vSimI );
    int nWords = Vec_WrdSize(vSimI) / Gia_ManCiNum(p);
    assert( Vec_WrdSize(vSimI) % Gia_ManCiNum(p) == 0 );
    if ( fPrint )
        Gia_ManSimInfoPrintOne( p, vSimI, vSimsOut, nWords, 6 );
    return vSimsOut;
}
int Gia_ManSimEvalOne( Gia_Man_t * p, Vec_Wrd_t * vSimO, Vec_Wrd_t * vSimO_new )
{
    int i, Count = 0, nWords = Vec_WrdSize(vSimO) / Gia_ManCoNum(p);
    word * pSim0 = ABC_CALLOC( word, nWords );
    assert( Vec_WrdSize(vSimO) == Vec_WrdSize(vSimO_new) );
    for ( i = 0; i < Gia_ManCoNum(p); i++ )
    {
        word * pSimGold = Vec_WrdEntryP( vSimO,     i * nWords );
        word * pSimImpl = Vec_WrdEntryP( vSimO_new, i * nWords );
        Abc_TtOrXor( pSim0, pSimImpl, pSimGold, nWords );
    }
    Count = Abc_TtCountOnesVec( pSim0, nWords );
    printf( "Number of failed patterns is %d (%8.4f %% of %d). The first one is %d.\n", 
        Count, 100.0*Count/(64*nWords), 64*nWords, Abc_TtFindFirstBit2(pSim0, nWords) );
    ABC_FREE( pSim0 );
    return Count;
}
int Gia_ManSimEvalOne2( Gia_Man_t * p, Vec_Wrd_t * vSimO, Vec_Wrd_t * vSimO_new )
{
    int i, Count = 0, nWords = Vec_WrdSize(vSimO) / Gia_ManCoNum(p);
    word * pSim0 = ABC_CALLOC( word, nWords );
    assert( Vec_WrdSize(vSimO) == Vec_WrdSize(vSimO_new) );
    for ( i = 0; i < Gia_ManCoNum(p); i++ )
    {
        word * pSimGold = Vec_WrdEntryP( vSimO,     i * nWords );
        word * pSimImpl = Vec_WrdEntryP( vSimO_new, i * nWords );
        Abc_TtXor( pSim0, pSimImpl, pSimGold, nWords, 0 );
        Count += Abc_TtCountOnesVec( pSim0, nWords );
    }
    printf( "Number of failed patterns is %d (%8.4f %% of %d). The first one is %d.\n", 
        Count, 100.0*Count/(64*nWords*Gia_ManCoNum(p)), 64*nWords*Gia_ManCoNum(p), Abc_TtFindFirstBit2(pSim0, nWords) );
    ABC_FREE( pSim0 );
    return Count;
}
Vec_Wrd_t * Gia_ManSimInfoTry( Gia_Man_t * p, Vec_Wrd_t * vSimI )
{
    int nWords = Vec_WrdSize(vSimI) / Gia_ManCiNum(p);
    int w, nWordsOne = 200, nWordBatches = (nWords + nWordsOne - 1)/nWordsOne;
    Vec_Wrd_t * vSimO_new = Vec_WrdStart( nWords * Gia_ManCoNum(p) );
    for ( w = 0; w < nWordBatches; w++ )
    {
        //int Value = printf( "%3d / %3d : ", w, nWordBatches );
        Vec_Wrd_t * vSimI_ = Vec_WrdZoneExtract( nWords, vSimI, w*nWordsOne, nWordsOne );
        Vec_Wrd_t * vSimO_ = Gia_ManSimInfoTryOne( p, vSimI_, 0 );
        Vec_WrdZoneInsert( vSimO_new, nWords, vSimO_, w*nWordsOne, nWordsOne );
        Vec_WrdFree( vSimI_ );
        Vec_WrdFree( vSimO_ );
        //Value = 0;
    }
    return vSimO_new;
}
int Gia_ManSimInfoEval( Gia_Man_t * p, Vec_Wrd_t * vSimO, Vec_Wrd_t * vSimO_new )
{
    int nResult = Gia_ManSimEvalOne2(p, vSimO, vSimO_new);
    //Vec_WrdDumpBin( "temp.simo", vSimO_new, 1 );
    printf( "Total errors = %d.  ", nResult );
    printf( "Density of output patterns %8.4f.\n", (float)Abc_TtCountOnesVec(Vec_WrdArray(vSimO_new), Vec_WrdSize(vSimO_new))/(64*Vec_WrdSize(vSimO_new)) );
    return nResult;
}
void Gia_ManSimInfoPassTest( Gia_Man_t * p, char * pFileName, char * pFileName2, int fCompare )
{
    abctime clk = Abc_Clock();
    if ( fCompare )
    {
        Vec_Wrd_t * vSim1 = Vec_WrdReadBin( pFileName,  1 );
        Vec_Wrd_t * vSim2 = Vec_WrdReadBin( pFileName2, 1 );
        printf( "Density of input  patterns %8.4f.\n", (float)Abc_TtCountOnesVec(Vec_WrdArray(vSim1), Vec_WrdSize(vSim1))/(64*Vec_WrdSize(vSim1)) );
        printf( "Density of output patterns %8.4f.\n", (float)Abc_TtCountOnesVec(Vec_WrdArray(vSim2), Vec_WrdSize(vSim2))/(64*Vec_WrdSize(vSim2)) );
        Gia_ManSimInfoEval( p, vSim1, vSim2 );
        Vec_WrdFree( vSim1 );
        Vec_WrdFree( vSim2 );
    }
    else
    {
        Vec_Wrd_t * vSimI = Vec_WrdReadBin( pFileName, 1 );
        Vec_Wrd_t * vSimO = Gia_ManSimInfoTry( p, vSimI );
        printf( "Density of input  patterns %8.4f.\n", (float)Abc_TtCountOnesVec(Vec_WrdArray(vSimI), Vec_WrdSize(vSimI))/(64*Vec_WrdSize(vSimI)) );
        printf( "Density of output patterns %8.4f.\n", (float)Abc_TtCountOnesVec(Vec_WrdArray(vSimO), Vec_WrdSize(vSimO))/(64*Vec_WrdSize(vSimO)) );
        Vec_WrdDumpBin( pFileName2, vSimO, 1 );
        Vec_WrdFree( vSimI );
        Vec_WrdFree( vSimO );
    }
    Abc_PrintTime( 1, "Time", Abc_Clock() - clk );
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
word * Gia_ManCountFraction( Gia_Man_t * p, Vec_Wrd_t * vSimI, Vec_Int_t * vSupp, int Thresh, int fVerbose, int * pCare )
{
    Gia_Obj_t * pObj;
    int i, k, nUsed = 0, nGood = 0;
    int nWords = Vec_WrdSize(vSimI) / Gia_ManCiNum(p);
    int nMints = 1 << Vec_IntSize(vSupp);
    word ** pSims = ABC_ALLOC( word *, Vec_IntSize(vSupp) );
    word * pRes   = ABC_CALLOC( word, Abc_Truth6WordNum(Vec_IntSize(vSupp)) );
    int * pCounts = ABC_CALLOC( int, nMints );
    Gia_ManForEachObjVec( vSupp, p, pObj, i )
        pSims[i] = Vec_WrdEntryP( vSimI, Gia_ObjCioId(pObj) * nWords );
    for ( k = 0; k < 64*nWords; k++ )
    {
        int iMint = 0;
        for ( i = 0; i < Vec_IntSize(vSupp); i++ )
            if ( Abc_TtGetBit(pSims[i], k) )
                iMint |= 1 << i;
        assert( iMint < nMints );
        pCounts[iMint]++;
    }
    for ( k = 0; k < nMints; k++ )
    {
        nUsed += (pCounts[k] > 0);
        nGood += (pCounts[k] >= Thresh);
        if ( pCounts[k] >= Thresh )
            Abc_TtXorBit( pRes, k );
        //printf( "%d ", pCounts[k] );
    }
    //printf( "\n" );
    if ( fVerbose )
    printf( "Used %4d and good %4d (out of %4d).\n", nUsed, nGood, nMints ); 
    ABC_FREE( pSims );
    ABC_FREE( pCounts );
    *pCare = nGood;
    return pRes;
}
void Gia_ManCollectSupp_rec( Gia_Man_t * p, int iObj, Vec_Int_t * vSupp )
{
    Gia_Obj_t * pObj;
    if ( !iObj || Gia_ObjIsTravIdCurrentId(p, iObj) )
        return;
    Gia_ObjSetTravIdCurrentId(p, iObj);
    pObj = Gia_ManObj( p, iObj );
    if ( Gia_ObjIsCi(pObj) )
    {
        //Vec_IntPush( vSupp, Gia_ObjCioId(pObj) );
        Vec_IntPush( vSupp, Gia_ObjId(p, pObj) );
        return;
    }
    assert( Gia_ObjIsAnd(pObj) );
    Gia_ManCollectSupp_rec( p, Gia_ObjFaninId0(pObj, iObj), vSupp );
    Gia_ManCollectSupp_rec( p, Gia_ObjFaninId1(pObj, iObj), vSupp );
}
Vec_Int_t * Gia_ManCollectSupp( Gia_Man_t * p, int iOut, int nOuts )
{
    Vec_Int_t * vSupp = Vec_IntAlloc( 16 ); int i;
    Gia_ManIncrementTravId( p );
    for ( i = 0; i < nOuts; i++ )
        Gia_ManCollectSupp_rec( p, Gia_ObjFaninId0p(p, Gia_ManCo(p, iOut+i)), vSupp );
    return vSupp;
}
Gia_Man_t * Gia_ManPerformLNetOpt( Gia_Man_t * p, char * pFileName, int nIns, int nOuts, int Thresh, int fVerbose )
{
    extern int Kit_TruthToGia2( Gia_Man_t * p, unsigned * pTruth0, unsigned * pTruth1, int nVars, Vec_Int_t * vMemory, Vec_Int_t * vLeaves, int fHash );
    Gia_Man_t * pNew; Gia_Obj_t * pObj;
    Vec_Int_t * vMemory = Vec_IntAlloc( 1 << 18 );
    Vec_Int_t * vLeaves = Vec_IntAlloc( nIns );
    Vec_Wrd_t * vSimI = pFileName ? Vec_WrdReadBin( pFileName, 1 ) : NULL;  
    word * pTruth0 = ABC_CALLOC( word, Abc_Truth6WordNum(nIns) );
    word * pTruth1 = ABC_CALLOC( word, Abc_Truth6WordNum(nIns) ); int g, k; float CareAve = 0;
    if ( vSimI )
    {
        int nPats = 64*Vec_WrdSize(vSimI)/Gia_ManCiNum(p);
        printf( "Density of input  patterns %8.4f.\n", (float)Abc_TtCountOnesVec(Vec_WrdArray(vSimI), Vec_WrdSize(vSimI))/(64*Vec_WrdSize(vSimI)) );
        printf( "Using patterns with count %d and higher as cares (%8.4f %% of all patterns).\n", Thresh, 100.0*Thresh/nPats );
    }
    Gia_ManFillValue( p );
    pNew = Gia_ManStart( Gia_ManObjNum(p) );
    pNew->pName = Abc_UtilStrsav( p->pName );
    pNew->pSpec = Abc_UtilStrsav( p->pSpec );
    Gia_ManConst0(p)->Value = 0;
    Gia_ManForEachCi( p, pObj, k )
        pObj->Value = Gia_ManAppendCi(pNew);
    Gia_ObjComputeTruthTableStart( p, nIns );
    Gia_ManHashStart( pNew );
    for ( g = 0; g < Gia_ManCoNum(p); g += nOuts )
    {
        Vec_Int_t * vSupp = Gia_ManCollectSupp( p, g, nOuts );
        int Care, Temp = fVerbose ? printf( "Group %3d / %3d / %3d : Supp = %3d   %s", g, nOuts, Gia_ManCoNum(p), Vec_IntSize(vSupp), vSimI ? "":"\n" ) : 0;
        word * pCare = vSimI ? Gia_ManCountFraction( p, vSimI, vSupp, Thresh, fVerbose, &Care ) : ABC_FALLOC( word, Abc_Truth6WordNum(Vec_IntSize(vSupp)) );
        int nWords = Abc_Truth6WordNum( Vec_IntSize(vSupp) );
        CareAve += 100.0*Care/(1 << nIns);
        assert( Vec_IntSize(vSupp) <= nIns );
        Vec_IntClear( vLeaves );
        Gia_ManForEachObjVec( vSupp, p, pObj, k )
            Vec_IntPush( vLeaves, pObj->Value );        
        for ( k = 0; k < nOuts; k++ )
        {
            Gia_Obj_t * pObj = Gia_ManCo( p, g+k );
            if ( Gia_ObjIsAnd(Gia_ObjFanin0(pObj)) )
            {
                word * pTruth = Gia_ObjComputeTruthTableCut( p, Gia_ObjFanin0(pObj), vSupp );
                Abc_TtSharp( pTruth0, pCare, pTruth, nWords );
                Abc_TtAnd( pTruth1, pCare, pTruth, nWords, 0 );
                pObj->Value = Kit_TruthToGia2( pNew, (unsigned *)pTruth0, (unsigned *)pTruth1, Vec_IntSize(vLeaves), vMemory, vLeaves, 1 );
            }
            else
                pObj->Value = Gia_ObjFanin0(pObj)->Value;
            pObj->Value ^= Gia_ObjFaninC0(pObj);
        }
        ABC_FREE( pCare );
        Vec_IntFree( vSupp );
        Temp = 0;
    }
    CareAve /= Gia_ManCoNum(p)/nOuts;
    printf( "Average size of the care set = %8.4f %%.\n", CareAve );
    Gia_ManHashStop( pNew );
    Gia_ManForEachCo( p, pObj, k )
        pObj->Value = Gia_ManAppendCo( pNew, pObj->Value );
    Gia_ObjComputeTruthTableStop( p );
    ABC_FREE( pTruth0 );
    ABC_FREE( pTruth1 );
    Vec_IntFree( vLeaves );
    Vec_IntFree( vMemory );
    Vec_WrdFreeP( &vSimI );
    Gia_ManSetRegNum( pNew, Gia_ManRegNum(p) );
    return pNew;
}

#ifdef ABC_USE_CUDD

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
Gia_Man_t * Gia_ManDoMuxMapping( Gia_Man_t * p )
{
    extern Gia_Man_t * Gia_ManPerformMfs( Gia_Man_t * p, Sfm_Par_t * pPars );
    Gia_Man_t * pTemp, * pNew = Gia_ManDup( p );
    Jf_Par_t Pars, * pPars = &Pars; int c, nIters = 2;
    Sfm_Par_t Pars2, * pPars2 = &Pars2;
    Lf_ManSetDefaultPars( pPars );
    Sfm_ParSetDefault( pPars2 );
    pPars2->nTfoLevMax  =    5;
    pPars2->nDepthMax   =  100;
    pPars2->nWinSizeMax = 2000;
    for ( c = 0; c < nIters; c++ )
    {
        pNew = Lf_ManPerformMapping( pTemp = pNew, pPars );
        Gia_ManStop( pTemp );
        pNew = Gia_ManPerformMfs( pTemp = pNew, pPars2 );
        Gia_ManStop( pTemp );
        if ( c == nIters-1 )
            break;
        pNew = (Gia_Man_t *)Dsm_ManDeriveGia( pTemp = pNew, 0 );
        Gia_ManStop( pTemp );
    }
    return pNew;
}
Gia_Man_t * Gia_ManDoMuxTransform( Gia_Man_t * p, int fReorder )
{
    extern Gia_Man_t * Abc_NtkStrashToGia( Abc_Ntk_t * pNtk );
    extern int Abc_NtkBddToMuxesPerformGlo( Abc_Ntk_t * pNtk, Abc_Ntk_t * pNtkNew, int Limit, int fReorder, int fUseAdd );
    Gia_Man_t * pRes = NULL;
    Aig_Man_t * pMan = Gia_ManToAig( p, 0 );
    Abc_Ntk_t * pNtk = Abc_NtkFromAigPhase( pMan );
    Abc_Ntk_t * pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_SOP );
    pNtk->pName = Extra_UtilStrsav( pMan->pName );
    Aig_ManStop( pMan );
    //pNtkNew = Abc_NtkBddToMuxes( pNtk, 1, 1000000, 1 );
    if ( Abc_NtkBddToMuxesPerformGlo( pNtk, pNtkNew, 1000000, fReorder, 0 ) )
    {
        Abc_Ntk_t * pStrash = Abc_NtkStrash( pNtkNew, 1, 1, 0 );
        pRes = Abc_NtkStrashToGia( pStrash );
        Abc_NtkDelete( pStrash );
    }
    Abc_NtkDelete( pNtkNew );
    Abc_NtkDelete( pNtk );
    return pRes;
}
int Gia_ManDoTest1( Gia_Man_t * p, int fReorder )
{
    Gia_Man_t * pTemp, * pNew; int Res;
    pNew = Gia_ManDoMuxTransform( p, fReorder );
    pNew = Gia_ManDoMuxMapping( pTemp = pNew );
    Gia_ManStop( pTemp );
    Res = Gia_ManLutNum( pNew );
    Gia_ManStop( pNew );
    return Res;
}
Abc_Ntk_t * Gia_ManDoTest2( Gia_Man_t * p, int fReorder )
{
    extern Abc_Ntk_t * Abc_NtkFromMappedGia( Gia_Man_t * p, int fFindEnables, int fUseBuffs );
    Abc_Ntk_t * pNtkNew;
    Gia_Man_t * pTemp, * pNew;
    pNew = Gia_ManDoMuxTransform( p, fReorder );
    pNew = Gia_ManDoMuxMapping( pTemp = pNew );
    Gia_ManStop( pTemp );
    pNtkNew = Abc_NtkFromMappedGia( pNew, 0, 0 );
    pNtkNew->pName = Extra_UtilStrsav(p->pName);
    Gia_ManStop( pNew );
    Abc_NtkToSop( pNtkNew, 1, ABC_INFINITY );
    return pNtkNew;
}
Abc_Ntk_t * Abc_NtkMapTransform( Gia_Man_t * p, int nOuts, int fUseFixed, int fVerbose )
{
    extern Abc_Ntk_t * Abc_NtkSpecialMapping( Abc_Ntk_t * pNtk, int fVerbose );
    int i, k, g, nGroups = Gia_ManCoNum(p) / nOuts, CountsAll[3] = {0}; 
    Abc_Obj_t * pObjNew, * pFaninNew;  Gia_Obj_t * pObj;
    Abc_Ntk_t * pNtkNew = Abc_NtkAlloc( ABC_NTK_LOGIC, ABC_FUNC_SOP, 1 );
    assert( Gia_ManCoNum(p) % nOuts == 0 );
    pNtkNew->pName = Extra_UtilStrsav(p->pName);
    pNtkNew->pSpec = Extra_UtilStrsav(p->pSpec);
    Gia_ManFillValue( p );
    Gia_ManForEachPi( p, pObj, i )
        Abc_NtkCreatePi( pNtkNew );
    Gia_ManForEachPo( p, pObj, i )
        Abc_NtkCreatePo( pNtkNew );
    assert( nOuts <= 64 );
    for ( g = 0; g < nGroups; g++ )
    {
        Gia_Man_t * pNew;   Aig_Man_t * pMan;
        Abc_Ntk_t * pNtk, * pNtkRes, * pNtkMap;
        int pPos[64], Counter = 0, Counts[3] = {0};
        for ( i = 0; i < nOuts; i++ )
            pPos[i] = g*nOuts+i;
        pNew = Gia_ManDupCones( p, pPos, nOuts, 1 );
        if ( !fUseFixed )
            pNtkMap = Gia_ManDoTest2( pNew, 1 );
        else
        {
            pMan = Gia_ManToAig( pNew, 0 );
            pNtk = Abc_NtkFromAigPhase( pMan );
            Aig_ManStop( pMan );
            pNtkRes = Abc_NtkBddToMuxes( pNtk, 1, 1000000, 1 );
            Abc_NtkDelete( pNtk );
            pNtkMap = Abc_NtkSpecialMapping( pNtkRes, 0 );
            Abc_NtkDelete( pNtkRes );
        }
        Gia_ManStop( pNew );
        Gia_ManForEachCi( p, pObj, i )
            if ( ~pObj->Value )
                Abc_NtkCi(pNtkMap, Counter++)->pCopy = Abc_NtkCi(pNtkNew, i);
        assert( Counter == Abc_NtkCiNum(pNtkMap) );
        Abc_NtkForEachNode( pNtkMap, pObjNew, i )
        {
            pObjNew->pCopy = Abc_NtkDupObj( pNtkNew, pObjNew, 0 );
            pObjNew->pCopy->fPersist = pObjNew->fPersist;
            if ( pObjNew->fPersist )
                Counts[1]++;
            else
                Counts[0]++;
            Abc_ObjForEachFanin( pObjNew, pFaninNew, k )
                Abc_ObjAddFanin( pObjNew->pCopy, pFaninNew->pCopy );
        }
        Abc_NtkForEachCo( pNtkMap, pObjNew, i )
            Abc_ObjAddFanin( Abc_NtkCo(pNtkNew, g*nOuts+i), Abc_ObjFanin0(pObjNew)->pCopy );
        Abc_NtkDelete( pNtkMap );

        if ( fVerbose )
        {
        printf( "%3d / %3d :  ", g, nGroups );
        printf( "Test   = %4d   ", Counts[0] );
        printf( "MarkA  = %4d   ", Counts[1] );
        printf( "MarkB  = %4d   ", Counts[2] );
        printf( "\n" );
        }

        CountsAll[0] += Counts[0];
        CountsAll[1] += Counts[1];
        CountsAll[2] += Counts[2];
    }
    if ( fVerbose )
    printf( "Total LUT count = %5d.  MarkA = %5d. MarkB = %5d.\n", CountsAll[0], CountsAll[1], CountsAll[2] );
    // create names
    Abc_NtkAddDummyPiNames( pNtkNew );
    Abc_NtkAddDummyPoNames( pNtkNew );
    Abc_NtkAddDummyBoxNames( pNtkNew );
    // check the resulting AIG
    if ( !Abc_NtkCheck( pNtkNew ) )
        Abc_Print( 1, "Abc_NtkFromMappedGia(): Network check has failed.\n" );
    return pNtkNew;
}

Abc_Ntk_t * Gia_ManPerformLNetMap( Gia_Man_t * p, int GroupSize, int fUseFixed, int fVerbose )
{
    int fPrintOnly = 0;
    int Res1, Res2, Result = 0;
    int g, nGroups = Gia_ManCoNum(p) / GroupSize;
    assert( Gia_ManCoNum(p) % GroupSize == 0 );
    assert( GroupSize <= 64 );
    if ( fPrintOnly )
    {
        for ( g = 0; g < nGroups; g++ )
        {
            Gia_Man_t * pNew;
            int o, pPos[64];
            for ( o = 0; o < GroupSize; o++ )
                pPos[o] = g*GroupSize+o;
            pNew = Gia_ManDupCones( p, pPos, GroupSize, 0 );
            printf( "%3d / %3d :  ", g, nGroups );
            printf( "Test1 = %4d   ", Res1 = Gia_ManDoTest1(pNew, 0) );
            printf( "Test2 = %4d   ", Res2 = Gia_ManDoTest1(pNew, 1) );
            printf( "Test  = %4d   ", Abc_MinInt(Res1, Res2) );
            printf( "\n" );
            Result += Abc_MinInt(Res1, Res2);
            //Gia_ManPrintStats( pNew, NULL );
            Gia_ManStop( pNew );
        }
        printf( "Total LUT count = %d.\n", Result );
        return NULL;

    }
    return Abc_NtkMapTransform( p, GroupSize, fUseFixed, fVerbose );
}

#else

Gia_Man_t * Gia_ManPerformLNetMap( Gia_Man_t * p, int GroupSize, int fUseFixed, int fVerbose )
{
    return NULL;
}

#endif

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END

