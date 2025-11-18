#include "ergm_MHproposal.h"
#include "ergm_MHproposal_bd.h"
#include "ergm_dyadgen.h"
#include "ergm_MHstorage.h"
#include "ergm.multi_changestat_multilayer.h"

typedef struct{DyadGen *gen; DegreeBound *bd;} StoreDyadGenAndDegreeBound;

/*********************
 void MH_randomtoggleFixL

 Default MH algorithm with dyad generator API that enforces a layer constraint.
*********************/
MH_I_FN(Mi_randomtoggleFixL){
  ALLOC_STORAGE(1, StoreDyadGenAndDegreeBound, storage);
  storage->gen = DyadGenInitializeR(MHp->R, nwp, FALSE);
  storage->bd = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles = storage->gen->ndyads!=0 ? 1 : MH_FAILED;
}

MH_P_FN(MH_randomtoggleFixL){
  GET_STORAGE(StoreDyadGenAndDegreeBound, storage);
  BD_LOOP(storage->bd, {
      DyadGenRandDyad(Mtail, Mhead, storage->gen);
    });
  for(unsigned int ml=0; ml < MHp->n_aux; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    if(ergm_LayerLogic_affects(Mtail[0], Mhead[0], ll, (LayerLogicTask) TRUE, NULL, NULL)){
      Mtail[0]=MH_FAILED;
      Mhead[0]=MH_CONSTRAINT;
      return;
    }
  }
}

MH_F_FN(Mf_randomtoggleFixL){
  GET_STORAGE(StoreDyadGenAndDegreeBound, storage);
  DyadGenDestroy(storage->gen);
  DegreeBoundDestroy(storage->bd);
}

/********************
   void MH_TNT

   A TNT sampler that that enforces a layer constraint.

   A "intersect" network is constructed that is the intersection of
   dyads on the static list and the edges present in nwp. Then,
   standard TNT procedure is followed, but the dyad space (and the
   number of dyads) is the number of dyads in the static list and the
   network for the ties is the ties in the discord network.
***********************/

MH_I_FN(Mi_TNTFixL){
  ALLOC_STORAGE(1, StoreDyadGenAndDegreeBound, storage);
  storage->gen = DyadGenInitializeR(MHp->R, nwp, TRUE);
  storage->bd = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles = storage->gen->ndyads!=0 ? 1 : MH_FAILED;
}

MH_P_FN(Mp_TNTFixL){
  GET_STORAGE(StoreDyadGenAndDegreeBound, storage);

  const double P=0.5, Q=1-P;
  double DP = P*storage->gen->ndyads, DO = DP/Q;

  Edge nedges = DyadGenEdgecount(storage->gen);
  double logratio=0;
  BD_LOOP(storage->bd, {
      if (unif_rand() < P && nedges > 0) { /* Select a tie at random from the network of eligibles */
        DyadGenRandEdge(Mtail, Mhead, storage->gen);
	logratio = TNT_LR_E(nedges, Q, DP, DO);
      }else{ /* Select a dyad at random from the list */
	DyadGenRandDyad(Mtail, Mhead, storage->gen);

	if(IS_OUTEDGE(Mtail[0],Mhead[0])){
	  logratio = TNT_LR_DE(nedges, Q, DP, DO);
	}else{
	  logratio = TNT_LR_DN(nedges, Q, DP, DO);
	}
      }
    });
  for(unsigned int ml=0; ml < MHp->n_aux; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    if(ergm_LayerLogic_affects(Mtail[0], Mhead[0], ll, (LayerLogicTask)TRUE, NULL, NULL)){
      Mtail[0]=MH_FAILED;
      Mhead[0]=MH_CONSTRAINT;
      return;
    }
  }
  MHp->logratio += logratio;
}


MH_F_FN(Mf_TNTFixL){
  GET_STORAGE(StoreDyadGenAndDegreeBound, storage);
  DyadGenDestroy(storage->gen);
  DegreeBoundDestroy(storage->bd);
}
