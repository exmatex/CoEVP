#include "siloDump.h"


#ifdef SILO

#ifdef __cplusplus
extern "C" {
#endif
#include "silo.h"
#ifdef __cplusplus
}
#endif

static void
DumpDomainToVisit(DBfile *db, Domain& domain, int myRank,
                  int beginPlane, int endPlane,
                  int sampling, int debug_topology)
{
   int ok = 0;
   int numElem    = domain.numElem() ;
   int numNode    = domain.numNode() ;
   int planeElems = numElem / domain.sliceHeight() ;
   int planeNodes = numNode / (domain.sliceHeight()+1) ;

   int *inodeMap = new int[numNode] ; /* maps actual nodeID to 'local' nodeID */

   if ((endPlane-beginPlane) != domain.sliceHeight() ) {
      numElem = planeElems * (endPlane-beginPlane) ;
      numNode = planeNodes * ((endPlane+1)-beginPlane) ;
   }

   int *elemMap = new int[numElem] ; /* maps 'local' elemID to actual elemID */
   int *nodeMap = new int[numNode] ; /* maps 'local' nodeID to actual nodeID */

   {
      int eOffset = 0 ;
      int ei = 0 ;
      for (int e=0; e<planeElems; ++e) {
         for (int p=beginPlane; p<endPlane; ++p) {
            elemMap[ei++] = eOffset + p ;
         }
         eOffset += domain.sliceHeight() ;
      }

      // if (ei != numElem) {
      //    exit(-1) ;
      // }
   }

   {
      int nOffset = 0 ;
      int ni = 0 ;
      for (int n=0; n<planeNodes; ++n) {
         for (int p=beginPlane; p<endPlane+1; ++p) {
            nodeMap[ni] = nOffset + p ;
            inodeMap[nOffset + p] = ni ; 
            ++ni ;
         }
         nOffset += domain.sliceHeight()+1 ;
      }

      // if (ni != numNode) {
      //    exit(-1) ;
      // }
   }


   /* Create an option list that will give some hints to VisIt for
      printing out the cycle and time in the annotations */
   DBoptlist *optlist;


   /* Write out the mesh connectivity in fully unstructured format */
   int shapetype[1] = {DB_ZONETYPE_HEX};
   int shapesize[1] = {8};
   int shapecnt[1] = {numElem};
   int *conn = new int[numElem*8] ;
   int ci = 0 ;
   for (int ei=0; ei < numElem; ++ei) {
      Index_t *elemToNode = domain.nodelist(elemMap[ei]) ;
      for (int ni=0; ni < 8; ++ni) {
         conn[ci++] = inodeMap[elemToNode[ni]] ;
      }
   }
   ok += DBPutZonelist2(db, "connectivity", numElem, 3,
                        conn, numElem*8,
                        0,0,0, /* Not carrying ghost zones */
                        shapetype, shapesize, shapecnt,
                        1, NULL);
   delete [] conn ;

   /* Write out the mesh coordinates associated with the mesh */
   const char* coordnames[3] = {"X", "Y", "Z"};
   Real_t *coords[3] ;
   coords[0] = new double[numNode] ;
   coords[1] = new double[numNode] ;
   coords[2] = new double[numNode] ;
   for (int ni=0; ni < numNode ; ++ni) {
      coords[0][ni] = domain.x(nodeMap[ni]) ;
      coords[1][ni] = domain.y(nodeMap[ni]) ;
      coords[2][ni] = domain.z(nodeMap[ni]) ;
   }
   optlist = DBMakeOptlist(2);
   DBClearOptlist(optlist) ;
   ok += DBAddOption(optlist, DBOPT_DTIME, &domain.time());
   ok += DBAddOption(optlist, DBOPT_CYCLE, &domain.cycle());
   ok += DBPutUcdmesh(db, "mesh", 3, (char**)&coordnames[0], (float**)coords,
                      numNode, numElem, "connectivity",
                      0, DB_DOUBLE, optlist);
   ok += DBFreeOptlist(optlist);
   delete [] coords[2] ;
   delete [] coords[1] ;
   delete [] coords[0] ;

   /* Write out the materials */
   int matnums = 1 ;
   int dims = numElem ; // No mixed elements
   int *regNumList = new int[numElem] ;

   for (int ei=0; ei<numElem; ++ei) {
      regNumList[ei] = 1 ;
   }

   ok += DBPutMaterial(db, "regions", "mesh", 1,
                       &matnums, regNumList, &dims, 1,
                       NULL, NULL, NULL, NULL, 0, DB_DOUBLE, NULL);
   delete [] regNumList;

   /* Write out pressure, energy, relvol, q */

   Real_t *e = new double[numElem] ;
   for (int ei=0; ei < numElem; ++ei) {
      e[ei] = domain.e(elemMap[ei]) ;
   }
   ok += DBPutUcdvar1(db, "e", "mesh", (float*) e,
                      numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                      NULL);
   delete [] e ;


   Real_t *p = new double[numElem] ;
   for (int ei=0; ei < numElem; ++ei) {
      p[ei] = domain.p(elemMap[ei]) ;
   }
   ok += DBPutUcdvar1(db, "p", "mesh", (float*) p,
                      numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                      NULL);
   delete [] p ;

   Real_t *mises = new double[numElem] ;
   for (int ei=0; ei < numElem; ++ei) {
      mises[ei] = domain.mises(elemMap[ei]) ;
   }
   ok += DBPutUcdvar1(db, "mises", "mesh", (float*) mises,
                      numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                      NULL);
   delete [] mises ;

   Real_t *v = new double[numElem] ;
   for (int ei=0; ei < numElem; ++ei) {
      v[ei] = domain.v(elemMap[ei]) ;
   }
   ok += DBPutUcdvar1(db, "v", "mesh", (float*) v,
                      numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                      NULL);
   delete [] v ;

   Real_t *volo = new double[numElem] ;
   for (int ei=0; ei < numElem; ++ei) {
      volo[ei] = domain.volo(elemMap[ei]) ;
   }
   ok += DBPutUcdvar1(db, "volo", "mesh", (float*) volo,
                      numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                      NULL);
   delete [] volo ;

   Real_t *q = new double[numElem] ;
   for (int ei=0; ei < numElem; ++ei) {
      q[ei] = domain.q(elemMap[ei]) ;
   }
   ok += DBPutUcdvar1(db, "q", "mesh", (float*) q,
                      numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                      NULL);
   delete [] q ;

   /* Write out nodal speed, velocities */
   Real_t *zd    = new double[numNode];
   Real_t *yd    = new double[numNode];
   Real_t *xd    = new double[numNode];
   Real_t *speed = new double[numNode];
   Real_t *nodalmass = new double[numNode];
   for(int ni=0 ; ni < numNode ; ++ni) {
      xd[ni]    = domain.xd(nodeMap[ni]);
      yd[ni]    = domain.yd(nodeMap[ni]);
      zd[ni]    = domain.zd(nodeMap[ni]);
      speed[ni] = sqrt((xd[ni]*xd[ni])+(yd[ni]*yd[ni])+(zd[ni]*zd[ni]));
      nodalmass[ni] = domain.nodalMass(nodeMap[ni]) ;
   }

   ok += DBPutUcdvar1(db, "speed", "mesh", (float*)speed,
                      numNode, NULL, 0, DB_DOUBLE, DB_NODECENT,
                      NULL);
   delete [] speed;


   ok += DBPutUcdvar1(db, "xd", "mesh", (float*) xd,
                      numNode, NULL, 0, DB_DOUBLE, DB_NODECENT,
                      NULL);
   delete [] xd ;

   ok += DBPutUcdvar1(db, "yd", "mesh", (float*) yd,
                      numNode, NULL, 0, DB_DOUBLE, DB_NODECENT,
                      NULL);
   delete [] yd ;

   ok += DBPutUcdvar1(db, "zd", "mesh", (float*) zd,
                      numNode, NULL, 0, DB_DOUBLE, DB_NODECENT,
                      NULL);
   delete [] zd ;

   ok += DBPutUcdvar1(db, "nodalmass", "mesh", (float*) nodalmass,
                      numNode, NULL, 0, DB_DOUBLE, DB_NODECENT,
                      NULL);
   delete [] nodalmass ;

   if (sampling) {

      Real_t *num_as_models = new double[numElem] ;
      Int_t numModels, numPairs;
      for (int ei=0; ei < numElem; ++ei) {
         domain.cm(elemMap[ei])->getModelInfo(numModels, numPairs);
         num_as_models[ei] = Real_t(numModels) ;
      }
      ok += DBPutUcdvar1(db, "num_as_models", "mesh", (float*) num_as_models,
                         numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                         NULL);
      delete [] num_as_models ;

      Real_t *as_efficiency = new double[numElem] ;
      for (int ei=0; ei < numElem; ++ei) {
         Int_t numSuccessful = domain.cm(elemMap[ei])->getNumSuccessfulInterpolations() ;
         Int_t numSamples = domain.cm(elemMap[ei])->getNumSamples() ;
         if ( numSamples > 0 ) {
            as_efficiency[ei] = Real_t(numSuccessful) / Real_t(numSamples) ;
         }
         else {
            as_efficiency[ei] = Real_t(1.) ;
         }
      }
      ok += DBPutUcdvar1(db, "as_efficiency", "mesh", (float*) as_efficiency,
                         numElem, NULL, 0, DB_DOUBLE, DB_ZONECENT,
                         NULL);

      delete [] as_efficiency;
   }

   if (debug_topology) {
      Index_t *nodeList = new int[numElem] ;

      for (int i=0 ; i<8; ++i) {
         char fieldName[40] ;
         sprintf(fieldName, "node_conn%d", i) ;
         for (Index_t j=0; j<numElem; ++j) {
            nodeList[j] = inodeMap[domain.nodelist(elemMap[j])[i]] ;
         }
         ok += DBPutUcdvar1(db, fieldName, "mesh", (int*) nodeList,
                         numElem, NULL, 0, DB_INT, DB_ZONECENT,
                         NULL);
      }

      delete [] nodeList ;


      char fname[][10] = { "lxim", "lxip", "letam", "letap", "lzetam", "lzetap" } ;
      Index_t *eList = new int[numElem] ;

      for (int i=0; i<6; ++i) {
         for (Index_t j=0; j<numElem; ++j) {
            switch (i) {
               case 0:
                  eList[j] = domain.lxim(elemMap[j]) ;
                  break ;
               case 1:
                  eList[j] = domain.lxip(elemMap[j]) ;
                  break ;
               case 2:
                  eList[j] = domain.letam(elemMap[j]) ;
                  break ;
               case 3:
                  eList[j] = domain.letap(elemMap[j]) ;
                  break ;
               case 4:
                  eList[j] = domain.lzetam(elemMap[j]) ;
                  break ;
               case 5:
                  eList[j] = domain.lzetap(elemMap[j]) ;
                  break ;
            }
         }
         ok += DBPutUcdvar1(db, fname[i], "mesh", (int*) eList,
                            numElem, NULL, 0, DB_INT, DB_ZONECENT,
                            NULL);
      }
      delete [] eList ;

      Index_t *eBC = new int[numElem] ;

      for (Index_t i=0; i<numElem; ++i) {
         eBC[i] = domain.elemBC(elemMap[i]) ;
      }
      ok += DBPutUcdvar1(db, "elembc", "mesh", (int*) eBC,
                         numElem, NULL, 0, DB_INT, DB_ZONECENT,
                         NULL);

      delete [] eBC ;
   }

   delete [] nodeMap ;
   delete [] elemMap ;
   delete [] inodeMap ;

   if (ok != 0) {
      printf("Error writing out viz file - rank %d\n", myRank);
   }
}

void DumpMultiblockObjects(DBfile *db, char basename[], int numRanks,
  int sampling, int debug_topology)
{
   /* MULTIBLOCK objects to tie together multiple files */
  char **multimeshObjs;
  char **multimatObjs;
  char ***multivarObjs;
  int *blockTypes;
  int *varTypes;
  int ok = 0;
  // Make sure this list matches what's written out above
  // All variables related to adaptive samplig MUST come AFTER the others
  char const *vars[] = {"p","e","mises","v","volo","q","speed","xd","yd","zd",
                        "num_as_models","as_efficiency",
                        "node_conn0", "node_conn1", "node_conn2", "node_conn3",
                        "node_conn4", "node_conn5", "node_conn6", "node_conn7",
                        "lxim", "lxip", "letam", "letap", "lzetam", "lzetap",
                        "elembc"
  };
  
  //  This is kinda hacky--find a cleaner way to handle this.
  int numvars = 10 ;
  if (sampling) {
    numvars += 2 ;
  }
  if (debug_topology) {
     /* hack */
     if (!sampling) {
        for (int i=0; i<15; ++i) {
           vars[numvars+i] = vars[numvars+i+2] ;
        }
     }
     numvars += 15 ;
  }

  // Reset to the root directory of the silo file
  DBSetDir(db, "/");

  // Allocate a bunch of space for building up the string names
  multimeshObjs = new char*[numRanks];
  multimatObjs = new char*[numRanks];
  multivarObjs = new char**[numvars];
  blockTypes = new int[numRanks];
  varTypes = new int[numRanks];

  for(int v=0 ; v<numvars ; ++v) {
     multivarObjs[v] = new char*[numRanks];
  }

  for(int i=0 ; i<numRanks ; ++i) {
     multimeshObjs[i] = new char[64];
     multimatObjs[i] = new char[64];
     for(int v=0 ; v<numvars ; ++v) {
        multivarObjs[v][i] = new char[64];
     }
     blockTypes[i] = DB_UCDMESH;
     varTypes[i] = DB_UCDVAR;
  }

  // Build up the multiobject names
  for(int i=0 ; i<numRanks ; ++i) {
    int iorank = i;

    //delete multivarObjs[i];
    if (iorank == 0) {
      snprintf(multimeshObjs[i], 64, "/data_%d/mesh", i);
      snprintf(multimatObjs[i], 64, "/data_%d/regions",i);
      for(int v=0 ; v<numvars ; ++v) {
        snprintf(multivarObjs[v][i], 64, "/data_%d/%s", i, vars[v]);
      }

    }
    else {
      snprintf(multimeshObjs[i], 64, "%s.%03d:/data_%d/mesh",
               basename, iorank, i);
      snprintf(multimatObjs[i], 64, "%s.%03d:/data_%d/regions",
               basename, iorank, i);
      for(int v=0 ; v<numvars ; ++v) {
         snprintf(multivarObjs[v][i], 64, "%s.%03d:/data_%d/%s",
                  basename, iorank, i, vars[v]);
      }
    }
  }

  // Now write out the objects
  ok += DBPutMultimesh(db, "mesh", numRanks,
                       (char**)multimeshObjs, blockTypes, NULL);
  ok += DBPutMultimat(db, "regions", numRanks,
                      (char**)multimatObjs, NULL);
  for(int v=0 ; v<numvars ; ++v) {
     ok += DBPutMultivar(db, vars[v], numRanks,
                         (char**)multivarObjs[v], varTypes, NULL);
  }

  for(int v=0; v < numvars; ++v) {
    for(int i = 0; i < numRanks; i++) {
      delete [] multivarObjs[v][i];
    }
    delete [] multivarObjs[v];
  }

  // Clean up
  for(int i=0 ; i<numRanks ; i++) {
    delete [] multimeshObjs[i];
    delete [] multimatObjs[i];
  }

  delete [] multimeshObjs;
  delete [] multimatObjs;
  delete [] multivarObjs;
  delete [] blockTypes;
  delete [] varTypes;

  if (ok != 0) {
    printf("Error writing out multiXXX objs to viz file - rank 0\n");
  }
}


void DumpToVisit(Domain& domain, char *baseName, char *meshName,
                 int myRank, int numRanks, int beginPlane, int endPlane,
		 int sampling, int debug_topology)
{
  DBfile *db;

  db = (DBfile*)DBCreate(meshName, DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5X);

  if (db) {
     char subdirName[128];

     sprintf(subdirName, "data_%d", myRank);
     DBMkDir(db, subdirName);
     DBSetDir(db, subdirName);
     DumpDomainToVisit(db, domain, myRank, beginPlane, endPlane, sampling, debug_topology);
     if (myRank == 0) {
        DumpMultiblockObjects(db, baseName, numRanks, sampling, debug_topology);
     }
     DBClose(db) ;
  }
  else {
     printf("Error writing out viz file - rank %d\n", myRank);
  }
}

void DumpDomain(Domain *domain, int myRank, int numProcs, int fileParts,
     int sampling, int debug_topology)
{
   char baseName[64] ;
   char meshName[64] ;

   /* set default slice information */
   int beginRank  = myRank ;
   int endRank    = myRank + 1 ;
   int beginPlane = 0 ;
   int endPlane   = domain->sliceHeight() ;


   if (fileParts != 0) {
      beginRank = 0 ;
      endRank   = fileParts ;
      numProcs  = fileParts ;
   }

   for (int rank = beginRank ; rank <endRank; ++rank) {

      if (fileParts != 0) {
         Index_t chunkSize = domain->sliceHeight() / fileParts ;
         Index_t remainder = domain->sliceHeight() % fileParts ;
         if (rank < remainder) {
            beginPlane = (chunkSize+1)*rank ;
            endPlane   = beginPlane + (chunkSize+1) ;
         }
         else {
            beginPlane = (chunkSize+1)*remainder + (rank - remainder)*chunkSize ;
            endPlane   = beginPlane + chunkSize ;
         }
      }

      sprintf(baseName, "taylor_%04d.silo", int(domain->cycle())) ;

      if (rank == 0) {
         sprintf(meshName, "%s", baseName) ;
      }
      else {
         sprintf(meshName, "%s.%03d", baseName, rank) ;
      }

      DumpToVisit(*domain, baseName, meshName, rank, numProcs, beginPlane, endPlane, sampling, debug_topology) ;
   }
}

#endif
