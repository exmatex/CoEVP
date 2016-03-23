#ifndef __LULESH_H__
#define __LULESH_H__

#include "domain.h"
#include <list>
#include <chrono>
#define MAX_NONLINEAR_ITER 5
#define SEDOV_SYNC_POS_VEL_EARLY 1

class Lulesh {

private:

public:
  // Factor to be multiply the time step by to compensate
  // for fast time scales in the fine-scale model
  Real_t finescale_dt_modifier;
  int timer;
  int time_output;
  std::list<std::chrono::high_resolution_clock::time_point> timings;
  double write_timing;
//  std::chrono::duration<double> write_timing;
//  std::chrono::high_resolution_clock::time_point write_timing;
  std::ofstream timerfile;
  Domain domain;

Lulesh(){ finescale_dt_modifier = Real_t(1.); }
~Lulesh(){}

void CommRecv(Domain *domain, int msgType, Index_t xferFields, Index_t size,
         bool recvMin);
void CommSend(Domain *domain, int msgType,
         Index_t xferFields, Real_t **fieldData,
         Index_t *iset,  Index_t size, Index_t offset,
         bool sendMax);
void CommSBN(Domain *domain, int xferFields, Real_t **fieldData,
         Index_t *iset, Index_t size, Index_t offset);
void CommSyncPosVel(Domain *domain,
         Index_t *iset, Index_t size, Index_t offset);
void CommMonoQ(Domain *domain, Index_t *iset, Index_t size, Index_t offset);
void TimeIncrement();
void InitStressTermsForElems(Index_t numElem, 
         Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
         Real_t *sigxy, Real_t *sigxz, Real_t *sigyz);
void CalcElemShapeFunctionDerivatives( const Real_t* const x,
         const Real_t* const y,
         const Real_t* const z,
         Real_t b[][8],
         Real_t* const volume );                           
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
         Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
         Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
         Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
         const Real_t x0, const Real_t y0, const Real_t z0,
         const Real_t x1, const Real_t y1, const Real_t z1,
         const Real_t x2, const Real_t y2, const Real_t z2,
         const Real_t x3, const Real_t y3, const Real_t z3);
void CalcElemNodeNormals(Real_t pfx[8],
         Real_t pfy[8],
         Real_t pfz[8],
         const Real_t x[8],
         const Real_t y[8],
         const Real_t z[8]);
void SumElemStressesToNodeForces( const Real_t B[][8],
         const Real_t stress_xx,
         const Real_t stress_yy,
         const Real_t stress_zz,
         const Real_t stress_xy,
         const Real_t stress_xz,
         const Real_t stress_yz,
         Real_t* const fx,
         Real_t* const fy,
         Real_t* const fz );
void IntegrateStressForElems( Index_t numElem,
         Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
         Real_t *sigxy, Real_t *sigxz, Real_t *sigyz,
         Real_t *determ);
void CollectDomainNodesToElemNodes(const Index_t* elemToNode,
         Real_t elemX[8],
         Real_t elemY[8],
         Real_t elemZ[8]);
void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
         const Real_t x3, const Real_t x4, const Real_t x5,
         const Real_t y0, const Real_t y1, const Real_t y2,
         const Real_t y3, const Real_t y4, const Real_t y5,
         const Real_t z0, const Real_t z1, const Real_t z2,
         const Real_t z3, const Real_t z4, const Real_t z5,
         Real_t* dvdx, Real_t* dvdy, Real_t* dvdz);
void CalcElemVolumeDerivative(Real_t dvdx[8],
         Real_t dvdy[8],
         Real_t dvdz[8],
         const Real_t x[8],
         const Real_t y[8],
         const Real_t z[8]);
void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,  Real_t *hourgam0,
         Real_t *hourgam1, Real_t *hourgam2, Real_t *hourgam3,
         Real_t *hourgam4, Real_t *hourgam5, Real_t *hourgam6,
         Real_t *hourgam7, Real_t coefficient,
         Real_t *hgfx, Real_t *hgfy, Real_t *hgfz );
void CalcFBHourglassForceForElems(Real_t *determ,
         Real_t *x8n,      Real_t *y8n,      Real_t *z8n,
         Real_t *dvdx,     Real_t *dvdy,     Real_t *dvdz,
         Real_t hourg);
void CalcHourglassControlForElems(Real_t determ[], Real_t hgcoef);
void CalcVolumeForceForElems();
void CalcForceForNodes();
void CalcAccelerationForNodes();
void ApplyAccelerationBoundaryConditionsForNodes();
void CalcVelocityForNodes(const Real_t dt, const Real_t u_cut);
void CalcPositionForNodes(const Real_t dt);
void LagrangeNodal();
void LagrangeNodal1();
void LagrangeNodal2();
void CalcElemVelocityGradient( const Real_t* const xvel,
         const Real_t* const yvel,
         const Real_t* const zvel,
         const Real_t b[][8],
         const Real_t detJ,
         Real_t* const d, Real_t* const w );
void CalcKinematicsForElems( Index_t numElem, Real_t dt );
void CalcLagrangeElements(Real_t deltatime);
void CalcMonotonicQGradientsForElems();
void CalcMonotonicQRegionForElems(
         Real_t qlc_monoq,
         Real_t qqc_monoq,
         Real_t monoq_limiter_mult,
         Real_t monoq_max_slope,
         Real_t ptiny,

         // the elementset length
         Index_t elength );
void CalcMonotonicQForElems();
void CalcQForElems();
void CalcQForElems2();
void CalcPressureForElems(Real_t* p_new, Real_t* bvc,
         Real_t* pbvc, Real_t* e_old,
         Real_t* compression, Real_t *vnewc,
         Real_t pmin,
         Real_t p_cut, Real_t eosvmax,
         Index_t length);
void CalcEnergyForElems(Real_t* p_new, Real_t* e_new, Real_t* q_new,
         Real_t* bvc, Real_t* pbvc,
         Real_t* p_old, Real_t* e_old, Real_t* q_old,
         Real_t* compression, Real_t* compHalfStep,
         Real_t* vnewc, Real_t* work, Real_t* delvc, Real_t pmin,
         Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,
         Real_t* qq, Real_t* ql,
         Real_t rho0,
         Real_t eosvmax,
         Index_t length);
void CalcSoundSpeedForElems(Real_t *vnewc, Real_t rho0, Real_t *enewc,
         Real_t *pnewc, Real_t *pbvc,
         Real_t *bvc, Real_t ss4o3, Index_t nz);
void CalcWorkForElems(Real_t *vc, Real_t *work, Index_t length);
void EvalEOSForElems(Real_t *vnewc, Index_t length);
void ApplyMaterialPropertiesForElems();
void UpdateVolumesForElems();
void LagrangeElements();
void LagrangeElements2();
void CalcCourantConstraintForElems();
void CalcHydroConstraintForElems();
void CalcTimeConstraintsForElems();
void LagrangeLeapFrog();
void OutputTiming();
void FinalTime();
int UpdateStressForElems();
void UpdateStressForElems2(int reducedIters);
/*
void DumpDomainToVisit(DBfile *db, Domain& domain, int myRank);
void DumpMultiblockObjects(DBfile *db, char basename[], int numRanks);
void DumpToVisit(Domain& domain, char *baseName, char *meshName,
         int myRank, int numRanks);
void DumpSAMI(Domain *domain, char *name);
*/

void Initialize(int myRank, int numRanks, int edgeDim, int heightDim, double domainStopTime, int simStopCycle, int timerSampleRate);
void ConstructFineScaleModel(bool sampling,ModelDatabase * global_modelDB,ApproxNearestNeighbors* global_ann, int flanning, int flann_n_trees, int flann_n_checks, int global_ns);
void ExchangeNodalMass();
void go(int myRank, int numRanks, int sampling, int visit_data_interval,int file_parts, int debug_topology);

};

#endif
