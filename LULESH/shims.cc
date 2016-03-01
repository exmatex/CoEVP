#include "domain.h"
#include "ElastoViscoPlasticity.h"

#if defined(PROTOBUF)
#include "shims.h"
#include "advance.pb.h"
#include <string>
#include <iostream>

#define SET_SCALAR
#define SET_TENSOR2SYM
#define SET_TENSOR2GEN
#define SET_DARRAY


//  This code stolen from Milo's Constitutive.h file
template <class T>
inline void pushObjectToBuffer(void** buffer, const T& data)
{
  size_t object_size = sizeof(T);
  memcpy(*buffer, &data, object_size);
  *buffer = ((char*)(*buffer)) + object_size;
};

template <class T>
inline void popObjectFromBuffer(void** buffer, T& data)
{
  size_t object_size = sizeof(T);
  memcpy(&data, *buffer, object_size);
  *buffer = ((char*)(*buffer)) + object_size;
};



std::string  &shimmy(Domain &domain, int k, std::string &sent_args, unsigned int ssize) {
  //  De-serialize the arguments so we can call the real advance()
  advance::Args  args;
  args.ParseFromString(sent_args);
  std::cout << args.DebugString();

  double      deltatime;
  Tensor2Gen  cm_vel_grad;
  double      cm_vol_chng;

  deltatime = args.delta_t();

  advance::Tensor_2_Gen  L_new = args.l_new();
  for (int i=0; i<9; i++) {
    cm_vel_grad.a[i] = L_new.a(i);
  }
  
  cm_vol_chng = args.volume_change();
  
  void                  *cm_state = ::operator new(ssize * sizeof(char));
  void                  *buffer = cm_state;
  advance::Tensor_2_Gen  t_2_g;     // temporaries
  advance::Tensor_2_Sym  t_2_s;
  Tensor2Sym             t2s;
  Tensor2Gen             t2g;
  double                 dub;
  int                    tint;

  advance::CM_State      _cm_state = args.cm_state();

  pushObjectToBuffer<double>(&buffer, _cm_state.m_delta_max());
  pushObjectToBuffer<double>(&buffer, _cm_state.m_delta());

  for (int i=0; i<6; i++) {
    t2s.a[i] = _cm_state.m_d_old().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = _cm_state.m_w_old().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<6; i++) {
    t2s.a[i] = _cm_state.m_d_new().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = _cm_state.m_w_new().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<9; i++) {
    t2g.a[i] = _cm_state.m_r().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, _cm_state.m_j());
  pushObjectToBuffer<int>(&buffer, _cm_state.m_num_iters());

  for (int i=0; i<6; i++) {
    t2s.a[i] = _cm_state.m_vbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = _cm_state.m_vbar_prime_dot().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = _cm_state.m_dbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = _cm_state.m_wbar().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, _cm_state.m_volume_change());
  pushObjectToBuffer<double>(&buffer, _cm_state.m_tol());
  pushObjectToBuffer<double>(&buffer, _cm_state.m_k());
  pushObjectToBuffer<double>(&buffer, _cm_state.m_g());

  for (int i=0; i<6; i++) {
    pushObjectToBuffer<double>(&buffer, _cm_state.fs_state().a(i));
  }

  
  ConstitutiveData cm_data =
    domain.cm(k)->advance(deltatime, cm_vel_grad, cm_vol_chng, cm_state);

  exit(1);
  //  Serialize cm_data for return
}


//  Test constitutive model's advance() call with:
//  - serialization of arguments
//  - call shim real advance() call
//    - deserialize arguments
//    - execute real advance() call
//    - serialize return values and modified arguments
//  - deserialize return (and modified arguments) from shim advance() call
ConstitutiveData  shim_advance(Domain &domain, int k) {
  std::string     send_args;                // string buffer to serialize to
  advance::Args   args;                     // protobuf class

  //  ======== Serialize arguments to advance call
                                                              // advance()'s argument names...
  args.set_delta_t(domain.deltatime());                       // delta_t

  for (int i=0; i<9; i++) {                                   // L_new
    args.mutable_l_new()->add_a(domain.cm_vel_grad(k).a[i]);
  }

  args.set_volume_change(domain.cm_vol_chng(k));              // volume_change

  advance::CM_State  *_cm_state = args.mutable_cm_state();    // state
  void* buffer = domain.cm_state(k);
  double      dub;      // temporaries 
  int         tint;
  Tensor2Sym  t2s;
  Tensor2Gen  t2g;

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_delta_max(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_delta(dub);
  
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    _cm_state->mutable_m_d_old()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    _cm_state->mutable_m_w_old()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    _cm_state->mutable_m_d_new()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    _cm_state->mutable_m_w_new()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    _cm_state->mutable_m_r()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_j(dub);

  popObjectFromBuffer<int>(&buffer, tint);
  _cm_state->set_m_num_iters(tint);

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    _cm_state->mutable_m_vbar_prime()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    _cm_state->mutable_m_vbar_prime_dot()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    _cm_state->mutable_m_dbar_prime()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    _cm_state->mutable_m_wbar()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_volume_change(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_tol(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_k(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  _cm_state->set_m_g(dub);

  for (int i=0; i<6; i++) {
    popObjectFromBuffer<double>(&buffer, dub);
    _cm_state->mutable_fs_state()->add_a(dub);
  }

  //  At this point, buffer is pointing to the end of the input state (except
  //  for fine-scale model state which we have yet to peel off. So, we can find the
  //  size of the passes on state and appropriately allocate it in shimme().

  unsigned int state_size = (char *)buffer - (char *)domain.cm_state(k);
  std::cout << "state size: " << state_size << endl;

  args.SerializeToString(&send_args);
  //std::cout << args.DebugString();

  //  ======== Call shim for advance() call
  std::string &retval = shimmy(domain, k, send_args, state_size);
  
  //  ======== Deserialize return values from shim advance() call
  //  build cm_data
  ConstitutiveData cm_data;
  return cm_data;
}

#endif
