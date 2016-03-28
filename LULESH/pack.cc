//  This file contains two external entry points:
//  packAdvance()
//  unpackAdvance()
//  unpackAdvanceAndCall()
//
//  that use Google's protoc buffer serialization tool to pack and unpack
//  arguments to/from the CoEVP advance() call. The arguments are packed
//  to one long string.

#include "domain.h"
#include "ElastoViscoPlasticity.h"

#if defined(PROTOBUF)
#include "pack.pb.h"
#include <string>
#include <iostream>


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



//  Pack all arguments for advance() call
std::string  packAdvance(int k, ConstitutiveData &cm_data, double deltatime,
                         Tensor2Gen &cm_vel_grad, double cm_vol_chng, void *state) {
  pack::Args    args;              // protobuf object
  std::string   packed_string;     // serialize into this

  //  k
  args.set_k(k);
  //  deltatime
  args.set_delta_t(deltatime);
  //  cm_vel_grad
  for (int i=0; i<9; i++) {
    args.mutable_l_new()->add_a(cm_vel_grad.a[i]);
  }
  //  cm_vol_chng
  args.set_volume_change(cm_vol_chng);
  //  state
  double      dub;             // temporaries 
  int         tint;
  Tensor2Sym  t2s;
  Tensor2Gen  t2g;
  void       *buffer = state;  // so we can calculate length of state buffer
  //  state: m_delta_max
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_delta_max(dub);
  //  state: m_delta
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_delta(dub);
  //  state: m_d_old
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    args.mutable_cm_state()->mutable_m_d_old()->add_a(t2s.a[i]);
  }
  //  state: m_w_old
  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    args.mutable_cm_state()->mutable_m_w_old()->add_a(t2g.a[i]);
  }
  //  state: m_d_new
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    args.mutable_cm_state()->mutable_m_d_new()->add_a(t2s.a[i]);
  }
  //  state: m_w_new
  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    args.mutable_cm_state()->mutable_m_w_new()->add_a(t2g.a[i]);
  }
  //  state: m_r
  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    args.mutable_cm_state()->mutable_m_r()->add_a(t2g.a[i]);
  }
  //  state: m_j
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_j(dub);
  //  state: m_num_iters
  popObjectFromBuffer<int>(&buffer, tint);
  args.mutable_cm_state()->set_m_num_iters(tint);
  //  state: m_vbar_prime
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    args.mutable_cm_state()->mutable_m_vbar_prime()->add_a(t2s.a[i]);
  }
  //  state: m_vbar_prime_dot
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    args.mutable_cm_state()->mutable_m_vbar_prime_dot()->add_a(t2s.a[i]);
  }
  //  state: m_dbar_prime
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    args.mutable_cm_state()->mutable_m_dbar_prime()->add_a(t2s.a[i]);
  }
  //  state: m_wbar
  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    args.mutable_cm_state()->mutable_m_wbar()->add_a(t2g.a[i]);
  }
  //  state: m_volume_change
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_volume_change(dub);
  //  state: m_tol
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_tol(dub);
  //  state: m_k
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_k(dub);
  //  state: m_g
  popObjectFromBuffer<double>(&buffer, dub);
  args.mutable_cm_state()->set_m_g(dub);
  //  state: fs_state
  //  Fine-scale state (6 doubles for Taylor) is just tacked on to the end(!)
  for (int i=0; i<6; i++) {
    popObjectFromBuffer<double>(&buffer, dub);
    args.mutable_cm_state()->mutable_fs_state()->add_a(dub);
  }
  //  state_size
  //  At this point, buffer is pointing to the end of the input state.
  //  We use it to calculate the size of the buffer AND PASS IT BACK
  //  so subsequent unpacking can properly allocate a buffer to
  //  de-serialize into.
  unsigned int size = (char *)buffer - (char *)state;
  args.set_state_size(size);
  // std::cout << "state size: " << *size << endl;
  //  cm_data
  //  cm_data: sigma_prime
  for (int i=0; i<6; i++) {
    args.mutable_cm_data()->mutable_sigma_prime()->add_a(cm_data.sigma_prime.a[i]);
  }
  //  cm_data: num_models
  args.mutable_cm_data()->set_num_models(cm_data.num_models);
  //  cm_data: num_point_value_pairs
  args.mutable_cm_data()->set_num_point_value_pairs(cm_data.num_point_value_pairs);
  //  cm_data: num_newton_iters
  args.mutable_cm_data()->set_num_newton_iters(cm_data.num_Newton_iters);

  args.SerializeToString(&packed_string);
  return packed_string;
}



//  Unpack returned results from advance() call. Note that we are returning these
//  MULTIPLE results through pointers passed in to this function (k, cd, state).
//  They had better have been properly allocated and initialized in the caller!
void  unpackAdvance(std::string buf, int *k, ConstitutiveData *cd,
                    double *deltatime, Tensor2Gen *cm_vel_grad, double *cm_vol_chng, void *state) {
  //  First de-serialize the three (non-state) arguments
  pack::Args  args;
  args.ParseFromString(buf);
  // std::cout << args.DebugString();

  *k = args.k();
  *deltatime = args.delta_t();
  pack::Tensor_2_Gen  L_new = args.l_new();
  for (int i=0; i<9; i++) {
    cm_vel_grad->a[i] = L_new.a(i);
  }
  *cm_vol_chng = args.volume_change();

  //  WARNING: we are returning into the state pointer passed in. We
  //  reasonably assume that it is allocated to the correct size.
  //void  *state = ::operator new(state_size * sizeof(char));
  int    state_size = args.state_size();
  void  *buffer = state;

  pack::Tensor_2_Gen  t_2_g;     // temporaries
  pack::Tensor_2_Sym  t_2_s;
  Tensor2Sym             t2s;
  Tensor2Gen             t2g;
  double                 dub;
  int                    tint;

  pushObjectToBuffer<double>(&buffer, args.cm_state().m_delta_max());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_delta());

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_d_old().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_w_old().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_d_new().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_w_new().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_r().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, args.cm_state().m_j());
  pushObjectToBuffer<int>(&buffer, args.cm_state().m_num_iters());

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_vbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_vbar_prime_dot().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_dbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_wbar().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, args.cm_state().m_volume_change());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_tol());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_k());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_g());

  for (int i=0; i<6; i++) {
    pushObjectToBuffer<double>(&buffer, args.cm_state().fs_state().a(i));
  }

  unsigned int ssize = (char *)buffer - (char *)state;
  if (state_size != ssize) {
    std::cout << "Oops..." << state_size << " != " << ssize << std::endl;
    std::cout << "Something bound to fail!" << std::endl;
  }

  //  Unpack the ConstitutiveData portion of the return
  for (int i=0; i<6; i++) {
    cd->sigma_prime.a[i] = args.cm_data().sigma_prime().a(i);
  }
  cd->num_models = args.cm_data().num_models();
  cd->num_point_value_pairs = args.cm_data().num_point_value_pairs();
  cd->num_Newton_iters = args.cm_data().num_newton_iters();
}


//  This is a weird function because we need to test the operation of the pack
//  and unpack stuff by calling advance().  It's easiest to call advance() here,
//  rather than in LULESH itself.
//  So, while similar looking to the unpackAdvance function, this function will:
//
//  1) unpack the (previously packed--buf) arguments for advance()
//  2) call it using the additional domain structure parameter (containing the
//     function pointer to the advance() call
//  3) pack the results and return to LULESH
//
//  Other than the advance() call, and the pack at the end, it behaves much
//  as unpackAdvance() except that we don't pass anything back through the
//  argument list. We repack results into a string and let LULESH unpack it
//  just as if we'd called a remote instance of advance().
std::string  unpackAdvanceAndCall(Domain &domain, std::string buf, int *k,
                                  ConstitutiveData *cd, void *state) {
  //  First de-serialize the three (non-state) arguments
  pack::Args  args;
  args.ParseFromString(buf);
  // std::cout << args.DebugString();

  int         kk;
  double      deltatime;
  Tensor2Gen  cm_vel_grad;
  double      cm_vol_chng;

  kk = args.k();
  deltatime = args.delta_t();
  pack::Tensor_2_Gen  L_new = args.l_new();
  for (int i=0; i<9; i++) {
    cm_vel_grad.a[i] = L_new.a(i);
  }
  cm_vol_chng = args.volume_change();

  //  Allocate the void* buffer to a size determined a while back when
  //  it was being serialized.
  int    state_size = args.state_size();
  void  *sstate = ::operator new(state_size * sizeof(char));
  void  *buffer = sstate;

  pack::Tensor_2_Gen  t_2_g;     // temporaries
  pack::Tensor_2_Sym  t_2_s;
  Tensor2Sym             t2s;
  Tensor2Gen             t2g;
  double                 dub;
  int                    tint;

  pushObjectToBuffer<double>(&buffer, args.cm_state().m_delta_max());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_delta());

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_d_old().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_w_old().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_d_new().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_w_new().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_r().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, args.cm_state().m_j());
  pushObjectToBuffer<int>(&buffer, args.cm_state().m_num_iters());

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_vbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_vbar_prime_dot().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = args.cm_state().m_dbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = args.cm_state().m_wbar().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, args.cm_state().m_volume_change());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_tol());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_k());
  pushObjectToBuffer<double>(&buffer, args.cm_state().m_g());

  for (int i=0; i<6; i++) {
    pushObjectToBuffer<double>(&buffer, args.cm_state().fs_state().a(i));
  }

  unsigned int ssize = (char *)buffer - (char *)sstate;
  if (state_size != ssize) {
    std::cout << "Oops..." << state_size << " != " << ssize << std::endl;
    std::cout << "Something bound to fail!" << std::endl;
  }

  //  Unpack the ConstitutiveData portion of the return
  //  Don't really need to do this here since we're not using it for the call.
  ConstitutiveData  cm_data;
  for (int i=0; i<6; i++) {
    cm_data.sigma_prime.a[i] = args.cm_data().sigma_prime().a(i);
  }
  cm_data.num_models = args.cm_data().num_models();
  cm_data.num_point_value_pairs = args.cm_data().num_point_value_pairs();
  cm_data.num_Newton_iters = args.cm_data().num_newton_iters();

  ConstitutiveData ret_cm_data =
    domain.cm(kk)->advance(deltatime, cm_vel_grad, cm_vol_chng, sstate);

  return packAdvance(kk, ret_cm_data, deltatime, cm_vel_grad, cm_vol_chng, sstate);
}
#endif
