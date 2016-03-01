#include "domain.h"
#include "ElastoViscoPlasticity.h"

#if defined(PROTOBUF)
#include "shims.h"
#include "advance.pb.h"
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

//  Values returned from shimmy() are protobuf-serialized strings
struct ShimmyReturn {
  std::string  str_cm_data;
  std::string  str_cm_state;
};

//  Forward references
std::string           pack_state(void*, unsigned int*);
void                 *unpack_state(std::string, unsigned int);
struct ShimmyReturn   shimmy(Domain&, int, std::string&, std::string&, unsigned int);



//  Wrap constitutive model's advance() call with:
//  - serialization of arguments
//  - call shim real advance() call
//    - deserialize arguments
//    - execute real advance() call
//    - serialize return values and modified arguments
//  - deserialize return (and modified arguments) from shim advance() call
struct WrapReturn  *wrap_advance(Domain &domain, int k) {
  //  The real advance() call wants four arguments. We pack the first three into one
  //  string and the last (representing constitutive model and fine-scale model
  //  state) into another string. We use these two strings (selialized arguments)
  //  instead of one because the advance call passes a lot of data back through
  //  the state argument--so it's best to pack/unpack it independently.
  
  //  ======== Serialize first three arguments
  // advance()'s argument names: delta_t, L_new, volume_change
  std::string        args3;                   // string buffer to serialize to
  advance::Args      args;                     // protobuf class
  args.set_delta_t(domain.deltatime());
  for (int i=0; i<9; i++) {
    args.mutable_l_new()->add_a(domain.cm_vel_grad(k).a[i]);
  }
  args.set_volume_change(domain.cm_vol_chng(k));
  // std::cout << " ********* packed args in wrap_advance()" << std::endl;
  // std::cout << args.DebugString();
  args.SerializeToString(&args3);

  // ========= Serialize the state argument (big void* buffer)
  unsigned int  state_size;
  std::string  str_cm_state = pack_state(domain.cm_state(k), &state_size);
  //  state_size now has the length in bytes of the passed in void* buffer
  //  and it can be used in the remaining function calls--it does not change

  //  ======== Call real advance() via the shim
  struct ShimmyReturn shimmy_ret = shimmy(domain, k, args3, str_cm_state, state_size);

  //  ======== Deserialize return values from shimmy()
  //  We need to break the return data back into two pieces. The real
  //  advance() call returns a ConstitutiveData object AND additional data in
  //  the state's void* argument.
  struct WrapReturn *wrap_ret = new WrapReturn();
  wrap_ret->cm_data = new ConstitutiveData();

  //  Unpack the ConstitutiveData portion of the return
  advance::CM_Data cm_data;
  // std::cout << " ********* unpacking CM_Data in wrap_advance()" << std::endl;
  // std::cout << cm_data.DebugString();
  
  cm_data.ParseFromString(shimmy_ret.str_cm_data);
  for (int i=0; i<6; i++) {
    wrap_ret->cm_data->sigma_prime.a[i] = cm_data.sigma_prime().a(i);
  }
  wrap_ret->cm_data->num_models = cm_data.num_models();
  wrap_ret->cm_data->num_point_value_pairs = cm_data.num_point_value_pairs();
  wrap_ret->cm_data->num_Newton_iters = cm_data.num_newton_iters();
  
  //  Unpack the returned state into the required void*
  advance::CM_State cm_state;
  cm_state.ParseFromString(shimmy_ret.str_cm_state);
  // std::cout << " *** calling unpack with state_size: " << state_size << std::endl;
  wrap_ret->state = unpack_state(shimmy_ret.str_cm_state, state_size);

  // This is a nasty bit of code, but matches best with what actually happens
  // when LULESH calls advance(). The state, which is passed in as a void*
  // buffer is MODIFIED IN PLACE in advance (and hence "returned" to the
  // caller. So, we mimic that here by copying our void* buffer over the top
  // of the passed in buffer. We know the correct size from unpacking it.

  std::memcpy(wrap_ret->state, domain.cm_state(k), sizeof(char));
  
  return wrap_ret;
}


//  This is where we call the actual advance() function. You can imagine that this
//  function will eventually be in a different process, and likely executing on a
// different node than the caller. For that reason we need to have 1) serialized the
// arguments (which we have and are passed in here as two strings), 2) de-serialize
// them so that we 3) can call the actual advance() function with the argument as it
// expects them, and finally 4) pack up advance()'s return values (serialize) for the
// trip back the calling node where they will be 5) unpacked to the form that the
// original caller requires.  Note that in this function we only do items 2, 3, and 4.
struct ShimmyReturn shimmy(Domain &domain, int k, std::string &sent_args,
                           std::string &state, unsigned int ssize) {
  //  De-serialize the arguments so we can call the real advance()

  //  First de-serialize the three (non-state) arguments
  advance::Args  args;
  args.ParseFromString(sent_args);
  // std::cout << " ********* unpacking args in shimmy()" << std::endl;
  // std::cout << args.DebugString();

  double      deltatime;
  Tensor2Gen  cm_vel_grad;
  double      cm_vol_chng;

  deltatime = args.delta_t();
  advance::Tensor_2_Gen  L_new = args.l_new();
  for (int i=0; i<9; i++) {
    cm_vel_grad.a[i] = L_new.a(i);
  }
  cm_vol_chng = args.volume_change();
  
  //  Now, de-serialize the state arguments and stuff them in to the void*
  //  buffer that the real advance() call requires.
  // std::cout << " *** calling unpack with ssize: " << ssize << std::endl;
  void *cm_state = unpack_state(state, ssize);
  
  ConstitutiveData ret_cm_data =
    domain.cm(k)->advance(deltatime, cm_vel_grad, cm_vol_chng, cm_state);

  //  Pack cm_data for return.
  ShimmyReturn      shimmy_ret;
  advance::CM_Data  cm_data;
  unsigned int      tuint;          // temporaries
  Tensor2Sym        t2s;
  
  for (int i=0; i<6; i++) {
    cm_data.mutable_sigma_prime()->add_a(ret_cm_data.sigma_prime.a[i]);
  }
  cm_data.set_num_models(ret_cm_data.num_models);
  cm_data.set_num_point_value_pairs(ret_cm_data.num_point_value_pairs);
  cm_data.set_num_newton_iters(ret_cm_data.num_Newton_iters);
  cm_data.SerializeToString(&shimmy_ret.str_cm_data);
  // std::cout << " ********* packed CM_Data in shimmy()" << std::endl;
  // std::cout << cm_data.DebugString();

  //  Pack cm_state because lots of stuff is passed back through it.
  shimmy_ret.str_cm_state = pack_state(cm_state, &tuint);

  return shimmy_ret;
}




//  Serialize (pack) the state vector (in buffer) into a protobuf
//  encoded string. Return the size of the buffer (which is arrived at
//  using pointer arithmetic between its initial value and a pointer
//  that's left at the end of the buffer after numerous calls to
//  popObjectFromBuffer) so subsequent unpack functions can correctly
//  allocate a target void*.
std::string  pack_state(void *state, unsigned int *size) {  
  std::string         ret_state;     // serialize into this
  advance::CM_State   cm_state;      // protobuf object
  double      dub;      // temporaries 
  int         tint;
  Tensor2Sym  t2s;
  Tensor2Gen  t2g;

  void  *buffer = state;

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_delta_max(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_delta(dub);
  
  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    cm_state.mutable_m_d_old()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    cm_state.mutable_m_w_old()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    cm_state.mutable_m_d_new()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    cm_state.mutable_m_w_new()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    cm_state.mutable_m_r()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_j(dub);

  popObjectFromBuffer<int>(&buffer, tint);
  cm_state.set_m_num_iters(tint);

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    cm_state.mutable_m_vbar_prime()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    cm_state.mutable_m_vbar_prime_dot()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Sym>(&buffer, t2s);
  for (int i=0; i<6; i++) {
    cm_state.mutable_m_dbar_prime()->add_a(t2s.a[i]);
  }

  popObjectFromBuffer<Tensor2Gen>(&buffer, t2g);
  for (int i=0; i<9; i++) {
    cm_state.mutable_m_wbar()->add_a(t2g.a[i]);
  }

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_volume_change(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_tol(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_k(dub);

  popObjectFromBuffer<double>(&buffer, dub);
  cm_state.set_m_g(dub);

  //  Fine-scale state (6 doubles for Taylor) is just tacked on to the end(!)
  for (int i=0; i<6; i++) {
    popObjectFromBuffer<double>(&buffer, dub);
    cm_state.mutable_fs_state()->add_a(dub);
  }

  //  At this point, buffer is pointing to the end of the input state.
  //  We use it to calculate the size of the buffer AND PASS IT BACK
  //  so subsequent unpacking can properly allocate a buffer to
  //  de-serialize into.

  *size = (char *)buffer - (char *)state;
  // std::cout << "state size: " << *size << endl;
  
  // std::cout << " ********* packing state in pack_state()" << std::endl;
  // std::cout << cm_state.DebugString();
  cm_state.SerializeToString(&ret_state);
  
  return ret_state;
}



//  Unpack (de-serialize) state vector into a void* buffer (which is
//  how the real advance() call wants its state argument.
void  *unpack_state(std::string str_state, unsigned int ssize) {
  //  Allocate the void* buffer to a size determined a while back when
  //  it was being serialized.
  void                  *state = ::operator new(ssize * sizeof(char));
  void                  *buffer = state;

  advance::Tensor_2_Gen  t_2_g;     // temporaries
  advance::Tensor_2_Sym  t_2_s;
  Tensor2Sym             t2s;
  Tensor2Gen             t2g;
  double                 dub;
  int                    tint;

  advance::CM_State      cm_state;
  cm_state.ParseFromString(str_state);

  // std::cout << " ********* unpacking state in unpack_state()" << std::endl;
  // std::cout << cm_state.DebugString();

  pushObjectToBuffer<double>(&buffer, cm_state.m_delta_max());
  pushObjectToBuffer<double>(&buffer, cm_state.m_delta());

  for (int i=0; i<6; i++) {
    t2s.a[i] = cm_state.m_d_old().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = cm_state.m_w_old().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<6; i++) {
    t2s.a[i] = cm_state.m_d_new().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = cm_state.m_w_new().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  for (int i=0; i<9; i++) {
    t2g.a[i] = cm_state.m_r().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, cm_state.m_j());
  pushObjectToBuffer<int>(&buffer, cm_state.m_num_iters());

  for (int i=0; i<6; i++) {
    t2s.a[i] = cm_state.m_vbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = cm_state.m_vbar_prime_dot().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<6; i++) {
    t2s.a[i] = cm_state.m_dbar_prime().a(i);
  }
  pushObjectToBuffer<Tensor2Sym>(&buffer, t2s);

  for (int i=0; i<9; i++) {
    t2g.a[i] = cm_state.m_wbar().a(i);
  }
  pushObjectToBuffer<Tensor2Gen>(&buffer, t2g);

  pushObjectToBuffer<double>(&buffer, cm_state.m_volume_change());
  pushObjectToBuffer<double>(&buffer, cm_state.m_tol());
  pushObjectToBuffer<double>(&buffer, cm_state.m_k());
  pushObjectToBuffer<double>(&buffer, cm_state.m_g());

  for (int i=0; i<6; i++) {
    pushObjectToBuffer<double>(&buffer, cm_state.fs_state().a(i));
  }

  unsigned int state_size = (char *)buffer - (char *)state;
  if (state_size != ssize) {
    // std::cout << "Oops..." << state_size << " != " << ssize << std::endl;
    // std::cout << "Something bound to fail!" << std::endl;
  }

  return state;
}
#endif
