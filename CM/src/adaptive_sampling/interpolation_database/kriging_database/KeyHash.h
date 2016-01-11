#ifndef included_KeyHash_h
#define included_KeyHash_h

#define uint128_t unsigned __int128

namespace std {

   // Defining a hash function in order to use uint128_t as a key
   // for an std::unordered_map
   template <>
      struct hash<uint128_t>
      {
         std::size_t operator()(const uint128_t& in) const
            {
               return (size_t)in;
            }
      };


}

#endif // included_KeyHash_h
