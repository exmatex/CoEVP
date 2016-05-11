#ifndef included_KeyToString_h
#define included_KeyToString_h

#define uint128_t unsigned __int128
#include <string>

std::string to_string_hack(unsigned long long inVal);

std::string uint128_to_string(const uint128_t &in);

#endif
