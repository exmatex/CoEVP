#include "KeyToString.h"

///TODO: Find a way to use to_string() if available  because this is bad
std::string to_string_hack(unsigned long long inVal)
{
	char buf[256];
	sprintf(buf, "%llu", inVal);
	std::string retString(buf);
	return buf;
}

std::string uint128_to_string(const uint128_t &in)
{
   uint64_t *in64 = (uint64_t *)&in; 
   return to_string_hack((unsigned long long) *in64)+to_string_hack((unsigned long long)*(in64+1));
}


