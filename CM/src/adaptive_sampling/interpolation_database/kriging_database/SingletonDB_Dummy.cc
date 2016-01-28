#include "SingletonDB_Dummy.h"

#include <iostream>

SingletonDB_Dummy::SingletonDB_Dummy(bool dummyArg)
{
    throw std::runtime_error("Error: Instantiation Dummy DB Backend. Check your build flags"); 
}

SingletonDB_Dummy::~SingletonDB_Dummy()
{
    throw std::runtime_error("Error: Instantiation Dummy DB Backend. Check your build flags"); 
}

void  SingletonDB_Dummy::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length)
{
    throw std::runtime_error("Error: Instantiation Dummy DB Backend. Check your build flags"); 
}

void  SingletonDB_Dummy::erase(const uint128_t &key)
{
    throw std::runtime_error("Error: Instantiation Dummy DB Backend. Check your build flags"); 
}

std::vector<double> SingletonDB_Dummy::pull(const uint128_t &key)
{
    throw std::runtime_error("Error: Instantiation Dummy DB Backend. Check your build flags"); 
	std::vector<double> retVec;
	return retVec;
}

std::vector<double> SingletonDB_Dummy::pull_key(const uint128_t &key)
{
    throw std::runtime_error("Error: Instantiation Dummy DB Backend. Check your build flags"); 
	std::vector<double> retVec;
	return retVec;
}
