#include "SingletonDB_Dummy.h"

#include <iostream>

SingletonDB_Dummy::SingletonDB_Dummy(bool dummyArg)
{
	std::cerr << "Error: Instantiating Dummy DB Backend. Check your build flags" << std::endl;
}

SingletonDB_Dummy::~SingletonDB_Dummy()
{
	std::cerr << "Error: Deleting Dummy DB Backend. Check your build flags" << std::endl;
}

void  SingletonDB_Dummy::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length)
{
	std::cerr << "Error: Pushing to Dummy DB Backend. Check your build flags" << std::endl;

}

void  SingletonDB_Dummy::erase(const uint128_t &key)
{
	std::cerr << "Error: Erasing from Dummy DB Backend. Check your build flags" << std::endl;

}

std::vector<double> SingletonDB_Dummy::pull(const uint128_t &key)
{
	std::cerr << "Error: Reading from Dummy DB Backend. Check your build flags" << std::endl;
	std::vector<double> retVec;
	return retVec;
}

std::vector<double> SingletonDB_Dummy::pull_key(const uint128_t &key)
{
	std::cerr << "Error: Reading from Dummy DB Backend. Check your build flags" << std::endl;
	std::vector<double> retVec;
	return retVec;
}
