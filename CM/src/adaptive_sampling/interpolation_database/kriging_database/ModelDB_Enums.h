#ifndef included_ModelDB_Enums_h
#define included_ModelDB_Enums_h

enum SingletonDBBackendEnum
{
  REDIS_DB,
  HASHMAP_DB,
  POSIX_DB,
  HIO_DB,
  DIST_REDIS_DB
};

extern const char * SingletonDBBackendStrings[]; 


#endif
