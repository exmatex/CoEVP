from __future__ import print_function
import redis


def main():
    db = redis.StrictRedis(host='localhost', port=6379, db=0)
    all_keys = get_keys(db)
    for key in all_keys:
        print ('key: {:s}  val: {:s}'.format(key,db.get(key)))


# return a list of all keys in database
def get_keys(db):
    l = []
    for key in db.scan_iter():
        l.append(key)
    return l


# return a list of all elements in the set with argument key
def get_set(db, key):
    l = []
    for val in db.sscan_iter(key):
        l.append(val.decode("utf-8"))  # db sends bytes, we need strings
    return l


def get_all_sets(db):
    all_keys = get_keys(db)
    for key in all_keys:
        slist = get_set(db, key)


if __name__ == "__main__":
    main()
