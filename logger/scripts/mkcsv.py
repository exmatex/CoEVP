from __future__ import print_function
import redis


def main():
    db = redis.StrictRedis(host='localhost', port=6379, db=0)
    all_keys = get_keys(db)
    csv_file = open('dump.csv', 'w')
    init_csv_file(csv_file)
    for key in all_keys:
        add_csv_line(csv_file, key, db.get(key))
    csv_file.close()
    

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

        
def init_csv_file(f):
    f.write('Mode,Node,Rank,Timestep,Operation,Value,Executed\n')
    

def add_csv_line(f, k, v):
    k_list = k.split(',')
    for k_word in k_list:
        f.write(k_word + ', ')
    v_list = v.split(',')
    # [HACK] This relies on knowledge of the REDIS file format
    # Should just write to REDIS with comma separators ala CSV
    if k_list[0] == 'TIMER':
        f.write(v_list[0] + ', ' + v_list[1] + '\n')   # no timestamp!
    else:
        f.write(v_list[0] + ', '             + '\n')
        

if __name__ == "__main__":
    main()
