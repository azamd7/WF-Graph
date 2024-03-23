#include <iostream>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <string>
#include <signal.h>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <ctime> // std::time
#include <random>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <atomic>
#include <list>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>

using namespace std;

#define vntp "VERTEX NOT PRESENT"
#define entp "EDGE NOT PRESENT"
#define ventp "VERTEX OR EDGE NOT PRESENT"
#define ep "EDGE PRESENT"
#define eadd "EDGE ADDED"
#define er "EDGE REMOVED"
#define ef "EDGE FOUND"
#define vp "VERTEX PRESENT"
#define vadd "VERTEX ADDED"
#define vr "VERTEX REMOVED"

#include "common/recordmgr/record_manager.h"

int num_of_threads;
atomic<int> *shared_active_version_history = nullptr;

atomic<int> global_ts;

inline int is_marked_ref(long i)
{
    return (int)(i & 0x1L);
}

inline int is_freezed_ref(long i)
{
    return (int)(i & 0x2L);
}

inline long set_mark(long i)
{
    i |= 0x1L;
    return i;
}

inline long unset_mark(long i)
{
    i &= ~0x1L;
    return i;
}

inline long set_freezer(long i)
{
    i |= 0x2L;
    return i;
}

inline long unset_freezer(long i)
{
    i |= 0x2L;
    return i;
}

inline long get_marked_ref(long w)
{
    return w | 0x1L;
}

// both mark bit and freeze bits are removed
inline long get_unmarked_ref(long w)
{
    return w & ~0x3L;
}

inline long get_freezed_ref(long w)
{
    return w | 0x2L;
}

inline long get_unfreezed_ref(long w)
{
    return w & ~0x3L;
}

class e_version;
class v_version;

// Edge Node structure
class eNode
{
    public:
    int val;                                  // data
    atomic<e_version *> vhead; // pointer to its version head
    atomic<eNode*> enext;          // pointer to the next EdgeNode
    eNode()
    {
        vhead.store(nullptr);
        enext.store(nullptr);
    }
    eNode(int a)
    {
        val = a;
        vhead.store(nullptr);
        enext.store(nullptr);
    }
};

// (class LFListNode of FSet.h)
// Vertex Node structure
class vNode
{
public:
    int val;                                    // data
    vNode *vnext;                               // pointer to the next VertexNode
    atomic<v_version *> vhead; // version head
    vNode(int v)
    {
        val = v;
        vnext = nullptr;
        vhead.store(nullptr);
    }

    vNode(int v, vNode *vn, v_version *vh)
    {
        val = v;
        vnext = vn;
        vhead.store(vh);
    }
};

// version nodes for vertex
class v_version
{
    public:
    int val;
    atomic<int> insert_ts;
    atomic<int> delete_ts;
    atomic<eNode*> enext; // pointer to the EHead
    atomic<v_version*> next_ver;
    v_version()
    {
        enext.store(nullptr);
        next_ver.store(nullptr);
    }
};

// version nodes for edge
class e_version
{
    public:
    atomic<int> insert_ts;
    atomic<int> delete_ts;
    v_version *dest_version_node;
    int weight;
    atomic<e_version *> next_ver;
    e_version()
    {
        dest_version_node = (nullptr);
        next_ver.store(nullptr);
    }
};

// bucket class
class bucket
{
private:
    std::atomic<vNode *> head;

public:
    // default constructor
    bucket()
    {
        head = new vNode(-1);
    }

    // copy constructor
    bucket(vNode *h)
    {
        head.store(h);
    }

    // function to read head of the bucket
    vNode *get_bucket_head()
    {
        return head.load();
    }

    // function to cas the head
    bool casHead(vNode *a, vNode *b)
    {
        return head.compare_exchange_strong(a, b);
    }

    // helps in splitting the elements in the bucket to different buckets
    vNode *split(int size, int mod_remainder)
    {
        // requires this passed fset id immutable
        vNode *copy = new vNode(INT_MAX);
        // vNode *curr = head.load();
        vNode *before_v = head.load();
        vNode *curr = (vNode *)get_unfreezed_ref((long)head.load());
        vNode *after_v = curr;

        while (curr->val != INT_MAX) // until tail sentinel is reached
        {
            if (curr->val % size == mod_remainder)
            {
                vNode *n = new vNode(curr->val, copy, curr->vhead.load());
                copy = n;
            }
            curr = curr->vnext;
        }
        return copy;
    }

    // need to handle sen
    vNode *merge(bucket *b2)
    {
        // requires this and t2 are immutable
        vNode *copy = nullptr;
        for (int p = 0; p < 2; p++)
        {
            // std::vector<int> V;
            vNode *curr = (p == 0) ? head.load() : b2->head.load();

            while (curr->val != INT_MAX) // until tail sentinel is reached
            {
                vNode *n = new vNode(curr->val, copy, curr->vhead.load());
                copy = n;
                curr = curr->vnext;
            }
        }
        return copy;
    }
};

class BFS_node
{
public:
    int val; // value/key of the node
    BFS_node *parent = nullptr;
    BFS_node *child_head = nullptr;
    BFS_node *sibling = nullptr;

    BFS_node(int a, BFS_node *b)
    {
        val = a;
        parent = b;
        child_head = nullptr;
        sibling = nullptr;
    }

    BFS_node(int a, BFS_node *b, BFS_node *c)
    {
        val = a;
        parent = b;
        child_head = nullptr;
        sibling = c;
    }
};

class Queue_entry
{
public:
    v_version *vertex_version;
    BFS_node *vertex_bfs_node;

    Queue_entry(v_version *a, BFS_node *b)
    {
        vertex_version = a;
        vertex_bfs_node = b;
    }
};

// class for priority queue elements
template <typename T>
class SSSP_PQ_node
{
public:
    int val;
    long distance;
    T *dest_vertex;

    void operator=(const SSSP_PQ_node &C)
    {
        val = C.val;
        distance = C.distance;
        dest_vertex = C.dest_vertex;
    }

    SSSP_PQ_node(int a, long b, T *c)
    {
        this->val = a;
        this->distance = b;
        this->dest_vertex = c;
    }
};

// operator overloading for priority queue (returns true if v2 is closer than v1)
template <typename T>
class Compare_Distance
{
public:
    bool operator()(SSSP_PQ_node<T> const &v1, SSSP_PQ_node<T> const &v2)
    {
        return v1.distance > v2.distance;
    }
};

class LFGraph
{
    // hash class
    struct HNode
    {
        std::atomic<HNode *> old;
        bucket **buckets_array = nullptr;
        // std::vector<bucket *> buckets_array;
        int size;

        HNode(HNode *o, int s)
        {
            this->old.store(o);
            this->size = s;
            this->buckets_array = new bucket *[s];
            // this->buckets_array.resize(s);
            for (int i = 0; i < s; i++)
            {
                vNode *Tail = (vNode *)malloc(sizeof(vNode));
                Tail->val = INT_MAX;
                Tail->vnext = nullptr;
                Tail->vhead.store(nullptr);
                this->buckets_array[i] = new bucket(Tail);
            }
        }
    };

    static const int MIN_BUCKET_NUM = 1;
    static const int MAX_BUCKET_NUM = 1 << 16;

public:
    std::atomic<HNode *> hash_head;

    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, e_version>* myRecManager_e_ver;
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, v_version>* myRecManager_v_ver;
    //record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, vNode>* myRecManager_vertex;
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, eNode>* myRecManager_edge;

    LFGraph(record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, e_version> * myRecManager1,
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, v_version> * myRecManager2,
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, eNode> * myRecManager4){
        hash_head.store(new HNode(nullptr, MIN_BUCKET_NUM));
        // hash_head.load()->buckets_array[0] = (new bucket());
        myRecManager_e_ver = myRecManager1;
        myRecManager_v_ver = myRecManager2;
        myRecManager_edge = myRecManager4;
        global_ts.store(1);
    }

    HNode *get_hash_head()
    {
        return hash_head.load();
    }

    bool casHead(HNode *o, HNode *n)
    {
        return hash_head.compare_exchange_strong(o, n);
    }

    // creation of new Edge Node
    eNode *createE(int key, int tid)
    {
        eNode *newe = myRecManager_edge->template allocate<eNode>(tid);
        newe->val = key;
        newe->vhead.store(nullptr);
        newe->enext.store(nullptr);
        return newe;
    }

    // creation of new version node for vertices
    v_version *create_vertex_Ver(int key, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        eNode *EHead = createE(INT_MIN, tid); // create Edge Head
        eNode *ETail = createE(INT_MAX, tid); // create Edge Tail
        EHead->enext.store(ETail);
        //v_version *newv = (v_version *)malloc(sizeof(v_version));
        v_version *newv = myRecManager_v_ver->template allocate<v_version>(tid);
        newv->val = key;
        newv->delete_ts = INT_MAX;
        newv->insert_ts = -1;
        newv->enext.store(EHead);
        newv->next_ver.store(nullptr);
        return newv;
    }

    e_version *create_edge_Ver(v_version *dest, int weight, int tid)
    {
        //e_version *newv = (e_version *)malloc(sizeof(e_version));
        e_version *newv = myRecManager_e_ver->template allocate<e_version>(tid);
        newv->delete_ts = INT_MAX;
        newv->insert_ts = -1;
        newv->weight = weight;
        newv->next_ver.store(nullptr);
        newv->dest_version_node = dest;
        return newv;
    }

    // creation of new Vertex Node
    vNode *createV(int key, int tid)
    {
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        vNode *newv = (vNode *)malloc(sizeof(vNode));
        newv->val = key;
        newv->vnext = nullptr;
        newv->vhead = create_vertex_Ver(key, tid);
        return newv;
    }

    // Find pred and curr for eNode(key)
    void locateE(eNode *startE, eNode **n1, eNode **n2, int key)
    {
        eNode *curre, *prede;
        vNode *tv;
        prede = startE;
        curre = prede->enext.load();
        while (curre)
        {
            if (curre->val >= key)
            {
                (*n1) = prede;
                (*n2) = curre;
                return;
            }
            prede = curre;
            curre = curre->enext.load();
        }
    }

    // sets timestamp on the vertex's version node
    void setTs(v_version *ver, int tid)
    {
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        int t = -1;
        if (ver == nullptr)
            return;
        if (ver->insert_ts == -1)
            ver->insert_ts.compare_exchange_strong(t, global_ts);
        t = -1;
        if (ver->delete_ts == -1)
            ver->delete_ts.compare_exchange_strong(t, global_ts);
    }

    // sets timestamp on the edge's version node
    void setTs(e_version *ver, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        int t = -1;
        if (ver == nullptr)
            return;
        if (ver->insert_ts == -1)
            ver->insert_ts.compare_exchange_strong(t, global_ts);
        t = -1;
        if (ver->delete_ts == -1)
            ver->delete_ts.compare_exchange_strong(t, global_ts);
    }

    int getMyTs()
    {
        int temp1 = global_ts.load();
        int temp2 = temp1;
        global_ts.compare_exchange_strong(temp1, temp1 + 1);
        return temp2;
    }

    // read the version head of a vertex
    v_version *read(vNode *V, int tid)
    {
        v_version *curr_version = V->vhead.load();
        setTs(curr_version, tid);
        return curr_version;
    }

    // read the version head of an edge
    e_version *read(eNode *V, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        e_version *curr_version = V->vhead.load();
        if (curr_version == nullptr)
            return curr_version;
        setTs(curr_version, tid);
        return curr_version;
    }

    // find smallest active timestamp among all threads (for GARBAGE COLLECTION)
    int find_min_TS(fstream *logfile, bool debug, int tid)
    {
        int min_TS = INT_MAX, temp;
        for(int i=0;i<num_of_threads;i++){
            temp = shared_active_version_history[i].load();
            if(temp == -1){
                int new_TS = getMyTs();
                int exp_val = -1;
                shared_active_version_history[i].compare_exchange_strong(exp_val, new_TS);
                temp = shared_active_version_history[i].load();
            }
            if(temp < min_TS){
                min_TS = temp;
            }
        }
        int curr_global_TS = global_ts.load();
        return ((curr_global_TS < min_TS) ? curr_global_TS : min_TS);
    }

    // returns the appropriate version by storing it in return_version if it exists for the asked key
    bool ContainsV(int key, v_version **return_version, fstream *logfile, bool debug, int tid)
    {
        HNode *t = hash_head.load();
        vNode *b = t->buckets_array[key % t->size]->get_bucket_head();
        int temp = b->val;
        if (b->val == INT_MAX)
        {
            HNode *s = t->old.load();
            // helpResize(t, key % t->size, 0);
            b = (s == nullptr) ? t->buckets_array[key % t->size]->get_bucket_head() : s->buckets_array[key % s->size]->get_bucket_head();
        }
        b = (vNode *)get_unfreezed_ref((long)b);
        while (b->val != INT_MAX)
        {
            // cout<<"From ContainsV b = : "<<b<<endl;
            temp = b->val;
            if (b->val == key)
            {
                *return_version = read(b, tid);
                // cout<<"From ContainsV return_version = : "<<return_version<<endl;
                if ((*return_version)->delete_ts == INT_MAX)
                {
                    if (debug)
                        (*logfile) << "Vertex found." << endl;
                    return true;
                }
                else
                {
                    setTs(b->vhead.load(), tid);
                    return_version = nullptr;
                    return false;
                }
            }
            b = b->vnext;
        }
        if (debug)
            (*logfile) << "Vertex NOT found." << endl;
        *return_version = nullptr;
        return false;
    }

    // Returns version of both the vertices(in n1 and n2), if both the vertices have active version
    bool ContainsVPlus(v_version **n1, v_version **n2, int key1, int key2, fstream *logfile, bool debug, int tid)
    {
        v_version *curr_ver1, *curr_ver2;
        (*n1) = nullptr;
        (*n2) = nullptr;
        if (ContainsV(key1, &curr_ver1, logfile, debug, tid) && ContainsV(key2, &curr_ver2, logfile, debug, tid))
        {
            if (curr_ver1->delete_ts == INT_MAX && curr_ver2->delete_ts == INT_MAX)
            {
                *n1 = curr_ver1;
                *n2 = curr_ver2;
                return true;
            }
        }
        return false;
    }

//    // GC routine for edge versions (GET IT CHECKED BY """""SIR""""")
//    void do_GC_edge_ver(e_version *curr, int my_valid_GC_timestamp, fstream *logfile, bool debug, int tid){
//       auto guard1 = myRecManager_e_ver->getGuard(tid);
//       e_version *prev = curr;
//       curr = curr->next_ver.load();
//       e_version *curr_E_ver_temp = curr;
//        while(curr)
//        {
//            curr_E_ver_temp = curr;
//            curr = curr->next_ver.load();
//            if((curr_E_ver_temp->delete_ts != -1) && (curr_E_ver_temp->delete_ts < my_valid_GC_timestamp))
//            {
//                prev->next_ver.store(nullptr);
//                myRecManager_e_ver->template retire(tid, curr_E_ver_temp);
//                while(curr)     // to avoid the IF condtion from here on
//                {
//                    curr_E_ver_temp = curr;
//                    prev = curr;
//                    curr = curr->next_ver.load();
//                    prev->next_ver.store(nullptr);
//                    myRecManager_e_ver->template retire(tid, curr_E_ver_temp);
//                }
//            }
//        }
//    }

    // add a new edge in the edge-list
    bool AddE(int key1, int key2, int weight, fstream *logfile, bool debug, int tid)
    {
        
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        int my_valid_GC_timestamp = find_min_TS(logfile, debug, tid);
        eNode *prede, *curre;
        v_version *u, *v;
        e_version *curr_ver, *prev_ver;
        while (true)
        {
            bool flag = ContainsVPlus(&u, &v, key1, key2, logfile, debug, tid);
            if (flag == false)
            { // either of the vertex is not present
                if (debug)
                    (*logfile) << key1 << " - " << key2 << " Either of the vertex is NOT present." << endl;
                return false;
            }

            eNode *abc1 = u->enext.load();
            eNode *abc2 = v->enext.load();
            locateE(u->enext.load(), &prede, &curre, key2);
            if (u->delete_ts != INT_MAX || v->delete_ts != INT_MAX)
            {
                setTs(u, tid);
                setTs(v, tid);
                if (debug)
                    (*logfile) << key1 << " - " << key2 << " Either of the vertex is NOT present." << endl;
                return false;
            }
            if (curre->val != key2)
            { // if edge was never present for this version of vertex
                eNode *node = createE(key2, tid);
                node->enext.store(curre);
                node->vhead.store(create_edge_Ver(v, weight, tid));
                e_version *curr_version = node->vhead;
                if (prede->enext.compare_exchange_strong(curre, node))
                {
                    setTs(curr_version, tid);
                    int ets = curr_version->insert_ts;
                    int sts = u->delete_ts; // source vertex deleted concurrently as indicated by source vertex's version u
                    if (sts == -1)
                    {
                        setTs(u, tid);
                        sts = u->delete_ts;
                    }
                    int dts = v->delete_ts;
                    if (dts == -1)
                    { // dest vertex deleted concurrently as indicated by dest vertex's version v
                        setTs(v, tid);
                        dts = v->delete_ts;
                    }
                    if (ets < sts && ets < dts)
                    { // IF edge's insert_ts is smaller than source's/dest's delete_ts [sorted]
                        if (debug)
                            (*logfile) << key1 << " - " << key2 << " Edge ADDED." << endl;
                        return true;
                    }
                    else
                    {
                        setTs(u, tid);
                        setTs(v, tid);
                        if (debug)
                            (*logfile) << key1 << " - " << key2 << " Either of the vertex GOT DELETED concurrently." << endl;
                        return false;
                    }
                }
            }
            else
            { // if edge was present at some point for this version of vertex
                curr_ver = read(curre, tid);
                if (curr_ver->delete_ts == INT_MAX && curr_ver->weight == weight)
                { // edge already present
                    if (debug)
                        (*logfile) << key1 << " - " << key2 << " Edge ALREADY present." << endl;
                    // GARBAGE COLLECTION
                    //(curr_ver, my_valid_GC_timestamp, logfile, debug, tid);
                    return false;
                }
                else
                {
                    setTs(curr_ver, tid);
                    e_version *new_version = create_edge_Ver(v, weight, tid);
                    new_version->next_ver.store(curr_ver);
                    if (curre->vhead.compare_exchange_strong(curr_ver, new_version))
                    {
                        setTs(new_version, tid);
                        int ets = new_version->insert_ts;
                        int sts = u->delete_ts;
                        if (sts == -1)
                        {
                            setTs(u, tid);
                            sts = u->delete_ts;
                        }
                        int dts = v->delete_ts;
                        if (dts == -1)
                        {
                            setTs(v, tid);
                            dts = v->delete_ts;
                        }
                        if (ets < sts && ets < dts)
                        {
                            if (debug)
                                (*logfile) << key1 << " - " << key2 << " Edge ADDED." << endl;
                            // GARBAGE COLLECTION
                            //do_GC_edge_ver(curr_ver, my_valid_GC_timestamp, logfile, debug, tid);
                            return true;
                        }
                        
                    }
                    if (debug)
                        (*logfile) << key1 << " - " << key2 << " Either of the vertex GOT DELETED concurrently." << endl;
                    // GARBAGE COLLECTION
                    //do_GC_edge_ver(curr_ver, my_valid_GC_timestamp, logfile, debug, tid);
                    return false;
                }
            }
        }
    }

    // Deletes an edge from the edge-list if present
    bool RemoveE(int key1, int key2, fstream *logfile, bool debug, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        int my_valid_GC_timestamp = find_min_TS(logfile, debug, tid);
        eNode *prede, *curre;
        v_version *u, *v;
        e_version *curr_ver;
        bool flag = ContainsVPlus(&u, &v, key1, key2, logfile, debug, tid);
        if (flag == false)
        {
            if (debug)
                (*logfile) << "Either of the vertex is NOT present." << endl;
            return false; // either of the vertex is not present
        }
        locateE(u->enext.load(), &prede, &curre, key2);
        if (u->delete_ts != INT_MAX || v->delete_ts != INT_MAX)
        {
            setTs(u, tid);
            setTs(v, tid);
            if (debug)
                (*logfile) << "Either of the vertex is NOT present." << endl;
            return false;
        }
        if (curre->val != key2)
        {
            if (debug)
                (*logfile) << "Edge NOT present." << endl;
            return false;
        }
        else
        {
            curr_ver = read(curre, tid);
            if (curr_ver->delete_ts != INT_MAX)
            {
                setTs(curr_ver, tid);
                if (debug)
                    (*logfile) << "Edge ALREADY removed (-1 beforehand)." << endl;
                // GARBAGE COLLECTION
                //do_GC_edge_ver(curr_ver, my_valid_GC_timestamp, logfile, debug, tid);
                return false;
            }
            else
            {
                int t = INT_MAX;
                if (curr_ver->delete_ts.compare_exchange_strong(t, -1))
                {
                    setTs(curr_ver, tid);
                    int ets = curr_ver->delete_ts;
                    int sts = u->delete_ts;
                    if (sts == -1)
                    {
                        setTs(u, tid);
                        sts = u->delete_ts;
                    }
                    int dts = v->delete_ts;
                    if (dts == -1)
                    {
                        setTs(v, tid);
                        dts = v->delete_ts;
                    }
                    if (ets < sts && ets < dts)
                    {
                        if (debug)
                            (*logfile) << "Edge REMOVED." << endl;
                        // GARBAGE COLLECTION
                        //do_GC_edge_ver(curr_ver, my_valid_GC_timestamp, logfile, debug, tid);
                        return true;
                    }
                }
                setTs(curr_ver, tid);
                // GARBAGE COLLECTION
                //do_GC_edge_ver(curr_ver, my_valid_GC_timestamp, logfile, debug, tid);
                return false;
            }
        }
    }

    // Contains Edge Node
    bool ContainsE(int key1, int key2, fstream *logfile, bool debug, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        eNode *prede, *curre;
        v_version *u, *v;
        e_version *curr;
        bool flag = ContainsVPlus(&u, &v, key1, key2, logfile, debug, tid);
        if (flag == false)
        {
            if (debug)
                (*logfile) << "Either of the vertex is NOT present(1)." << endl;
            return false; // either of the vertex is not present
        }
        locateE(u->enext.load(), &prede, &curre, key2);
        if (u->delete_ts != INT_MAX || v->delete_ts != INT_MAX)
        {
            setTs(u, tid);
            setTs(v, tid);
            if (debug)
                (*logfile) << "Either of the vertex is NOT present(2)." << endl;
            return false;
        }
        if (curre->val != key2)
        {
            if (debug)
                (*logfile) << "Edge NOT present(1)." << endl;
            return false;
        }
        else
        {
            curr = read(curre, tid);
            int ets = curr->insert_ts;
            if (curr->delete_ts != INT_MAX)
            {
                setTs(curr, tid);
                if (debug)
                    (*logfile) << "Edge DELETED concurrently, therefore NOT present(2)." << endl;
                return false;
            }
            int sts = u->delete_ts;
            if (sts == -1)
            {
                setTs(u, tid);
                sts = u->delete_ts;
            }
            int dts = v->delete_ts;
            if (dts == -1)
            {
                setTs(v, tid);
                dts = v->delete_ts;
            }
            if (ets < sts && ets < dts)
            {
                if (debug)
                    (*logfile) << "Edge FOUND." << endl;
                return true;
            }
        }
        return false;
    }

    void initGraph(int n)
    {
        int i, j;
        for (i = 1; i <= n; i++)
        {
            AddV(i, nullptr, false, 0);
        }
        for (i = 1; i <= n; i = i + 2)
        {
            // for( j=i+1; j<=n; j=j+2){
            int u = rand() % n + 1;
            int v = rand() % n + 1;
            int weight = rand() % n + 1;
            if (u != v)
                AddE(u, v, weight, nullptr, false, 0);
        }
    }

    // return the appropriate version of a vertex
    v_version *find_my_version(v_version *V, int myTS, int tid)
    {
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        setTs(V, tid);
        v_version *curr_version = V;
        while (curr_version && curr_version->insert_ts.load() > myTS)
        {
            curr_version = curr_version->next_ver.load();
        }
        if (curr_version && (curr_version->delete_ts.load() > myTS || curr_version->delete_ts.load() == -1))
            return curr_version;
        else
            return nullptr;
    }

    // return the appropriate version of an edge // CHECK CHECK
    e_version *find_my_version(e_version *V, int myTS, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        setTs(V, tid);
        e_version *curr_version = V;
        while (curr_version && curr_version->insert_ts.load() > myTS)
        {
            curr_version = curr_version->next_ver.load();
        }
        if (curr_version && (curr_version->delete_ts.load() > myTS || curr_version->delete_ts.load() == -1) &&
            (curr_version->dest_version_node->delete_ts.load() > myTS || curr_version->dest_version_node->delete_ts.load() == -1))
            return curr_version;
        else
            return nullptr;
    }

    void helpResize(HNode *t, int i, int tid)
    {
        bucket *b = t->buckets_array[i];
        // vNode *b = t->buckets_array[i];
        HNode *s = t->old.load();
        if (b->get_bucket_head()->val == INT_MAX && s != nullptr)
        {
            vNode *set = nullptr;
            if (s->size * 2 == t->size) // growing
            {
                bucket *b_old = s->buckets_array[i % s->size];
                vNode *bucket_head_pointer = (vNode *)get_unfreezed_ref((long)b_old->get_bucket_head()); // verify get_unfreezed_ref correctness once
                bool ans = b_old->casHead(bucket_head_pointer, (vNode *)get_freezed_ref((long)bucket_head_pointer));
                set = b_old->split(t->size, i);
            }
            /*else // shrinking (WILL ONLY HAPPEN IF WE DO GARBAGE COLLECTION and remove the node from list)
            {
                bucket *b_old1 = s->buckets_array[i];
                bucket *b_old2 = s->buckets_array[i + t->size];
                vNode *bucket_head_pointer1 = (vNode *)get_unfreezed_ref((long)b_old1->get_bucket_head());
                vNode *bucket_head_pointer2 = (vNode *)get_unfreezed_ref((long)b_old2->get_bucket_head());
                b_old1->casHead(bucket_head_pointer1, (vNode *)get_freezed_ref((long)bucket_head_pointer1));
                b_old2->casHead(bucket_head_pointer2, (vNode *)get_freezed_ref((long)bucket_head_pointer2));
                set = b_old1->merge(b_old2);
            }*/
            // Assuming buckets is a std::atomic array or similar thread-safe structure (CHECK THIS)
            vNode *bucket_head = t->buckets_array[i]->get_bucket_head();
            if (bucket_head->val == INT_MAX)
            {
                t->buckets_array[i]->casHead(bucket_head, set);
            }
        }
    }

    bool hash_resize(HNode *t, bool grow, int tid)
    {
        if ((t->size == 1 << 16 && grow) || (t->size == 1 && !grow))
            return false;

        // cout<<"\n\n RESIZING the current head "<< t <<"\n\n";
        // for(int inx = 0; inx < t->size ; inx++)
        // {
        //     vNode* curr = hash_head.load()->buckets_array[inx]->get_bucket_head();
        //     while(curr){
        //         cout<<curr->val<<" ";
        //         curr = curr->vnext;
        //     }
        //     cout<<endl;
        // }

        if (t == hash_head.load())
        {
            // make sure we can deprecate t's predecessor
            for (int i = 0; i < t->size; i++)
            {
                if (t->buckets_array[i]->get_bucket_head()->val == INT_MAX)
                {
                    // cout<<"Calling helpResze from hash_resize function\n";
                    helpResize(t, i, tid);
                }
            }
            // deprecate t's predecessor
            t->old.store(nullptr);

            // cout<<"\n\n RESIZING the current head "<< t <<"\n\n";
            //  for(int inx = 0; inx < t->size ; inx++)
            //  {
            //      vNode* curr = hash_head.load()->buckets_array[inx]->get_bucket_head();
            //      while(curr){
            //          cout<<curr->val<<" ";
            //          curr = curr->vnext;
            //      }
            //      cout<<endl;
            //  }

            // switch to a new bucket array
            if (t == hash_head)
            {
                HNode *n = new HNode(t, grow ? t->size * 2 : t->size / 2);
                int log_Check = grow ? t->size * 2 : t->size / 2;
                // cout<< "Resized from " << t->size<<" to " << log_Check << endl;
                // return casHead(t, n);
                if (casHead(t, n))
                {
                    // cout << "Resized by " << tid << " from " << t->size <<"("<< t <<")"<< " to " << log_Check<<"("<< n <<")" << endl;
                    return true;
                }
            }
        }
        return false;
    }

    // will remove all the edge nodes and edge versions, since the corresponding v_version is completely obsolete
    // ayaz-GC-debug-comment-13 and 14 all 'outFileStream' lines
    //void clean_ALL_edges_for_me(v_version* curr_v_ver_temp, int tid, std::ofstream& outFileStream){
    void clean_ALL_edges_for_me(v_version* curr_v_ver_temp, int tid){
        //outFileStream<<"\n\n                        clean_ALL_edges_for_me\n"<<endl;
       auto guard1 = myRecManager_e_ver->getGuard(tid);
       auto guard2 = myRecManager_v_ver->getGuard(tid);
       auto guard3 = myRecManager_edge->getGuard(tid);
       eNode* curr_e_node = curr_v_ver_temp->enext.load(), *prev_e_node = curr_e_node, *curr_e_node_temp;
       curr_v_ver_temp->enext.store(nullptr);
       while(curr_e_node){
           curr_e_node_temp = curr_e_node;
           curr_e_node = curr_e_node->enext.load();
           // retire all the edge versions
           //outFileStream << "CHECK-1: Working on edge node with dest: "<< curr_e_node_temp->val<< " & address : " << curr_e_node_temp  <<endl;
           e_version* curr_e_ver = curr_e_node_temp->vhead.load(), *curr_e_ver_temp;
           curr_e_node_temp->vhead.store(nullptr);
           while(curr_e_ver){
               curr_e_ver_temp = curr_e_ver;
               curr_e_ver = curr_e_ver->next_ver.load();
               //curr_e_ver_temp->next_ver.store(nullptr);
              //outFileStream << "       -> CHECK-2: Retiring edge_version with address : " << curr_e_ver_temp  <<endl;
               myRecManager_e_ver->template retire(tid, curr_e_ver_temp);
           }
           //outFileStream << "CHECK-3: Retiring edge NODE with address : " << curr_e_node_temp  <<endl;
           myRecManager_edge->template retire(tid, curr_e_node_temp);
       }
       //curr_v_ver_temp->enext.store(nullptr);
   }

    // will remove ONLY some of the edge versions, and none of the edge nodes
   // ayaz-GC-debug-comment-12 and 13 all 'outFileStream' lines
   // void clean_SOME_edges_for_me(v_version* curr_v_ver_temp, int tid, int my_valid_GC_timestamp, std::ofstream& outFileStream)
//    void clean_SOME_edges_for_me(v_version* curr_v_ver_temp, int tid, int my_valid_GC_timestamp){
//         //outFileStream<<"\n\n                    clean_SOME_edges_for_me\n"<<endl;
//        auto guard1 = myRecManager_e_ver->getGuard(tid);
//        auto guard2 = myRecManager_v_ver->getGuard(tid);
//        auto guard3 = myRecManager_edge->getGuard(tid);
//        eNode* curr_e_node = curr_v_ver_temp->enext.load(), *curr_e_node_temp;
//        curr_e_node = curr_e_node->enext.load();
//        while(curr_e_node->val != INT_MAX){
//             //outFileStream << "Point-1: Edge Node value : "<< curr_e_node->val<< " address : " << curr_e_node   <<endl;
//             curr_e_node_temp = curr_e_node;
//             curr_e_node = curr_e_node->enext.load();
//             // retire few of the edge versions
//             e_version* curr_e_ver = curr_e_node_temp->vhead.load(), *prev_e_ver = curr_e_ver, *curr_e_ver_temp = curr_e_ver;
//            // CASE: where we can't free the first e_node version, so we will set the delete timestamp if its destination vertex is gonna get removed
//                int del_ts_of_dest_ver = curr_e_ver_temp->dest_version_node->delete_ts ;
//                if(del_ts_of_dest_ver == -1){
//                     setTs(curr_e_ver_temp->dest_version_node, tid);
//                     del_ts_of_dest_ver = curr_e_ver_temp->dest_version_node->delete_ts ;
//                     //outFileStream << "Point-2: Set destination vertex "<< curr_e_ver_temp->dest_version_node << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
//                }
//                if((del_ts_of_dest_ver < my_valid_GC_timestamp)){
//                     int abc = curr_e_ver_temp->delete_ts.load();
//                     while(!curr_e_ver_temp->delete_ts.compare_exchange_strong(abc, del_ts_of_dest_ver)){
//                         abc = curr_e_ver_temp->delete_ts.load();
//                         //outFileStream << "Point-3: Retrying to set edge version "<< curr_e_ver_temp << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
//                     }
//                     //outFileStream << "Point-4: edge version "<< curr_e_ver_temp << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
//                }
//             //prev_e_ver = curr_e_ver_temp;
//             bool flag = 1;   // to hold the prev pointer at the same place
//             curr_e_ver = curr_e_ver->next_ver.load();
//             while(curr_e_ver){
//                 if(flag)
//                     prev_e_ver = curr_e_ver_temp;
//                 flag = 1;
//                 curr_e_ver_temp = curr_e_ver;
//                 curr_e_ver = curr_e_ver->next_ver.load();
//                 int del_ts_of_dest_ver_2 = curr_e_ver_temp->dest_version_node->delete_ts ;
//                 if(del_ts_of_dest_ver_2 == -1){
//                     setTs(curr_e_ver_temp->dest_version_node, tid);
//                     del_ts_of_dest_ver_2 = curr_e_ver_temp->dest_version_node->delete_ts ;
//                     //outFileStream << "Point-5: Helped setting vertex version "<< curr_e_ver_temp->dest_version_node << " del timestamp as" <<  del_ts_of_dest_ver_2 <<endl;
//                 }
//                 if(((del_ts_of_dest_ver_2 < my_valid_GC_timestamp)) ){
//                     curr_e_ver_temp->dest_version_node = nullptr;       /// ONLY CHANGE @ 14:24
//                     prev_e_ver->next_ver.store(curr_e_ver); // problem if more than 1 thread does GC
//                     //curr
//                     //outFileStream << "Point-6: Set e_ver "<< prev_e_ver << " next version as e_ver " <<  curr_e_ver <<endl;
//                     //curr_e_ver_temp->next_ver.store(nullptr);
//                     //outFileStream << "Point-7: Retiring e_ver"<< curr_e_ver_temp <<endl;
//                     myRecManager_e_ver->template retire(tid, curr_e_ver_temp);
//                     flag = 0;
//                 }
//            }
//            //myRecManager_edge->template retire(tid, curr_e_node_temp);
//        }
//    }


    // NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW APPROACH (not freeing up edge versions)
    // will remove ONLY some of the edge versions, and none of the edge nodes
    // ayaz-GC-debug-comment-12 and 13 all 'outFileStream' lines
    //void clean_SOME_edges_for_me(v_version* curr_v_ver_temp, int tid, int my_valid_GC_timestamp, std::ofstream& outFileStream){
    void clean_SOME_edges_for_me(v_version* curr_v_ver_temp, int tid, int my_valid_GC_timestamp){
        //outFileStream<<"\n\n                    clean_SOME_edges_for_me\n"<<endl;
       auto guard1 = myRecManager_e_ver->getGuard(tid);
       auto guard2 = myRecManager_v_ver->getGuard(tid);
       auto guard3 = myRecManager_edge->getGuard(tid);
       eNode* curr_e_node = curr_v_ver_temp->enext.load(), *curr_e_node_temp;
       curr_e_node = curr_e_node->enext.load();

       while(curr_e_node->val != INT_MAX){
            //outFileStream << "Point-1: Edge Node value : "<< curr_e_node->val<< " address : " << curr_e_node   <<endl;
            curr_e_node_temp = curr_e_node;
            curr_e_node = curr_e_node->enext.load();
            // retire few of the edge versions
            e_version* curr_e_ver = curr_e_node_temp->vhead.load(), *prev_e_ver = curr_e_ver, *curr_e_ver_temp = curr_e_ver;
           
           // ***********************************************************************************************************
           // CASE: where we can't free the first e_node version even if it should be, so we will set its delete timestamp, if its destination vertex is gonna get removed
           // SO HANDLE FIRST E_VERSION_NODE DIFFERENTLY ONLY    
           // ***********************************************************************************************************
               int del_ts_of_dest_ver = curr_e_ver_temp->dest_version_node->delete_ts ;
               if(del_ts_of_dest_ver == -1){
                    // help even though no action needed
                    setTs(curr_e_ver_temp->dest_version_node, tid);
                    del_ts_of_dest_ver = curr_e_ver_temp->dest_version_node->delete_ts ;
                    //outFileStream << "Point-2: Set destination vertex "<< curr_e_ver_temp->dest_version_node << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
               }
               if((del_ts_of_dest_ver < my_valid_GC_timestamp)){
                    int abc = curr_e_ver_temp->delete_ts.load();
                    while(!curr_e_ver_temp->delete_ts.compare_exchange_strong(abc, del_ts_of_dest_ver)){
                        abc = curr_e_ver_temp->delete_ts.load();
                        //outFileStream << "Point-3: Retrying to set edge version "<< curr_e_ver_temp << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
                    }
                    //outFileStream << "Point-4: edge version "<< curr_e_ver_temp << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
                    curr_e_ver = curr_e_ver->next_ver.load();
                    curr_e_ver_temp->next_ver.store(nullptr);
                    // retire all other edge_Versions
                    while(curr_e_ver){
                        curr_e_ver_temp = curr_e_ver;
                        curr_e_ver = curr_e_ver->next_ver.load();
                        myRecManager_e_ver->template retire(tid, curr_e_ver_temp);
                    }
                    continue;
               }
               else{
                // move to next version
                    curr_e_ver = curr_e_ver->next_ver.load();
               }
            bool flag = 1;   // to hold the prev pointer at the same place
            while(curr_e_ver){
                curr_e_ver_temp = curr_e_ver;
                curr_e_ver = curr_e_ver->next_ver.load();
                del_ts_of_dest_ver = curr_e_ver_temp->dest_version_node->delete_ts ;
                if(del_ts_of_dest_ver == -1){
                    setTs(curr_e_ver_temp->dest_version_node, tid);
                    del_ts_of_dest_ver = curr_e_ver_temp->dest_version_node->delete_ts ;
                    //outFileStream << "Point-2: Set destination vertex "<< curr_e_ver_temp->dest_version_node << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
                }
                if((del_ts_of_dest_ver < my_valid_GC_timestamp) || (curr_e_ver_temp->delete_ts < my_valid_GC_timestamp)){
                    
                    if(del_ts_of_dest_ver < my_valid_GC_timestamp){
                        int abc = curr_e_ver_temp->delete_ts.load();
                        while(!curr_e_ver_temp->delete_ts.compare_exchange_strong(abc, del_ts_of_dest_ver)){
                            abc = curr_e_ver_temp->delete_ts.load();
                            //outFileStream << "Point-3: Retrying to set edge version "<< curr_e_ver_temp << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
                        }
                    }
                    prev_e_ver->next_ver.store(nullptr);
                    curr_e_ver_temp->dest_version_node = nullptr;
                    myRecManager_e_ver->template retire(tid, curr_e_ver_temp);
                    curr_e_ver = curr_e_ver->next_ver.load();
                    while(curr_e_ver){
                        curr_e_ver_temp = curr_e_ver;
                        curr_e_ver = curr_e_ver->next_ver.load();
                        myRecManager_e_ver->template retire(tid, curr_e_ver_temp);
                    }
                    //outFileStream << "Point-4: edge version "<< curr_e_ver_temp << " del timestamp as" <<  del_ts_of_dest_ver <<endl;
                }
                prev_e_ver = curr_e_ver_temp; // should not reach here if curr_e_ver_temp has been retired above
           }
       }
   }




    // // GC routine for vertex versions (GET IT CHECKED BY """""SIR""""") (add removal of edge nodes as well as edge version nodes)
    // void do_GC_vertex_ver(v_version *curr_v_version, int my_valid_GC_timestamp, fstream *logfile, bool debug, int tid){
    //     auto guard1 = myRecManager_e_ver->getGuard(tid);
    //     auto guard2 = myRecManager_v_ver->getGuard(tid);
    //     auto guard3 = myRecManager_edge->getGuard(tid);
    //     //myRecManager_v_ver
    //     v_version *prev_v_version = curr_v_version;
    //     curr_v_version = curr_v_version->next_ver.load();
    //     v_version *curr_v_ver_temp = curr_v_version;
    //     while(curr_v_version)
    //     {
    //         curr_v_ver_temp = curr_v_version;
    //         curr_v_version = curr_v_version->next_ver.load();
    //         if((curr_v_ver_temp->delete_ts != -1) && (curr_v_ver_temp->delete_ts < my_valid_GC_timestamp))
    //         {
    //             prev_v_version->next_ver.store(nullptr);
    //             // retire all the edges
    //             clean_ALL_edges_for_me(curr_v_ver_temp, tid);
    //             myRecManager_v_ver->template retire(tid, curr_v_ver_temp);
    //             while(curr_v_version)     // to avoid the IF condtion from here on
    //             {
    //                 curr_v_ver_temp = curr_v_version;
    //                 prev_v_version = curr_v_version;
    //                 curr_v_version = curr_v_version->next_ver.load();
    //                 prev_v_version->next_ver.store(nullptr);
    //                 clean_ALL_edges_for_me(curr_v_ver_temp, tid);
    //                 myRecManager_v_ver->template retire(tid, curr_v_ver_temp);
    //             }
    //         }
    //     }
    // }

    void do_GC(fstream *logfile, bool debug, int tid, int graph_gc_cnt){
        int my_valid_GC_timestamp = find_min_TS(logfile, debug, tid);
        HNode *start_hash_head;
        start_hash_head = hash_head.load();
        //cout<<"hi from GC TS "<< my_valid_GC_timestamp<<endl;
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        auto guard3 = myRecManager_edge->getGuard(tid);

        //ayaz-GC-debug-comment-1
        cout<<"GC-number-"<<graph_gc_cnt<<" with GC_ts as = "<<my_valid_GC_timestamp<<endl;
        // print the entire graph (logging)
        /*
        string file_name = "/home/student/2022/ayaz/test_dec_18/log_ayaz/graph_status_before_GC_"+to_string(graph_gc_cnt)+".txt";
        ofstream MyFile(file_name);
        HNode *log_t;
        log_t = hash_head.load();
        vNode* log_curr_v;
        v_version *log_curr_v_version = nullptr;
        for(int i = 0;  i<log_t->size; i++) {
            MyFile<<"\n\n\n\n***************************************************\n **************BUCKET_NUMBER "<<i<<"************** \n***************************************************\n\n\n"<<endl;
            log_curr_v = (vNode*)get_unfreezed_ref((long)log_t->buckets_array[i]->get_bucket_head());

            while(log_curr_v != nullptr){
                MyFile<<"\n\n ##### Vertex Node Address "<<log_curr_v <<" Vertex Value "<<log_curr_v->val<<"##### \n"<<endl;
                log_curr_v_version = log_curr_v->vhead;
                if(log_curr_v_version == nullptr){  // should never enter
                    log_curr_v = log_curr_v->vnext;
                    continue;
                }
                int v_version_cnt = 1;
                while(log_curr_v_version){
                    MyFile<<"\n        $$$$ Vertex version - "<<v_version_cnt++ <<" Address "<<log_curr_v_version<<" with value "<< log_curr_v_version->val<<" ins-del "<<log_curr_v_version->insert_ts<<"-"<<log_curr_v_version->delete_ts <<"$$$$ \n"<<endl;
                    eNode* log_curr_e_node = log_curr_v_version->enext.load();
                    while(log_curr_e_node){
                        MyFile<<"       ~~~~~~ Edge Node addr - "<<log_curr_e_node <<" with dest vertex value as "<<log_curr_e_node->val<<" ~~~~~~ "<<endl;
                        e_version* log_curr_e_ver = log_curr_e_node->vhead.load();
                        while(log_curr_e_ver){
                            MyFile<<"           &&&&&& Edge version addr- "<<log_curr_e_ver <<" with dest vertex version addr as "<<log_curr_e_ver->dest_version_node <<" ins-del "<<log_curr_e_ver->insert_ts<<"-"<<log_curr_e_ver->delete_ts<<" &&&&&&  "<<endl;
                            log_curr_e_ver = log_curr_e_ver->next_ver.load();
                        }
                        log_curr_e_node = log_curr_e_node->enext;
                    }
                    log_curr_v_version = log_curr_v_version->next_ver;
                }
                log_curr_v = log_curr_v->vnext;
            }
        }
        MyFile.close();

        file_name = "/home/student/2022/ayaz/test_dec_18/log_ayaz_2/GC_log_"+to_string(graph_gc_cnt)+".txt";
        ofstream MyFile_2(file_name);
        
        
        if(!MyFile_2.is_open()){
            cout<<"coudln't open my file 2"<<endl;
        }
        */

        vNode* curr_v;
        v_version *curr_v_version = nullptr;
        v_version *prev_v_version = nullptr;

        HNode *t;
        int vertex_counter = 0, version_counter = 0;
        // PART-1: Cleaning up all edge nodes and edge versions
        t = hash_head.load();
        for(int i = 0;  i<t->size; i++){
            //ayaz-GC-debug-comment-2
            //MyFile_2<<"\n\n             Edge node and version cleanup, Working on Bucket "<<i<<"\n"<<endl;
            curr_v = t->buckets_array[i]->get_bucket_head();
            while(is_freezed_ref((long)curr_v)){
                t = hash_head.load();
                curr_v = t->buckets_array[i]->get_bucket_head();
            }
            if (curr_v->val == INT_MAX){
                helpResize(t, i, tid);
                curr_v =(vNode*)get_unfreezed_ref((long)t->buckets_array[i]->get_bucket_head());
            }

            // PART-1a: iterating all the vertices of this bucket
            while(curr_v != nullptr && curr_v->val != INT_MAX)
            {
                //MyFile_2<<" Bucket "<<i<<"'s vertex node address "<<curr_v<<" and value = "<<curr_v->val<<"is being cleaned"<<endl;
                vertex_counter++; // for debug purpose
                curr_v_version = curr_v->vhead;
                prev_v_version = nullptr;
                if(curr_v_version == nullptr){
                    curr_v = curr_v->vnext;
                    continue;
                }
                //curr_v_version = curr_v_version->next_ver.load();
                v_version *curr_v_ver_temp = curr_v_version;
                
                while(curr_v_version)
                {
                    //MyFile_2<<" ~~~~~~ Vertex node version address "<<curr_v_version<<"is being cleaned"<<endl;
                    version_counter++;  // for debug purpose
                    //cout<<"ayaz "<<my_valid_GC_timestamp<<endl;
                    curr_v_ver_temp = curr_v_version;
                    curr_v_version = curr_v_version->next_ver.load();
                    // CASE-1: where we will be cleaning only some of the edge nodes
                    int del_TS_v_ver = curr_v_ver_temp->delete_ts;
                    if( ((del_TS_v_ver == -1) || (del_TS_v_ver >= my_valid_GC_timestamp))){
                        //prev_v_version->next_ver.store(nullptr);
                        //cout<<"ayaz-2 "<<vertex_//            while(is_freezed_ref((long)curr_v)){
//                t = hash_head.load();
//                curr_v = t->buckets_array[i]->get_bucket_head();
//            }counter<<" "<<version_counter<<endl;
                        setTs(curr_v_ver_temp, tid);
                        // retire some of the edges
                        // ayaz-GC-debug-comment-3
                        //MyFile_2<<"Calling CLEAN_SOME_edges_for_me(from PARTIAL part) for vertex_version = "<<curr_v_ver_temp<<endl;
                        //clean_SOME_edges_for_me(curr_v_ver_temp, tid, my_valid_GC_timestamp, MyFile_2);
                        clean_SOME_edges_for_me(curr_v_ver_temp, tid, my_valid_GC_timestamp);
                        while(curr_v_version)     // to avoid the IF condition from here on
                        {
                            curr_v_ver_temp = curr_v_version;
                            curr_v_version = curr_v_version->next_ver.load();
                            // ayaz-GC-debug-comment-4
                            //MyFile_2<<"calling CLEAN_ALL_edges_for_me(from PARTIAL part)  for vertex_version = "<<curr_v_ver_temp<<endl;
                            //clean_ALL_edges_for_me(curr_v_ver_temp, tid, MyFile_2);
                            clean_ALL_edges_for_me(curr_v_ver_temp, tid);
                        }
                    }
                    // CASE-2: where we will be cleaning "ALL" of the edge nodes and their edge versions
                    else if((del_TS_v_ver != -1) && (del_TS_v_ver < my_valid_GC_timestamp)){
                        // retire all the edges
                        //cout<<"ayaz-3 "<<vertex_counter<<" "<<version_counter<<endl;
                        // ayaz-GC-debug-comment-5
                        //MyFile_2<<"calling (a)CLEAN_ALL_edges_for_me(from ALL part)  for vertex_version = "<<curr_v_ver_temp<<endl;
                        //clean_ALL_edges_for_me(curr_v_ver_temp, tid, MyFile_2);
                        clean_ALL_edges_for_me(curr_v_ver_temp, tid);
                        while(curr_v_version)     // to avoid the IF condition from here on
                        {
                            curr_v_ver_temp = curr_v_version;
                            //prev_v_version = curr_v_version;
                            curr_v_version = curr_v_version->next_ver.load();
                            // ayaz-GC-debug-comment-6
                            //MyFile_2<<"calling (b)CLEAN_ALL_edges_for_me(from ALL part)  for vertex_version = "<<curr_v_ver_temp<<endl;
                            //clean_ALL_edges_for_me(curr_v_ver_temp, tid, MyFile_2);
                            clean_ALL_edges_for_me(curr_v_ver_temp, tid);
                        }
                    }
                }
                curr_v = curr_v->vnext;
            }
        }
        t = hash_head.load();
        int size_hash = t->size;
        // PART-2: Cleaning up all vertex versions
        for(int i = 0; i < size_hash; i++){
            // ayaz-GC-debug-comment-7
            //MyFile_2<<"\n\n\n\n\n\n\n             Vertex version cleanup, Working on Bucket "<<i<<"\n"<<endl;
            //t = hash_head.load();
            curr_v = t->buckets_array[i]->get_bucket_head();
            while(is_freezed_ref((long)curr_v)){
                t = hash_head.load();
                curr_v = t->buckets_array[i]->get_bucket_head();
            }
            if (curr_v->val == INT_MAX){
                    helpResize(t, i, tid);
                    curr_v =(vNode*)get_unfreezed_ref((long)t->buckets_array[i]->get_bucket_head());
            }

            while(curr_v != nullptr && curr_v->val != INT_MAX)
            {
                // ayaz-GC-debug-comment-8
                //MyFile_2<<"\n\n     Vertex version cleanup, Working on vertex "<<curr_v<<"\n"<<endl;
                curr_v_version = curr_v->vhead;
                prev_v_version = curr_v_version;

                if(curr_v_version == nullptr){
                    curr_v = curr_v->vnext;
                    continue;
                }
                
                v_version *curr_v_ver_temp = curr_v_version;
                curr_v_version = curr_v_version->next_ver.load();
                
                
                while(curr_v_version)
                {
                    prev_v_version = curr_v_ver_temp;
                    curr_v_ver_temp = curr_v_version;
                    curr_v_version = curr_v_version->next_ver.load();
                    int del_TS_v_ver = curr_v_ver_temp->delete_ts;
                    if((del_TS_v_ver != -1) && (del_TS_v_ver < my_valid_GC_timestamp)){
                        prev_v_version->next_ver.store(nullptr);
                        // ayaz-GC-debug-comment-9
                        //MyFile_2<<"\n\n                 -> Vertex version(a) cleanup, Working on vertex version "<<curr_v_ver_temp<<"\n"<<endl;
                        myRecManager_v_ver->template retire(tid, curr_v_ver_temp);
                        while(curr_v_version)     // to avoid the IF condtion from here on
                        {

//                            prev_v_version = curr_v_ver_temp;
                            curr_v_ver_temp = curr_v_version;
                            //prev_v_version = curr_v_version;
                            curr_v_version = curr_v_version->next_ver.load();
//                            prev_v_version->next_ver.store(nullptr);
                            // ayaz-GC-debug-comment-10
                            //MyFile_2<<"\n\n                 -> Vertex version(b) cleanup, retiring vertex version "<<curr_v_ver_temp<<"\n"<<endl;
                           myRecManager_v_ver->template retire(tid, curr_v_ver_temp);
                        }
                    }
                }
                curr_v = curr_v->vnext;
            }
            // ayaz-GC-debug-comment-11

        }
        //MyFile_2.close();
    }

    //// does the real operation of adding a new vertex and removing it
    int invoke(int key, bool is_insert_op, bool *failed_bcz_head_freezed, fstream *logfile, bool debug, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        auto guard3 = myRecManager_edge->getGuard(tid);
        // vector<vNode*> test_vector;
        // vector<int> test_vector_log;
        int my_valid_GC_timestamp = find_min_TS(logfile, debug, tid);

        int itr_number = 0;
        int nos_of_elem_in_bucket = 0;
        HNode *h = hash_head.load();
        bucket *b = h->buckets_array[key % h->size];
        if (b->get_bucket_head()->val == INT_MAX)
        {
            // cout<<"Calling helpResze from invoke function by key "<<key <<"\n";
            helpResize(h, key % h->size, tid);
        }
        vNode *bucket_head = b->get_bucket_head(); // bucket's head
        int marker = INT_MAX;                      // stopper for preventing complete iteration
        while (!is_freezed_ref((long)bucket_head))
        {
            // Step-1: search for the key until you reach the stop point, else keep trying until you successfully add your node.
            // from here on start integrating LOCATE V etc.....
            vNode *curr = b->get_bucket_head();
            vNode *curr_copy = curr;
            while (curr->val != marker) // iterate until the marker point
            {
                // test_vector.push_back(curr);
                // test_vector_log.push_back(itr_number);
                nos_of_elem_in_bucket++;
                if (curr->val == key)
                    break;
                curr = curr->vnext;
            }
            marker = curr_copy->val; // have to stop searching here the next time

            if (curr && curr->val == key) // if found the key, check the operation (no need to REITERATE the while loop)
            {
                if (is_insert_op)
                {
                    v_version *curr_version = read(curr, tid);
                    if (curr_version->delete_ts == INT_MAX)
                    {
                        if (debug)
                            (*logfile) << "Vertex exists with key : " << key << "(" << curr_version << ")" << endl;
                        //do_GC_vertex_ver(curr_version, my_valid_GC_timestamp, logfile, debug, tid);
                        return -nos_of_elem_in_bucket;
                    }
                    else
                    {
                        setTs(curr_version, tid);
                        v_version *new_version = create_vertex_Ver(key, tid);
                        new_version->next_ver.store(curr_version);
                        if (curr->vhead.compare_exchange_strong(curr_version, new_version))
                        {
                            setTs(new_version, tid);
                            if (debug)
                                (*logfile) << "Vertex added with key : " << key << "(" << new_version << ")" << endl;
                            //do_GC_vertex_ver(new_version, my_valid_GC_timestamp, logfile, debug, tid);
                            return nos_of_elem_in_bucket;
                        }
                        else
                        {
                            if (debug)
                                (*logfile) << "Vertex added by other thread, therefore thread CAS failed." << endl;
                            return -nos_of_elem_in_bucket;
                        }
                    }
                }
                else
                {
                    v_version *curr_version = read(curr, tid);
                    int t = INT_MAX;
                    if (curr_version->delete_ts != INT_MAX)
                    {
                        setTs(curr_version, tid);
                        if (debug)
                            (*logfile) << "Vertex ALREADY deleted." << endl;
                        //do_GC_vertex_ver(curr_version, my_valid_GC_timestamp, logfile, debug, tid);
                        return -nos_of_elem_in_bucket;
                    }
                    else if (curr_version->delete_ts.compare_exchange_strong(t, -1))
                    {
                        setTs(curr_version, tid);
                        if (debug)
                            (*logfile) << "Vertex successfully deleted." << endl;
                        //do_GC_vertex_ver(curr_version, my_valid_GC_timestamp, logfile, debug, tid);
                        return nos_of_elem_in_bucket;
                    }
                    else
                    {
                        if (debug)
                            (*logfile) << "Vertex CONCURRENTLY deleted by other thread." << endl;
                        return -nos_of_elem_in_bucket; // deleted by some other thread
                    }
                }
            }
            else if (curr && curr->val != key) // this case will only arise if the node never existed (there may or maynot have been a try to insert the node)
            {
                if (is_insert_op)
                {
                    vNode *newv = createV(key, tid); // create a new vertex node
                    newv->vnext = curr_copy;
                    if (b->casHead(curr_copy, newv))
                    {
                        setTs(newv->vhead, tid);
                        if (debug)
                            (*logfile) << key << " Vertex ADDED." << endl;
                        return nos_of_elem_in_bucket;
                    }
                }
                else
                {
                    return -nos_of_elem_in_bucket; // key not present
                }
            }
            else
            {
                if (is_insert_op)
                {
                    vNode *newv = createV(key, tid); // create a new vertex node
                    newv->vnext = curr_copy;
                    if (b->casHead(curr_copy, newv))
                    {
                        setTs(newv->vhead, tid);
                        if (debug)
                            (*logfile) << key << " Vertex ADDED." << endl;
                        return nos_of_elem_in_bucket;
                    }
                }
                else
                {
                    return -nos_of_elem_in_bucket; // key not present
                }
            }
            // bucket_head = head.load();
            bucket_head = b->get_bucket_head();
            itr_number++;
        }
        *failed_bcz_head_freezed = true;
        return -nos_of_elem_in_bucket;
    }

    // add a new vertex in the vertex-list
    bool AddV(int key, fstream *logfile, bool debug, int tid)
    {
        bool failed_bcz_head_freezed;
        int result;
        HNode *h;
        // remove start
        h = hash_head.load();
        //cout<<"Hash size "<<h->size<<endl;
        // remove end
        do
        {
            h = hash_head.load();
            failed_bcz_head_freezed = false;
            result = invoke(key, true, &failed_bcz_head_freezed, logfile, debug, tid);
            if (failed_bcz_head_freezed)
            {
                hash_resize(h, true, tid);
            }
        } while (failed_bcz_head_freezed);

      if (std::abs(result) > 8)
      {
          //cout<<"********************************************\n";
          hash_resize(h, true, tid);
      }
        return result > 0;
    }

    bool AddV2(int key, fstream *logfile, bool debug, int tid)
    {
        bool failed_bcz_head_freezed;
        int result;
        HNode *h;
        do
        {
            h = hash_head.load();
            failed_bcz_head_freezed = false;
            result = invoke(key, true, &failed_bcz_head_freezed, logfile, debug, tid);
            if (failed_bcz_head_freezed)
            {
                hash_resize(h, true, tid);
            }
        } while (failed_bcz_head_freezed);

        if (std::abs(result) > 8)
        {
            hash_resize(h, true, tid);
        }
        return result > 0;
    }

    // Deletes the vertex from the vertex-list
    bool RemoveV(int key, fstream *logfile, bool debug, int tid)
    {
        bool failed_bcz_head_freezed;
        int result;
        HNode *h;
        do
        {
            h = hash_head.load();
            failed_bcz_head_freezed = false;
            result = invoke(key, false, &failed_bcz_head_freezed, logfile, debug, tid);
            if (failed_bcz_head_freezed)
            {
                hash_resize(h, true, tid);
            }
        } while (failed_bcz_head_freezed);

    //   if (std::abs(result) > 8)
    //   {
    //       //cout<<"********************************************\n";
    //       hash_resize(h, true, tid);
    //   }
        return result > 0;
    }

    // Deletes the vertex from the vertex-list
    bool RemoveV2(int key, fstream *logfile, bool debug, int tid)
    {
        bool failed_bcz_head_freezed;
        int result;
        HNode *h;
        do
        {
            h = hash_head.load();
            failed_bcz_head_freezed = false;
            result = invoke(key, false, &failed_bcz_head_freezed, logfile, debug, tid);
            if (failed_bcz_head_freezed)
            {
                hash_resize(h, true, tid);
            }
        } while (failed_bcz_head_freezed);

        if (std::abs(result) > 8)
        {
            hash_resize(h, true, tid);
        }
        return result > 0;
    }

    // partial
    BFS_node *getBFS(int src, fstream *logfile, bool debug, int tid)
    {
        
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        auto guard3 = myRecManager_edge->getGuard(tid);

        shared_active_version_history[tid].store(-1);
        int my_TS = getMyTs(); // get the timestamp and increment the global timestamp
        int expected_val = -1;
        if(!shared_active_version_history[tid].compare_exchange_strong(expected_val, my_TS)){
            my_TS = shared_active_version_history[tid];
        }
        shared_active_version_history[tid].store(my_TS);
        
        BFS_node *head = nullptr; // to track back the path
        v_version *u, *v;
        bool flag = ContainsV(src, &u, logfile, debug, tid);
        if (!flag)
        {
            if (debug)
                (*logfile) << "Returning from Getpath because u or v not found" << endl;
            return nullptr;
        }
        u = find_my_version(u, my_TS, tid);
        if (u == nullptr) // no appropriate version found
        {
            if (debug)
                (*logfile) << "Returning from Getpath because no appropriate version found" << endl;
            return nullptr;
        }
        queue<Queue_entry> Q;
        unordered_set<int> visited;
        head = new BFS_node(u->val, nullptr);
        Queue_entry temp_Q_entry(u, head);
        Q.push(temp_Q_entry);
        if (debug)
            (*logfile) << "pushing " << temp_Q_entry.vertex_bfs_node->val << " and Q.size = " << Q.size() << endl;
        visited.insert(u->val);
        while (!Q.empty())
        {
            Queue_entry top(Q.front().vertex_version, Q.front().vertex_bfs_node);
            Q.pop();
            if (debug)
            {
                (*logfile) << "Poping " << top.vertex_bfs_node->val << " and Q.size = " << Q.size() << endl;
                (*logfile) << top.vertex_bfs_node->val << " edges are as follows : ";

                eNode *temp_curr_edge = top.vertex_version->enext.load();
                do{
                    (*logfile) << temp_curr_edge->val << " ";
                    temp_curr_edge = temp_curr_edge->enext.load();
                } while (temp_curr_edge);
                (*logfile) << endl;
                (*logfile) << "Entering while loop" << endl;
            }
            //eNode *curr_edge = top.vertex_version->enext.load()->enext.load();

            eNode *curr_edge = top.vertex_version->enext.load();
            if(curr_edge == nullptr){
                continue;
            }
            curr_edge = curr_edge->enext.load();

            while (curr_edge && curr_edge->val != INT_MAX)
            {
                if (debug)
                    (*logfile) << curr_edge->val << endl;
                e_version *live_version = find_my_version(curr_edge->vhead, my_TS, tid);

                if (live_version != nullptr)
                {
                    if (debug)
                        (*logfile) << "version found with insert timestamp as " << live_version->insert_ts << "delete timestamp as " << live_version->delete_ts << " and my timestamp as " << my_TS << endl;
                    v_version *destination_version = live_version->dest_version_node;

                    if (destination_version && visited.find(destination_version->val) == visited.end()) // if not already visited
                    {
                        if (debug)
                        {
                            (*logfile) << "Entered IF" << endl;
                            (*logfile) << "Dest_Version node " << destination_version << endl;
                        }
                        temp_Q_entry.vertex_version = live_version->dest_version_node;
                        BFS_node *new_bfs_node = new BFS_node(destination_version->val, top.vertex_bfs_node, top.vertex_bfs_node->child_head);
                        top.vertex_bfs_node->child_head = new_bfs_node;
                        temp_Q_entry.vertex_bfs_node = new_bfs_node;
                        Q.push(temp_Q_entry);
                        if (debug)
                            (*logfile) << "pushing " << temp_Q_entry.vertex_bfs_node->val << " and Q.size = " << Q.size() << endl;
                        visited.insert(destination_version->val);
                    }
                }
                curr_edge = curr_edge->enext.load();
            }
        }
        return head; // can easily back track the path if not NULL
    }

    // full snapshot
    unordered_map<int, std::list<std::pair<int, int>>> snapshot(fstream *logfile, bool debug, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        auto guard3 = myRecManager_edge->getGuard(tid);

        shared_active_version_history[tid].store(-1);
        int my_TS = getMyTs(); // get the timestamp and increment the global timestamp
        int expected_val = -1;
        if(!shared_active_version_history[tid].compare_exchange_strong(expected_val, my_TS)){
            my_TS = shared_active_version_history[tid];
        }
        shared_active_version_history[tid].store(my_TS);

        unordered_map<int, list<pair<int, int>>> my_graph;
        HNode *h = hash_head.load();
        if (h == nullptr)
            return my_graph;

        for (int i = 0; i < h->size; i++)
        {
            //cout<<"HIHIHI-"<<tid<<"\n";
            // get correct bucket to collect from
            vNode *bucket_head = h->buckets_array[i]->get_bucket_head();
            while(is_freezed_ref((long)bucket_head))
            {
                h = hash_head.load();
                bucket_head = h->buckets_array[i]->get_bucket_head();
                //cout<<"HIHIHI\n";
            }
            // collect this bucket's vertices
            while (bucket_head->val != INT_MAX)
            {
                //cout<<"Inside While loop \n";
                if (my_graph.find(bucket_head->val) == my_graph.end())
                {
                    v_version *curr_version_vertex = read(bucket_head, tid);                               // read the cuurent version head
                    v_version *my_ver_vertex = find_my_version(curr_version_vertex, my_TS, tid); // now find the appropriate version as per my_TS time
                    if (my_ver_vertex == nullptr)
                    {
                        bucket_head = bucket_head->vnext;
                        continue;
                    }
                    // iterate its edges
                    list<pair<int, int>> edge_list;
                    //eNode *e = my_ver_vertex->enext.load()->enext.load(); // initialised to head sentinel of orig vertex node (edge list is already filled) (my_Ver does not has head, therefore TWO re-direction needed)
                    
                    eNode *e = my_ver_vertex->enext.load();
                    if(e == nullptr){
                        continue;
                    }
                    e = e->enext.load();
                    
                    
                    while (e->val != INT_MAX)
                    {
                        e_version *curr_version_edge = read(e, tid);                             // read the current version head
                        e_version *my_ver_edge = find_my_version(curr_version_edge, my_TS, tid); // now find the appropriate version as per my_TS time
                        if (my_ver_edge == nullptr)
                        {
                            e = e->enext.load();
                            continue;
                        }
                        edge_list.push_back(make_pair(my_ver_edge->dest_version_node->val, my_ver_edge->weight));
                        e = e->enext.load();
                    }

                    my_graph.insert(make_pair(bucket_head->val, edge_list));
                    bucket_head = bucket_head->vnext;
                }
                else{       // no need to explore rest of the elements in this bucket
                    break;
                }
            }
        }
        shared_active_version_history[tid].store(INT_MAX);
        return my_graph;
    }


    // Dijkstra's also for single source shortest path
    vector<pair<long, bool>> getSSSP(int source, fstream *logfile, bool debug, int tid)
    {
        auto guard1 = myRecManager_e_ver->getGuard(tid);
        auto guard2 = myRecManager_v_ver->getGuard(tid);
        auto guard3 = myRecManager_edge->getGuard(tid);

        shared_active_version_history[tid].store(-1);
        int my_TS = getMyTs(); // get the timestamp and increment the global timestamp
        int expected_val = -1;
        if(!shared_active_version_history[tid].compare_exchange_strong(expected_val, my_TS)){
            my_TS = shared_active_version_history[tid];
        }
        shared_active_version_history[tid].store(my_TS);
                                                                                                                       // get the timestamp and increment the global timestamp
        priority_queue<SSSP_PQ_node<v_version>, vector<SSSP_PQ_node<v_version>>, Compare_Distance<v_version>> PQ; // priority queue for Dijkstra's algorithm
        pair<long, bool> init_pair = make_pair(LONG_MAX, false);                                                                                // PAIR: < distance INF, whether its Neighbours were iterated after this vertex was finalised>
        vector<pair<long, bool>> distance(source + 1, init_pair);                                                                               // distance vector
        v_version *u;
        bool flag = ContainsV(source, &u, logfile, debug, tid);
        if (!flag)
        {
            vector<pair<long, bool>> zero_size_vector(1, init_pair);
            return zero_size_vector;
        }
        u = find_my_version(u, my_TS, tid);
        if (u == nullptr)
        { // no appropriate version found
            vector<pair<long, bool>> zero_size_vector(1, init_pair);
            return zero_size_vector;
        }

        // push the source into PQ and update the distance vector accordingly
        distance[source].first = 0;
        distance[source].second = false;
        SSSP_PQ_node<v_version> new_PQ_node(source, 0, u);
        PQ.push(new_PQ_node);

        while (!PQ.empty())
        {
            SSSP_PQ_node<v_version> top = PQ.top();
            PQ.pop();
            if (distance[top.val].second) // if node's neighbours were already iterated after this node was finalised (since PQ doesn't allow updation)
            {
                continue;
            }
            // finalise the top.val vertex
            distance[top.val].first = top.distance;
            distance[top.val].second = true; // next time we will ignore it
            // iterate its edges
            eNode *curr_edge;
            v_version *top_dest_vertex = top.dest_vertex;
            if(top_dest_vertex)
                curr_edge =    top_dest_vertex->enext.load();
            if(curr_edge == nullptr){
                continue;
            }
            curr_edge = curr_edge->enext.load();
            while(curr_edge->val != INT_MAX)
            {
                e_version *live_version = find_my_version(curr_edge->vhead, my_TS, tid);
                if (live_version != nullptr)
                {
                    v_version *destination_version = live_version->dest_version_node;
                    if (destination_version)
                    {
                        int dest = destination_version->val;
                        if (distance.size() <= dest) // distance vector resizing
                        {
                            for (int i = distance.size(); i <= dest; i++)
                            {
                                distance.push_back(init_pair);
                            }
                        }
                        // whether path with less distance is present
                        if (distance[dest].first > top.distance + live_version->weight)
                        {
                            distance[dest].first = top.distance + live_version->weight;
                            new_PQ_node.val = dest;
                            new_PQ_node.distance = distance[dest].first;
                            new_PQ_node.dest_vertex = destination_version;
                            PQ.push(new_PQ_node);
                        }
                    }
                }
                curr_edge = curr_edge->enext.load();
            }
        }
        return distance;
    }

    
};
