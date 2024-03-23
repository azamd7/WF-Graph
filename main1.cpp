// g++ -std=c++11 -pthread -O3 main.cpp
// ./a.out logfile 128 10 6000 12000 /home/student/2022/ayaz/July23/HASH_OCT/Share_Hash/input/datasets/p2p-Gnutella06.txt 3 2 3 2 40 40 10
#include <chrono>
#include <unistd.h>
//#include <sys/resource.h>
#include <iostream>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <signal.h>
#include <sys/time.h>
#include <time.h>
#include <ctime>
#include <random>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <atomic>
#include <list>
#include <queue>
#include <stack>
#include <random>
#include <thread>
#include <map>
#include <mutex>
#include <initializer_list>
#include <unistd.h>

using namespace std;
#include "PARTIAL-LFGraph-List-List.h"
atomic<int> cnt;
time_t start1, end1;
atomic<long> vertexID;
double seconds;
// int initial_vertices = 0;
typedef struct
{
    int secs;
    int usecs;
} TIME_DIFF;



atomic<int> ops;
atomic<int> opsUpdate;
atomic<int> opsLookup;
atomic<int> opsapp;
vector<double> op;
int nops;
int itr = 0;
double add_e = 0.03;
double add_v = 0.03;
double con_e = 0.39;
double con_v = 0.39;
double rem_e = 0.02;
double rem_v = 0.02;
double BFS = 0;
double SSSP = 0;
double SS = 0;
double do_bg_GC = 0.02;

// float sp = 0.01;

int test_duration;
int initial_vertices = 10000;
// int sp = 0;
vector<double> dist_prob1 = {1, 1, 1, 1, 1, 1, 1};
atomic<bool> continue_exec;
atomic<bool> flag_gc;
atomic<int> graph_writer;
std::mutex g_num_mutex;



/**
 * @brief paramteter that are to be passed on to the threads
 *
 */
struct thread_args
{
    LFGraph *graph;
    string logfilename;
    int thread_num;
    bool debug;
    int max_nodes;
    int max_threads;
    //    double *ops;
    double *max_times;
    double *avg_times;

    int *insertV_cnt;
    int *removeV_cnt;
    int *conatinsV_cnt;
    int *insertE_cnt;
    int *removeE_cnt;
    int *conatinsE_cnt;
    int *main_operation_cnt;
    // vector<double> dist_prob ;
};

/**
 * @brief
 *
 * prob_arr will denote prob of different operations
 * 0->add vertex
 * 1->delete vertex
 * 2->add edge
 * 3->delete edge
 * 4->contains edge
 * 5->contains vertex
 * 6->snapshot
 *
 * @param t_args
 * @return ** void*
 */



void *thread_funct(void *t_args)
{

    //std::uniform_real_distribution<double> ratio_dis(0, 1.000001);
    std::uniform_real_distribution<double> ratio_dis(0, 1.001);
    //std::uniform_real_distribution<double> ratio_dis(0, 1.2);

    std::random_device rd;
    std::mt19937 gen(rd());
    LFGraph *graph = ((struct thread_args *)t_args)->graph;
    string logFileName_ip = ((struct thread_args *)t_args)->logfilename;
    int thread_num = ((struct thread_args *)t_args)->thread_num;
    bool debug = ((struct thread_args *)t_args)->debug;
    int max_nodes = ((struct thread_args *)t_args)->max_nodes + 1;
    int max_threads = ((struct thread_args *)t_args)->max_threads;
    //    double *ops = ((struct thread_args *)t_args)->ops;
    double *max_times = ((struct thread_args *)t_args)->max_times;
    double *avg_times = ((struct thread_args *)t_args)->avg_times;

    int *insertV_cnt = ((struct thread_args *)t_args)->insertV_cnt;
    insertV_cnt[thread_num] = 0;
    int *removeV_cnt = ((struct thread_args *)t_args)->removeV_cnt;
    removeV_cnt[thread_num] = 0;
    int *conatinsV_cnt = ((struct thread_args *)t_args)->conatinsV_cnt;
    conatinsV_cnt[thread_num] = 0;

    int *insertE_cnt = ((struct thread_args *)t_args)->insertE_cnt;
    insertE_cnt[thread_num] = 0;
    int *removeE_cnt = ((struct thread_args *)t_args)->removeE_cnt;
    removeE_cnt[thread_num] = 0;
    int *conatinsE_cnt = ((struct thread_args *)t_args)->conatinsE_cnt;
    conatinsE_cnt[thread_num] = 0;

    int *main_operation_cnt = ((struct thread_args *)t_args)->main_operation_cnt;
    main_operation_cnt[thread_num] = 0;

    // vector<double> dist_prob = ((struct thread_args *)t_args)->dist_prob;
    bool main_op_print = true;
    vector<double> tts; // list of time taken for snapshot
    fstream logfile_th;
    if (debug)
    {
        string logFileName = logFileName_ip + "/" + to_string(thread_num) + ".txt";
        logfile_th.open(logFileName, fstream::out);
        logfile_th << "Hi debugger from " << thread_num << endl;
    }
    //    random_device rd;
    //    mt19937 gen(rd());
    discrete_distribution<> d(dist_prob1.begin(), dist_prob1.end());
    // if(thread_num == 1)
    //     cout<<thread_num<<" "<<dist_prob[0]<<" "<<dist_prob[1]<<" "<<dist_prob[2]<<" "<<dist_prob[3]<<" "<<dist_prob[4]<<" "<<dist_prob[5]<<" "<<dist_prob[6]<<endl;
    // cout<<"Hi from thread "<<thread_num<<endl;
    int op_index;
    std::uniform_int_distribution<int64_t> index_dis(0, max_nodes);
    std::uniform_int_distribution<int64_t> index_dis_10(0, max_nodes / 10);
    std::uniform_int_distribution<int64_t> index_dis_3(0, max_nodes / 3);

    cnt++;
    while (!continue_exec);
    while (continue_exec)
    {
        double d = (double)ratio_dis(gen);
        //        cout<<d<<"\n";
        if (d <= add_v)
        {
            int rand_node_id = index_dis(gen);
            graph->AddV(rand_node_id, &logfile_th, debug, thread_num);
            insertV_cnt[thread_num]++;
        }
        else if (d <= add_v + rem_v)
        {
            int rand_node_id = index_dis_10(gen);
            graph->RemoveV(rand_node_id, &logfile_th, debug, thread_num);
            removeV_cnt[thread_num]++;
        }
        else if (d <= add_v + rem_v + add_e)
        {
            int rand_source = index_dis_3(gen);
            int rand_dest = index_dis_3(gen);
            int weight = index_dis(gen);
            while (rand_dest == rand_source)
            {
                rand_dest = index_dis(gen);
            }
            graph->AddE(rand_source, rand_dest, weight, &logfile_th, debug, thread_num);
            insertE_cnt[thread_num]++;
        }
        else if (d <= add_v + rem_v + add_e + rem_e)
        {
            int rand_source = index_dis_10(gen);
            int rand_dest = index_dis_10(gen);
            while (rand_dest == rand_source)
            {
                rand_dest = index_dis(gen);
            }
            graph->RemoveE(rand_source, rand_dest, &logfile_th, debug, thread_num);
            removeE_cnt[thread_num]++;
        }
        else if (d <= add_v + rem_v + add_e + rem_e + con_e)
        {
            int rand_source = index_dis(gen);
            int rand_dest = index_dis(gen);
            while (rand_dest == rand_source)
            {
                rand_dest = index_dis(gen);
            }
            graph->ContainsE(rand_source, rand_dest, &logfile_th, debug, thread_num);
            conatinsE_cnt[thread_num]++;
        }
        else if (d <= add_v + rem_v + add_e + rem_e + con_v + con_e)
        {
            if(d==1)
                cout<<"d double "<<d <<" by thread "<<thread_num<<endl;
            int node_id = index_dis(gen);
            v_version* ver_ptr = nullptr;
            graph->ContainsV(node_id, &ver_ptr, &logfile_th, debug, thread_num);
            conatinsV_cnt[thread_num]++;
        }
        else if (d <= add_v + rem_v + add_e + rem_e + con_v + con_e + SS)
        {
            graph->snapshot(&logfile_th, debug, thread_num);
            //cout<<"GC-hi-thr_"<<thread_num<<"\n";

        }
        else if (d <= add_v + rem_v + add_e + rem_e + con_v + con_e + SS + BFS)
        {
            int node_id = index_dis(gen);
            graph->getBFS(node_id, &logfile_th, debug, thread_num);
        }
        else if (d <= add_v + rem_v + add_e + rem_e + con_v + con_e + SS + BFS + SSSP)
        {
            int node_id = index_dis(gen);
            graph->getSSSP(node_id, &logfile_th, debug, thread_num);
        }
        else
        {
            bool fl = false, tr = true;
            //cout<<"hahano-GC1\n";
            //if(flag_gc.load() ==  false)
            if (g_num_mutex.try_lock()){
                int local_gc_cnt_old = graph_writer.load();
                int local_gc_cnt_new = local_gc_cnt_old+1;
                //cout<<"haha-GC-"<<local_gc_cnt_new<<" by thread "<<thread_num<<"\n";

                graph_writer.compare_exchange_strong(local_gc_cnt_old, local_gc_cnt_new);
                graph->do_GC(&logfile_th, debug, thread_num, local_gc_cnt_new);
                g_num_mutex.unlock();
            }
            main_operation_cnt[thread_num]++;
        }
        op[thread_num]++;
    }
    if (debug)
        logfile_th.close();
    return nullptr;
}

LFGraph *create_graph_from_file(const std::string &file, fstream *logfile, bool debug, int &a, int &b,
                                record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, e_version> * myRecManager1,
                                record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, v_version> * myRecManager2,
                                record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, eNode> * myRecManager4){
    ifstream cinn;
    cinn.open(file);
    if (cinn.fail()){
        (*logfile) << "OPEN FAILE. Bye Bye \n";
    }
    long n, m;
    int u, v, w;
    cinn >> n >> m;
    //LFGraph *graph = nullptr;
    LFGraph *graph = new LFGraph(myRecManager1, myRecManager2, myRecManager4);
    a = n;
    b = m;
    (*logfile) << file << "bye \n";
    (*logfile) << "I read n and m: " << n << m << endl;
    int i, j, e = 0;

    for (i = 1; i <= n + 1; i++)
    {
        (*logfile) << "I am adding vertex @ ------------------------------ i = " << i << endl;
        graph->AddV2(i, logfile, true, 0);
    }

    //graph->hash_resize(graph->get_hash_head(), true, 0);    // resizing just to check the printing

    for (j = 1; j <= m; j = j + 1)
    {
        cinn >> u >> v;
        int weight = (rand() % (n + 1)) + 1;
        if (u == v)
            continue;
        bool ret = graph->AddE(u + 1, v + 1, weight, logfile, true, 0);
        e++;
    }
    cout<<"AYAZ-HII"<<endl;
    /*
    for(int inx = 0; inx < 1024 ; inx++)
    {
        vNode* curr = (vNode*) graph->hash_head.load()->buckets_array[inx]->get_bucket_head();
        while(curr){
            cout<<curr->val<<"("<<curr<<") ";
            curr = curr->vnext;
        }
        cout<<endl;
    }*/
    v_version *ans = nullptr;
    cout<<graph->ContainsV(4096, &ans, logfile, true, 0);
    cout<<endl;
    return graph;
    // cout<<"Edge:"<<e<<endl;
}

int main(int argc, char **argv)
{
    // abc
    cout<<"ayaz-hi\n";
    flag_gc.store(false);
    string logFileName = "../";
    // will be used in script
    graph_writer.store(0);

    int initial_edges = 2 * (int)pow(10, 4);

    bool debug = false;

    if (argc > 1)
    {
        logFileName = argv[1];
        num_of_threads = stoi(argv[2]);
        test_duration = stoi(argv[3]);
        initial_vertices = stoi(argv[4]);
        initial_edges = stoi(argv[5]);

        add_v = ((double)atoi(argv[7]))/100;
        rem_v = ((double)atoi(argv[8]))/100;
        con_v = ((double)atoi(argv[9]))/100;
        add_e = ((double)atoi(argv[10]))/100;
        rem_e = ((double)atoi(argv[11]))/100;
        con_e = ((double)atoi(argv[12]))/100;
        SS = ((double)atoi(argv[13]))/100;
        BFS = ((double)atoi(argv[14]))/100;
        SSSP = ((double)atoi(argv[15]))/100;
//        if (argc > 8)
//        {
//            // read dist probabilities
//            for (int i = 0; i < 7; i++)
//            {
//                dist_prob1[i] = stoi(argv[7 + i]);
//            }
//        }
//        if (argc >= 15)
//        {
//            if (argv[14][0] == 'D')
//            {
//                debug = true;
//            }
//        }
    }
    string opFileName = logFileName + "/main_output.txt";
    // cout<<opFileName<<endl;
    //    op = new double[num_of_threads];
    // LFGraph * graph = create_graph(initial_vertices ,initial_edges);
    fstream opfile;
    opfile.open(opFileName, ios::out | std::ios::trunc);
    opfile << "argv[6] = " << argv[6] << endl;
    int nos_of_vertex, nos_of_edges;
    op.resize(num_of_threads+10);

    shared_active_version_history = new atomic<int>[num_of_threads];
    for(int i=0;i<num_of_threads;i++){
        shared_active_version_history[i] = INT_MAX;
    }
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, e_version> * myRecManager_edge_ver
            = new record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, e_version>(num_of_threads);
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, v_version> * myRecManager_vertex_ver
            = new record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, v_version>(num_of_threads);
    record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, eNode> * myRecManager_edge
            = new record_manager<reclaimer_debra<int>, allocator_new<int>, pool_none<int>, eNode>(num_of_threads);

    for(int i = 0; i < num_of_threads; i++){
        myRecManager_edge_ver ->initThread(i);
        myRecManager_vertex_ver ->initThread(i);
        myRecManager_edge ->initThread(i);
    }

    LFGraph *graph = create_graph_from_file(argv[6], &opfile, true, nos_of_vertex, nos_of_edges,
                                            myRecManager_edge_ver, myRecManager_vertex_ver, myRecManager_edge);
    // cout<<"Successfully created the graph "<<endl;
    // sc->print_snap_graph(&logfile);
    // printf(graph->ContainsE(5,4,1) != 2? "False\n" : "True\n");
    struct thread_args t_args[num_of_threads];
    pthread_t threads[num_of_threads];

    // List of throughput from each thread

    continue_exec.store(false);
    // cout << "End snap Enode " << end_snap_Enode << endl;
    // cout << "Marked End snap Enode " << (Snap_Enode *)get_marked_ref((long) end_snap_Enode) << endl;
    // cout << "End snap Vnode " << end_snap_Vnode << endl;
    // cout << "Marked End snap Vnode " << (Snap_Vnode*) get_marked_ref((long)end_snap_Vnode) << endl;
    double *max_times = new double[num_of_threads];
    double *avg_times = new double[num_of_threads];

    int *insertV_cnt = new int[num_of_threads];
    int *removeV_cnt = new int[num_of_threads];
    int *conatinsV_cnt = new int[num_of_threads];
    int *insertE_cnt = new int[num_of_threads];
    int *removeE_cnt = new int[num_of_threads];
    int *conatinsE_cnt = new int[num_of_threads];
    int *main_operation_cnt = new int[num_of_threads];

    for (int i = 0; i < num_of_threads; i++)
    {
        t_args[i].logfilename = logFileName;
        t_args[i].graph = graph;
        t_args[i].debug = debug;
        t_args[i].thread_num = i;
        t_args[i].max_nodes = 100 * nos_of_vertex;
        t_args[i].max_threads = num_of_threads;
        t_args[i].debug = debug;
        //        t_args[i].ops = ops;
        t_args[i].max_times = max_times;
        t_args[i].avg_times = avg_times;

        t_args[i].insertV_cnt = insertV_cnt;
        t_args[i].removeV_cnt = removeV_cnt;
        t_args[i].conatinsV_cnt = conatinsV_cnt;
        t_args[i].insertE_cnt = insertE_cnt;
        t_args[i].removeE_cnt = removeE_cnt;
        t_args[i].conatinsE_cnt = conatinsE_cnt;
        t_args[i].main_operation_cnt = main_operation_cnt;

        // t_args[i].dist_prob = dist_prob;
        pthread_create(&threads[i], NULL, thread_funct, &t_args[i]);
    }
    while (cnt < num_of_threads);
    continue_exec.store(true);
    sleep(test_duration);
    continue_exec.store(false);


    for(int i=0 ; i< num_of_threads;i++){
        pthread_join(threads[i], NULL);
    }
    // cout<<"Successfully JOINED threads "<<endl;

    double max_time = 0;
    double avg_time = 0;

    uint64_t tot = 0;
    int total_insertV_cnt = 0;
    int total_removeV_cnt = 0;
    int total_conatinsV_cnt = 0;
    int total_insertE_cnt = 0;
    int total_removeE_cnt = 0;
    int total_conatinsE_cnt = 0;
    int total_main_operation_cnt = 0;

    for (int i = 0; i < num_of_threads; i++)
    {
        // check max
        // cout<<"MAX: "<<max_times[i]<<" AVG: "<<avg_times[i]<<endl;
        if (max_time < max_times[i])
        {
            max_time = max_times[i];
        }
        avg_time += avg_times[i];
        tot += op[i];

        total_insertV_cnt += insertV_cnt[i];
        total_removeV_cnt += removeV_cnt[i];
        total_conatinsV_cnt += conatinsV_cnt[i];
        total_insertE_cnt += insertE_cnt[i];
        total_removeE_cnt += removeE_cnt[i];
        total_conatinsE_cnt += conatinsE_cnt[i];
        total_main_operation_cnt += main_operation_cnt[i];
    }

//    for(int inx = 0; inx < graph->hash_head.load()->size ; inx++)
//    {
//        vNode* curr = (vNode*) graph->hash_head.load()->buckets_array[inx]->get_bucket_head();
//        while(curr){
//            cout<<curr->val<<"("<<curr<<") ";
//            curr = curr->vnext;
//        }
//        cout<<endl;
//    }

    avg_time = avg_time / num_of_threads;

    cout << tot / test_duration << fixed <<endl;
    cout << 0 << fixed << endl;             // just to avoid editing the script file
    cout << total_insertV_cnt << fixed << endl;
    cout << total_removeV_cnt << fixed << endl;
    cout << total_conatinsV_cnt << fixed << endl;
    cout << total_insertE_cnt << fixed << endl;
    cout << total_removeE_cnt << fixed << endl;
    cout << total_conatinsE_cnt << fixed << endl;
    cout << total_main_operation_cnt << fixed << endl;
    cout<< num_of_threads<<"    "<<(double)add_v <<"   "<<rem_v<<"   "<<con_v<<"   "<<add_e <<"   "<<rem_e<<"   "<<con_e<<"   "<<SS <<"   "<<BFS<<"   "<<SSSP<<"   "<<tot / test_duration<<"\n";
    cout<<"FINAL BYE \n";
    cout<<"Hash table size "<<graph->hash_head.load()->size<<endl;
    return 0;
}
// 64