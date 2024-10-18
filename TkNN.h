#include <vector>
using namespace std;

#define INT_MAX 999999999

struct HASH_NODE
{
	int a, b, c;
	int next;
	HASH_NODE() {}
	HASH_NODE(int _a, int _b, int _c)
	{
		a = _a;
		b = _b;
		c = _c;
	}
};
class HASH
{
public:
	static const int P = 100000007;

	vector<HASH_NODE> nodes;
	vector<int> start;
	HASH()
	{
		nodes.clear();
		nodes.push_back(HASH_NODE());
		start.resize(P);
		for (int i = 0; i < P; i++)
			start[i] = 0;
	}
	inline int get_id(int a, int b)
	{
		return ((long long)(a)*b) % P;
	}
	inline bool is_exist(int a, int b)
	{
		int t = get_id(a, b);
		int p = start[t];
		while (p != 0)
		{
			if (nodes[p].a == a && nodes[p].b == b)
				return true;
			p = nodes[p].next;
		}
		return false;
	}
	inline int get_value(int a, int b)
	{

		int t = get_id(a, b);
		//   if (a == 1476 && b == 30137)
		//      cout << "t: " << t << "  start[t]: " << start[t] << endl;
		int p = start[t];
		//   if (nodes[p].a == a && nodes[p].b == b)
		//      return nodes[p].c;
		//   else return -1;
		while (p != 0)
		{
			//    if (a == 1476 && b == 30137)
			//         cout << nodes[p].a << "  --  " << nodes[p].b << "  --  " << nodes[p].c << endl;
			if (nodes[p].a == a && nodes[p].b == b)
				return nodes[p].c;
			p = nodes[p].next;
		}
		return -1;
	}
	inline bool insert_node(int a, int b, int c)
	{

		if (is_exist(a, b) == true)
			return false;
		int t = get_id(a, b);
		//    if (a == 1476 && b == 30137){
		//        cout << "hehehe: " << c << endl;
		//       cout << "id: " << t << endl;
		//   }
		HASH_NODE hn = HASH_NODE(a, b, c);
		hn.next = start[t];
		nodes.push_back(hn);
		start[t] = nodes.size() - 1;
		return true;
	}
	inline bool delete_node(int a, int b)
	{
		int t = get_id(a, b);
		int p = start[t];
		int pre = -1;
		while (p != 0)
		{
			if (nodes[p].a == a && nodes[p].b == b)
			{
				if (pre < 0)
					start[t] = nodes[p].next;
				else
					nodes[pre].next = nodes[p].next;
				return true;
			}
			pre = p;
			p = nodes[p].next;
		}
		return false;
	}
};

class TkNN
{
public:
	struct list_node
	{
		int previous, next, key, dist;
		list_node()
		{
			previous = -1;
			next = -1;
			key = -1;  // vertex
			dist = -1; // distance
		}
	};
	struct object_saveing_structure
	{

		vector<list_node> a;
		vector<int> trash_can; // used to recover the ktnn in query_delay
		int current, size_num;
		object_saveing_structure()
		{
			a.clear();
			list_node _a;
			_a.key = -1;
			_a.dist = -1;
			_a.previous = 0;
			_a.next = 0;
			a.push_back(_a);
			trash_can.clear();
			current = 0;
			size_num = 0;
		}
	};
    void initialize_knn(char* index_file, char* obj_file);
    HASH hash;
	vector<double> times_period[20];
	int period = 20;
	vector<object_saveing_structure> OSS;
	vector<int> object_number;
	vector<int> user_number;
	vector<vector<pair<int, int>>> path_from_root;
    void create_kNN_index();
    void delete_element(int p, int x);
    inline void insert(int p, int x);
    inline void insert(int x);
    void get_all_object(int p, vector<int> &a);
    void OSS_push_back(int p, int key, int dist);
    void OSS_push_front(int p, int key, int dist);
    void get_subtree(int p, int key, vector<int> &a);
    void dfs_sort(int p, vector<int> &a);
    void join_subtree(int p, int x, vector<int> &a);
    void dfs_neighbor(int p);
    void traversal(int p);
    void compute_object_number();
    void initialize_object();
    int *query_mark;
	int query_mark_stamp;
	vector<double> time_save;
	vector<pair<int, int>> query(int x, int top_k, int limit_dis = INT_MAX);
	bool Verify(int x, int target, int top_k);
    vector<pair<int, int>> query_naive(int x, int top_k);
    vector<pair<int, int>> query_delay(int x, int top_k);
    void object_setting(int n);
    bool double_objects(int p);
    bool check_everyone();
    void print(int root);

};

int distanceQueryFull(int p, int q);
int distanceQuery(int p, int q);
// int SPNQuery(int x, int target);
int LCAQuery(int _p, int _q);
// int distanceQueryFull(int p, int q);