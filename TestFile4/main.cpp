#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <set>
#include <bitset>
#include <iomanip>
using namespace std;


//Please define MEOW in your compiler.
//Then you can compile on your computer first to see if the syntax is alright.
#define this it
#ifdef MEOW
#include "bc.h"
#else
#include <bc.h>
#endif
#undef this

/// See bc.h at xxx for the API you have access too.
/// Note: the API is not thread safe; don't use pthreads.
/// It's not very pretty, sorry. If you want a prettier low-level language, maybe consider rust?

// Any method in the API may set an error.
// Call check_errors() to get the most recent error.
bool check_errors(string s) {
    /// Check if we have an error...
    if (bc_has_err()) {
        char *err;
        /// Note: this clears the current global error.
        int8_t code = bc_get_last_err(&err);
        cout<<"Error occurs while "<<s<<endl;
        printf("Engine error code %d: %s\n", code, err);
        bc_free_string(err);
        return true;
    } else {
        return false;
    }
}

//Unify the names of the delete functions.
inline void release(bc_PlanetMap *planetmap){delete_bc_PlanetMap(planetmap);}
inline void release(bc_GameController *gc){delete_bc_GameController(gc);}
inline void release(bc_Unit *unit){delete_bc_Unit(unit);}
inline void release(bc_Location *location){delete_bc_Location(location);}
inline void release(bc_MapLocation *maplocation){delete_bc_MapLocation(maplocation);}
inline void release(bc_VecUnit *units){delete_bc_VecUnit(units);}
inline void release(bc_VecMapLocation *vecmaplocation){delete_bc_VecMapLocation(vecmaplocation);}
inline void release(bc_ResearchInfo *research_info){delete_bc_ResearchInfo(research_info);}
inline void release(bc_VecUnitID *units){delete_bc_VecUnitID(units);}
inline void release(bc_OrbitPattern *orbitpattern){delete_bc_OrbitPattern(orbitpattern);}
inline void release(bc_AsteroidPattern *asteroidpattern){delete_bc_AsteroidPattern(asteroidpattern);}
inline void release(bc_AsteroidStrike *asteroidstrike){delete_bc_AsteroidStrike(asteroidstrike);}

//Smart pointer. PLEASE ALWAYS USE SMART POINTER
//This prevents memory leaks and saves time :)
template <typename T>
struct Ptr
{
    struct Counter
    {
        T *ptr;
        uint16_t cnt;
        Counter(T *tmp):ptr(tmp), cnt(1){};
    };
    Counter *cntr;
    inline void inc(){cntr->cnt++;}
    inline void dec()
    {
        cntr->cnt--;
        if(!cntr->cnt && cntr->ptr) release(cntr->ptr), delete cntr;
    }
    Ptr(const Ptr& s):cntr(s.cntr){inc();}
    Ptr(T *s = 0){cntr = new Counter(s);}
    const Ptr& operator = (const Ptr& s)
    {
        if(cntr) dec();
        cntr = s.cntr;
        if(cntr) inc();
        return *this;
    }
    ~Ptr(){dec();}
    operator T*(){return cntr->ptr;}
    operator bool(){return cntr->ptr;}
    T& operator *(){return *cntr->ptr;}
    T* operator ->(){return cntr->ptr;}
};

//A smart pointer to bc_GameController type. It's the game controller that we'll use.
Ptr<bc_GameController> gc;

bc_Planet my_Planet; //Current planet
bc_Team my_Team; //Current team
int map_height[2], map_width[2]; //The size of the Earth map and the Mars map.
vector<bool> passable[2]; //The current map
int karbonite[2500];
vector<int> connected_square_num; //The size of connected components on Mars.
int max_connected_num;
int passable_count[2];
unsigned int shortest_distance[2500][2500]; //Shortest distance between squares. (x,y) corresponds to y*w+h.
unsigned int distance_to_wall[2500]; //Distance to impassable squares/the bound of the map.
Ptr<bc_PlanetMap> planetmap[2]; //The maps.
default_random_engine gen; //C++11 random engine
bool can_build_rocket;
map<int, int> building_rocket; //(id, loc)
map<int, pair<int,int>> built_rocket; //(id, (worker_num, total_num))
map<int, int> built_rocket_location;
map<int, int> not_free; //For Bug pathing. not_free[id] = (destination_loc, last_direction);
set<int> alive_rockets;
Ptr<bc_OrbitPattern> orbit_pattern;
int symmetry; //The guessed symmetry of this game. 0: identity(impossible). 1: reflection about y. 2: reflection about x. 3: rotation
//Cut the map into chunks
int chunk_num;
int chunk_label[2500]; //The number of chunk a square belongs to
vector<vector<pair<int,int>>> chunk_edge; //The connection between chunks
vector<int> chunk_size;
vector<int> chunk_karbonite;
vector<int> chunk_rep; //A representative of a chunk. Just for estimation.
vector<int> chunk_friend_fire;
vector<int> chunk_enemy_fire;
vector<int> chunk_next_chunk_to_karbonite;
int structure_wait_for_work_round[65536];
int worker_build_target[65536];
bool worker_is_building[65536];
int target_rocket[65536];
map<int, int> need_bot_rocket;
int poked_direction[65536];
int last_attack_round[65536];
int last_snipe_round[65536];
int idle_round[65536];
int enemy_direction[65536];
bool visible[2500];
Ptr<bc_Unit> units[2500];
int enemies_max_total_damage[2500];
int connected_comp_label[2500];
bool being_harvested[2500];
vector<pair<pair<int,int>,Ptr<bc_Unit>>> enemies;
vector<pair<pair<int,int>,Ptr<bc_Unit>>> teammates;
map<int, int> my_factories;
set<int> alive_factories;
vector<int> typecount; //Count the number of robots of a certain type.
typedef vector<pair<pair<int,int>,Ptr<bc_Unit>>>::iterator vit;
bool try_blueprint_rocket;
bool success_bluesprint_rocket;
bool need_worker;
bool first_enemy;
bool should_stay[65536];
int invisible_loc;
bool should_build_rocket;
bool rocket_just_blueprinted;
Ptr<bc_AsteroidPattern> asteroid_pattern;
int waiting_for_rocket;
int have_enemy_round;
vector<int> left_worker_num;
set<int> on_rocket;
set<int> overcharge_healer;


vit find_by_x(vit first, vit last, int x)
{
    int l = 0, r = last-first+1;
    while(r-l > 1)
    {
        int m = (l+r)/2;
        if((first+m-1)->first.first >= x) r = m;
        else l = m;
    }
    return first+l;
}

vit find_by_y(vit first, vit last, int y)
{
    int l = 0, r = last-first+1;
    while(r-l > 1)
    {
        int m = (l+r)/2;
        if((first+m-1)->first.second >= y) r = m;
        else l = m;
    }
    return first+l;
}

int corresponding_point(int loc, int sym) // Returns the location of the given location projected according to the given symmetry
{
    /*
    1 : Y-symmetry
    2 : X-symmetry
    3 : 180deg
    */
    int x = loc%map_width[my_Planet], y = loc/map_width[my_Planet];
    if(sym&1) x = map_width[my_Planet]-1-x;
    if(sym&2) y = map_height[my_Planet]-1-y;
    return map_width[my_Planet]*y+x;
}

bool out_of_bound(int loc, int dir, bc_Planet planet) //Test if going in the direction dir from (loc%w, loc/w) will go out the map
{
    int h = map_height[planet], w = map_width[planet];
    if(dir == 0) return loc/w == h-1;
    if(dir == 1) return loc/w == h-1 || loc%w == w-1;
    if(dir == 2) return loc%w == w-1;
    if(dir == 3) return loc%w == w-1 || loc/w == 0;
    if(dir == 4) return loc/w == 0;
    if(dir == 5) return loc/w == 0 || loc%w == 0;
    if(dir == 6) return loc%w == 0;
    if(dir == 7) return loc%w == 0 || loc/w == h-1;
    if(dir == 8) return false;
}

bool out_of_bound(int loc, int dir)
{
    return out_of_bound(loc, dir, my_Planet);
}

int go(int dir, bc_Planet planet) // Add what this function returns to a location variable to move in the given direction
{
    int w = map_width[planet];
    int a[9] = {w, w+1, 1, -w+1, -w, -w-1, -1, w-1, 0};
    return a[dir];
}

int go(int dir) // In case you don't specify the planet
{
    return go(dir, my_Planet);
}

int dfs(int f, int s, vector<bool>& p, vector<bool>& visited, vector<int>& fa) // Depth-first search (recursive) for Mars
{
    /*
    Performs a depth-first search, returning the total area of pockets on Mars
    s	    current node (position value)
    visited	array of positions recording whether or not dfs function has logged their areas
    fa	    father
    f       also father? Honestly I don't know what the two do
    p	    parent
    */
    visited[s] = 1;
    fa[s] = f;
    int ans = 1;
    for(int i = 0; i < 8; i++)
    {
        // Conditions to end recursion
        if(out_of_bound(s, i, Mars)) continue; // Out of bounds
        if(!p[s+go(i, Mars)] || visited[s+go(i, Mars)]) continue; // Impassable or visited
        ans += dfs(f, s+go(i, Mars), p, visited, fa); // Else, recur
    }
    return ans;
}

//Disjoint set
int _ds__p[2500];
int _ds__sz[2500];

inline void _ds__init(int n){for(int i = 0; i < n; i++) _ds__p[i] = i, _ds__sz[i] = 1;}

int _ds__find(int s){return (_ds__p[s] == s)?s:_ds__p[s]=_ds__find(_ds__p[s]);}

void _ds__union(int a, int b)
{
    a = _ds__find(a), b = _ds__find(b);
    if(a == b) return;
    if(_ds__sz[a] < _ds__sz[b])swap(a,b);
    _ds__sz[a] += _ds__sz[b];
    _ds__p[b] = a;
}

inline int abs(int a){return (a>0)?a:-a;}

void organize_map_info() //Organize all the map information
{
    for(int i = 0; i < 2; i++)
    {
        map_height[i] = bc_PlanetMap_height_get(planetmap[i]);
        map_width[i] = bc_PlanetMap_width_get(planetmap[i]);
        passable[i].resize(map_height[i]*map_width[i]);
        for(int x = 0; x < map_width[i]; x++) for(int y = 0; y < map_height[i]; y++) //Read in the map
        {
            Ptr<bc_MapLocation> tmp(new_bc_MapLocation(bc_Planet(i), x, y));
            passable_count[i] += passable[i][map_width[i]*y+x] = bc_PlanetMap_is_passable_terrain_at(planetmap[i], tmp);

            if(bc_Planet(i) == my_Planet) karbonite[map_width[i]*y+x] = bc_PlanetMap_initial_karbonite_at(planetmap[i], tmp);
        }
    }
    int h = map_height[my_Planet], w = map_width[my_Planet];
    vector<bool>& p = passable[my_Planet];
//    Guess the symmetry
    if(my_Planet == Earth)
    {
        bool possibly_1 = true, possibly_2 = true;

        // Check for symmetry = 2
        for(int loc = 0; loc < h*w/2; loc++){
            if(p[loc] != p[corresponding_point(loc, 2)])
            {
                possibly_2 = false;
                break;
            }
        }

        // Check for symmetry = 1
        for(int loc = 0; loc < h*w/2; loc++){ //Note this checks half the map, plus the other half through corresponding_point.
            int temp_x = loc/h, temp_y = loc%h; // Convert from row-major to column-major
            int inverted_loc = h*temp_y+temp_x;
            if(p[inverted_loc] != p[corresponding_point(inverted_loc, 1)])
            {
                possibly_1 = false;
                break;
            }
        }
        if (possibly_1 == possibly_2) symmetry = 3;
        else symmetry = possibly_1 ? 1 : 2;
    }
    else asteroid_pattern = bc_GameController_asteroid_pattern(gc);

//    BFS shortest path
    for(int i = 0; i < h*w; i++) for(int j = 0; j < h*w; j++) shortest_distance[i][j] = (i==j)?0:-1;// -1 = Very large
    for(int i = 0; i < h*w; i++) // BFS.
    {
        distance_to_wall[i] = -1; //-1 = Very large
        if(!p[i]) continue;
        queue<int> q; q.push(i);
        while(q.size())
        {
            int j = q.front(); q.pop();
            for(int k = 0; k < 8; k++)
            {
                if(!out_of_bound(j, k) && p[j+go(k)])
                {
                    if(shortest_distance[i][j+go(k)] == -1)
                    {
                        shortest_distance[i][j+go(k)] = shortest_distance[i][j]+1;
                        q.push(j+go(k));
                    }
                }
                else distance_to_wall[i] = min(distance_to_wall[i], shortest_distance[i][j]+1);
            }
        }
    }
//    Simplify the map
    vector<bitset<2500>> can_see(h*w);

//    for(int i = 0; i < h*w; i++)
//    {
//        if(!p[i]) continue;
//        int ix = i%w, iy = i/w;
//        for(int j = 0; j < h*w; j++)
//        {
//            if(!p[j]) continue;
//            int jx = j%w, jy = j/w;
//            if(shortest_distance[i][j] == max(abs(ix-jx), abs(iy-jy))) can_see[i][j] = 1;
//        }
//    }

    for(int i = 0; i < h*w; i++)
    {
        if(!p[i]) continue;
        for(int j = 0; j < 4; j++)
        {
            int ii = i;
            while(1)
            {
                int iii = ii;
                while(1)
                {
                    can_see[i][iii] = 1;
                    if(out_of_bound(iii, 2*j) || !p[iii+go(2*j)]) break;
                    iii += go(2*j);
                }
                if(out_of_bound(ii, (2*j+2)%8) || !p[ii+go((2*j+2)%8)]) break;
                ii += go((2*j+2)%8);
            }
            ii = i;
            while(1)
            {
                int iii = ii;
                while(1)
                {
                    can_see[i][iii] = 1;
                    if(out_of_bound(iii, 2*j) || !p[iii+go(2*j)]) break;
                    iii += go(2*j);
                }
                if(out_of_bound(ii, (2*j+6)%8) || !p[ii+go((2*j+6)%8)]) break;
                ii += go((2*j+6)%8);
            }
        }
    }

//    BFS again
    queue<int> qq;
    vector<bool> visited(h*w);
    for(int i = 0; i < h*w; i++) if(p[i])
    {
        qq.push(i); visited[i] = 1;
        break;
    }
    _ds__init(h*w);
    while(qq.size())
    {
        int j = qq.front(); qq.pop();
        for(int k = 0; k < 8; k++)
            if(!out_of_bound(j, k) && p[j+go(k)])
            {
                if(!visited[j+go(k)]) visited[j+go(k)] = 1, qq.push(j+go(k));
                else if(!(can_see[j]^can_see[j+go(k)]).count()) _ds__union(j, j+go(k));
            }
    }

    fill(visited.begin(), visited.end(), 0);
    for(int i = 0; i < h*w; i++) if(p[i])
    {
        qq.push(i); visited[i] = 1;
        break;
    }
    while(qq.size())
    {
        int j = qq.front(); qq.pop();
        for(int k = 0; k < 8; k++)
            if(!out_of_bound(j, k) && p[j+go(k)])
            {
                if(!visited[j+go(k)]) visited[j+go(k)] = 1, qq.push(j+go(k));
                else if((can_see[j]^can_see[j+go(k)]).count() <= h*w/10/min(_ds__sz[_ds__find(j)], _ds__sz[_ds__find(j+go(k))])) _ds__union(j, j+go(k));
            }
    }
    fill(chunk_label, chunk_label+h*w, -1);
    for(int i = 0; i < h*w; i++) if(p[i] && _ds__p[i] == i) chunk_label[i] = chunk_num++;
    for(int i = 0; i < h*w; i++) if(p[i]) chunk_label[i] = chunk_label[_ds__find(i)];
    chunk_edge.resize(chunk_num);
    chunk_size.resize(chunk_num);
    chunk_karbonite.resize(chunk_num);
    chunk_rep.resize(chunk_num);
    chunk_enemy_fire.resize(chunk_num);
    chunk_friend_fire.resize(chunk_num);
    chunk_next_chunk_to_karbonite.resize(chunk_num, -1);
    vector<pair<unsigned int,double>> chunk_dist_to_wall(chunk_num);
    vector<pair<double,double>> chunk_centroid(chunk_num);
    for(int i = 0; i < h*w; i++)
    {
        if(!p[i]) continue;
        chunk_size[chunk_label[i]]++;
        chunk_karbonite[chunk_label[i]] += karbonite[i];
        chunk_centroid[chunk_label[i]].first += i%w;
        chunk_centroid[chunk_label[i]].second += i/w;
        set<int> added;
        for(int j = 0; j < 8; j++) if(!out_of_bound(i,j) && p[i+go(j)] && chunk_label[i] != chunk_label[i+go(j)] && !added.count(chunk_label[i+go(j)]))
        {
            added.insert(chunk_label[i+go(j)]);
            bool found_it = 0;
            for(int k = 0; k < chunk_edge[chunk_label[i]].size(); k++)
                if(chunk_edge[chunk_label[i]][k].first == chunk_label[i+go(j)])
                {
                    chunk_edge[chunk_label[i]][k].second++;
                    found_it = 1;
                }
            if(!found_it) chunk_edge[chunk_label[i]].push_back(make_pair(chunk_label[i+go(j)], 1));
        }
    }
    for(int i = 0; i < chunk_num; i++) chunk_centroid[i].first /= chunk_size[i], chunk_centroid[i].second /= chunk_size[i];
    for(int i = 0; i < h*w; i++)
    {
        if(!p[i]) continue;
        int label = chunk_label[i];
        auto tmp = make_pair(distance_to_wall[i], -abs(i%w-chunk_centroid[label].first) - abs(i/w-chunk_centroid[label].second));
        if(tmp > chunk_dist_to_wall[label]) chunk_dist_to_wall[label] = tmp, chunk_rep[label] = i;
    }
//    DFS connected component on Mars for rockets landing
    if(my_Planet == Earth)
    {
        h = map_height[Mars], w = map_width[Mars];
        vector<bool>& q = passable[Mars];
        vector<int> fa(h*w);
        fill(visited.begin(), visited.end(), 0);
        visited.resize(h*w);
        connected_square_num.resize(h*w);
        for(int i = 0; i < h*w; i++) if(q[i])
        {
            if(!visited[i]) connected_square_num[i] = dfs(i, i, q, visited, fa);
            else connected_square_num[i] = connected_square_num[fa[i]];
            max_connected_num = max(max_connected_num, connected_square_num[i]);
        }
    }
    check_errors("Organizing");
    return;
}

bool is_robot(bc_UnitType s) //Check if a unit is a robot.
{
    return (s != Rocket) && (s != Factory);
}

inline int lowbit(int x){return x&(-x);}

vector<int> get_random_indices(vector<int>& weight) //Get random indices according probabilities proportional to weights
{
//    I use Fenwick tree to expedite from O(n^2) to O(nlgn) here. Hope this is faster.
    int sum = 0;
    vector<int> bit(weight.size()+1);
    for(int i = 0; i < weight.size(); i++) bit[i+1] = (sum += weight[i]);
    for(int i = weight.size(); i; i--)bit[i] -= bit[i-lowbit(i)];
    uniform_int_distribution<int> dis(0, 2147483647);
    vector<int> ret;
    int lg = 0, tmp = 1;
    while(tmp <= weight.size()) lg++, tmp *= 2;
    while(sum)
    {
        int random_number = dis(gen)%sum, ind = 0, tmpsum = 0;
        for(int i = lg; i >= 0; i--)
            if(ind + (1<<i) <= weight.size() && tmpsum + bit[ind+(1<<i)] <= random_number)
            {
                ind += (1<<i);
                tmpsum += bit[ind];
            }
        ret.push_back(ind);
        sum -= weight[ind];
        int s = ind+1;
        while(s <= weight.size())
        {
            bit[s] -= weight[ind];
            s += lowbit(s);
        }
    }
    return ret;
}

int get_total_damage(int loc)
{
    if(enemies_max_total_damage[loc] != -1) return enemies_max_total_damage[loc];
    enemies_max_total_damage[loc] = 0;
    int dist = 50;
    int sq = 0; while((sq+1)*(sq+1) <= dist) sq++;
    int now_x = loc%map_width[my_Planet], now_y = loc/map_width[my_Planet];
    auto first = find_by_x(enemies.begin(), enemies.end(), now_x-sq);
    auto last = find_by_x(enemies.begin(), enemies.end(), now_x+sq+1);
    if(last-first < 22*dist/7) for(auto it = first; it != last; it++)
    {
        int x = it->first.first, y = it->first.second;
        bc_UnitType type = bc_Unit_unit_type(it->second);
        if(!is_robot(type) || type == Worker || type == Healer) continue;
        if((x-now_x)*(x-now_x)+(y-now_y)*(y-now_y) <= bc_Unit_attack_range(it->second))
            enemies_max_total_damage[loc] += bc_Unit_damage(it->second);
    }
    else for(int i = 0; i <= sq; i++) for(int neg_x = -1; neg_x < (i?3:1); neg_x+=2)
        {
            if(now_x+neg_x*i < 0 || now_x+neg_x*i >= map_width[my_Planet]) continue;
            for(int j = 0; j*j+i*i <= dist; j++) for(int neg_y = -1; neg_y < (j?3:1); neg_y += 2)
            {
                if(now_y+neg_y*j < 0 || now_y+neg_y*j >= map_height[my_Planet]) continue;
                int new_loc = now_x+neg_x*i + (now_y+neg_y*j)*map_width[my_Planet];
                if(units[new_loc] && bc_Unit_team(units[new_loc]) != my_Team)
                {
                    bc_UnitType type = bc_Unit_unit_type(units[new_loc]);
                    if(!is_robot(type) || type == Worker || type == Healer) continue;
                    if(i*i+j*j <= bc_Unit_attack_range(units[new_loc]))
                        enemies_max_total_damage[loc] += bc_Unit_damage(units[new_loc]);
                    else if(type == Knight && i*i+j*j <= 10)
                        enemies_max_total_damage[loc] += bc_Unit_damage(units[new_loc]);
                }
            }
        }
    return enemies_max_total_damage[loc];
}

bool can_harvest(int id)
{
    for(int i = 0; i < 9; i++) if(bc_GameController_can_harvest(gc, id, bc_Direction(i))) return 1;
    return 0;
}

bool try_harvest(int id)
{
    for(int i = 0; i < 9; i++)
        if(bc_GameController_can_harvest(gc, id, bc_Direction(i)))
        {
            bc_GameController_harvest(gc, id, bc_Direction(i));
            return 1;
        }
    return 0;
}

bool try_build(int id, int loc)
{
    for(int i = 0; i < 8; i++)
        if(!out_of_bound(loc, i) && units[loc+go(i)] && !is_robot(bc_Unit_unit_type(units[loc+go(i)]))
            && bc_Unit_team(units[loc+go(i)]) == my_Team && bc_GameController_can_build(gc, id, bc_Unit_id(units[loc+go(i)])))
        {
            bc_GameController_build(gc, id, bc_Unit_id(units[loc+go(i)]));
            worker_build_target[id] = loc+go(i);
            worker_is_building[id] = 1;
            return 1;
        }
    return 0;
}

pair<bool,int> try_replicate(int id, int loc)
{
    for(int i = 0; i < 8; i++)
        if(bc_GameController_can_replicate(gc, id, bc_Direction(i)))
        {
            bc_GameController_replicate(gc, id, bc_Direction(i));
            int new_loc = loc+go(i);
            Ptr<bc_MapLocation> tmpmloc(new_bc_MapLocation(my_Planet, new_loc%map_width[my_Planet], new_loc/map_width[my_Planet]));
            Ptr<bc_Unit> tmp(bc_GameController_sense_unit_at_location(gc, tmpmloc));
            units[new_loc] = tmp;
            if(my_Planet == Earth) left_worker_num[connected_comp_label[loc]]--;
            return make_pair(1, new_loc);
        }
    return make_pair(0,-1);
}

bool try_blueprint(int id, int loc, bc_UnitType structure)
{
    int best_dir = -1, max_dist_to_wall = 0;
    if(structure == Rocket) try_blueprint_rocket = 1;
    for(int i = 0; i < 8; i++)
        if(bc_GameController_can_blueprint(gc, id, structure, bc_Direction(i)))
            if(distance_to_wall[loc+go(i)] > max_dist_to_wall) max_dist_to_wall = distance_to_wall[loc+go(i)], best_dir = i;
    if(best_dir == -1)
    {
        if(bc_GameController_karbonite(gc) >= bc_UnitType_blueprint_cost(structure))
        {
            for(int i = 0; i < 8; i++) if(!out_of_bound(loc, i) && passable[my_Planet][loc+go(i)])
            {
                Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, (loc+go(i))%map_width[my_Planet], (loc+go(i))/map_width[my_Planet]));
                if(bc_GameController_has_unit_at_location(gc, tmp))
                {
                    Ptr<bc_Unit> blocking_unit(bc_GameController_sense_unit_at_location(gc, tmp));
                    if(bc_Unit_team(blocking_unit) == my_Team && is_robot(bc_Unit_unit_type(blocking_unit)))
                        poked_direction[bc_Unit_id(blocking_unit)] = i;
                }
            }
        }
        return 0;
    }
    bc_GameController_blueprint(gc, id, structure, bc_Direction(best_dir));
    int new_loc = loc+go(best_dir);
    Ptr<bc_MapLocation> tmpmloc(new_bc_MapLocation(my_Planet, new_loc%map_width[my_Planet], new_loc/map_width[my_Planet]));
    Ptr<bc_Unit> tmp(bc_GameController_sense_unit_at_location(gc, tmpmloc));
    units[new_loc] = tmp;
    if(structure == Factory)
    {
        passable[my_Planet][new_loc] = 0;
        for(int i = 0; i <map_width[my_Planet]*map_height[my_Planet]; i++)
            shortest_distance[i][new_loc] = shortest_distance[new_loc][i] = -1;
        alive_factories.insert(bc_Unit_id(tmp));
        my_factories[bc_Unit_id(tmp)] = new_loc;
        for(int i = 0; i < chunk_num; i++) if(chunk_rep[i] == new_loc) for(int dir = 0; dir < 8; dir++)
            if(!out_of_bound(new_loc, dir) && passable[my_Planet][new_loc+go(dir)])
            {
                chunk_rep[i] += go(dir);
                break;
            }
    }
    if(structure == Rocket) success_bluesprint_rocket = 1, rocket_just_blueprinted = 1;
    worker_is_building[id] = 1, worker_build_target[id] = new_loc;
    return 1;
}

bool try_repair(int id, int loc)
{
    for(int i = 0; i < 8; i++)
        if(!out_of_bound(loc, i) && units[loc+go(i)] && !is_robot(bc_Unit_unit_type(units[loc+go(i)]))
           && bc_Unit_team(units[loc+go(i)]) == my_Team && bc_GameController_can_repair(gc, id, bc_Unit_id(units[loc+go(i)]))
           && bc_Unit_max_health(units[loc+go(i)]) != bc_Unit_health(units[loc+go(i)]))
       {
           bc_GameController_repair(gc, id, bc_Unit_id(units[loc+go(i)]));
           return 1;
       }
    return 0;
}

pair<bool, int> try_unload(int id, int loc)
{
    int tmp[8]; for(int i = 0; i < 8; i++) tmp[i] = i;
    random_shuffle(tmp, tmp+8);
    vector<int> high_priority;
    if(enemy_direction[id] != 8)
    {
        high_priority.push_back(enemy_direction[id]);
        int clockwise = 2*(rand()%2)-1;
        high_priority.push_back((enemy_direction[id]+clockwise+8)%8);
        high_priority.push_back((enemy_direction[id]-clockwise+8)%8);
    }
    for(int i = 0; i < high_priority.size(); i++) for(int j = 0; j < 8; j++) if(tmp[j] == high_priority[i]) swap(tmp[j], tmp[i]);
    for(int i = 0; i < 8; i++) if(bc_GameController_can_unload(gc, id, bc_Direction(tmp[i])))
    {
        bc_GameController_unload(gc, id, bc_Direction(tmp[i]));
        int new_loc = loc+go(tmp[i]);
        Ptr<bc_MapLocation> tmpmloc(new_bc_MapLocation(my_Planet, new_loc%map_width[my_Planet], new_loc/map_width[my_Planet]));
        Ptr<bc_Unit> tmp(bc_GameController_sense_unit_at_location(gc, tmpmloc));
        units[new_loc] = tmp;
        return make_pair(1, new_loc);
    }
    for(int i = 0; i < 8; i++) if(!out_of_bound(loc, i) && passable[my_Planet][loc+go(i)])
    {
        Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, (loc+go(i))%map_width[my_Planet], (loc+go(i))/map_width[my_Planet]));
        if(bc_GameController_has_unit_at_location(gc, tmp))
        {
            Ptr<bc_Unit> blocking_unit(bc_GameController_sense_unit_at_location(gc, tmp));
            if(bc_Unit_team(blocking_unit) == my_Team && is_robot(bc_Unit_unit_type(blocking_unit))
                && bc_Unit_unit_type(blocking_unit) != Worker)
                poked_direction[bc_Unit_id(blocking_unit)] = i;
        }
    }
    return make_pair(0, -1);
}

bool try_produce(int id, vector<int>& weight) //weight: in proportion to the probability of producing a certain type
{
    vector<int> tmp(get_random_indices(weight));
    for(auto i:tmp)
    {
        if(bc_GameController_can_produce_robot(gc, id, bc_UnitType(i)))
        {
            bc_GameController_produce_robot(gc, id, bc_UnitType(i));
            return 1;
        }
    }
    //    Impossible to produce any unit. :(
    return 0;
}

void get_nearby_enemies(int loc, int dist, vector<Ptr<bc_Unit>>& nearby_enemies)
{
    nearby_enemies.clear();
    int sq = 0; while((sq+1)*(sq+1) <= dist) sq++;
    int now_x = loc%map_width[my_Planet], now_y = loc/map_width[my_Planet];
    auto first = find_by_x(enemies.begin(), enemies.end(), now_x-sq);
    auto last = find_by_x(enemies.begin(), enemies.end(), now_x+sq+1);
    if(last-first < 22*dist/7) for(auto it = first; it != last; it++)
    {
        int x = it->first.first, y = it->first.second;
        if((x-now_x)*(x-now_x)+(y-now_y)*(y-now_y) <= dist) nearby_enemies.push_back(it->second);
    }
    else for(int i = 0; i <= sq; i++) for(int neg_x = -1; neg_x < (i?3:1); neg_x+=2)
        {
            if(now_x+neg_x*i < 0 || now_x+neg_x*i >= map_width[my_Planet]) continue;
            for(int j = 0; j*j+i*i <= dist; j++) for(int neg_y = -1; neg_y < (j?3:1); neg_y += 2)
            {
                if(now_y+neg_y*j < 0 || now_y+neg_y*j >= map_height[my_Planet]) continue;
                int new_loc = now_x+neg_x*i + (now_y+neg_y*j)*map_width[my_Planet];
                if(units[new_loc] && bc_Unit_team(units[new_loc]) != my_Team) nearby_enemies.push_back(units[new_loc]);
            }
        }
}

bool try_attack(int id, int loc, int dist, vector<Ptr<bc_Unit>>& nearby_enemies)
{
//    vector<int> weight(nearby_enemies.size());
//    for(int i = 0; i < nearby_enemies.size(); i++)
//    {
//        bc_UnitType type = bc_Unit_unit_type(nearby_enemies[i]);
//        if(type == Worker) weight[i] = 1;
//        else if(type == Knight) weight[i] = 100;
//        else if(type == Mage) weight[i] = 10000;
//        else if(type == Ranger) weight[i] = 2000;
//        else if(type == Healer) weight[i] = 500;
//        else if(type == Factory) weight[i] = 10000;
//        else if(type == Rocket) weight[i] = 700;
//    }
//    vector<int> tmp(get_random_indices(weight));
//    for(int i = 0; i < tmp.size(); i++)
//    {
//        if(bc_GameController_can_attack(gc, id, bc_Unit_id(nearby_enemies[tmp[i]])))
//        {
//            bc_GameController_attack(gc, id, bc_Unit_id(nearby_enemies[tmp[i]]));
//            return 1;
//        }
//    }
//    return 0;
    if(!bc_GameController_is_attack_ready(gc, id)) return 0;
    int max_priority = -1, min_heallth = 300, att_id;
    for(int i = 0; i < nearby_enemies.size(); i++)
    {
        int ene_id = bc_Unit_id(nearby_enemies[i]);
        if(!bc_GameController_can_sense_unit(gc, ene_id)) continue;
        if(!bc_GameController_can_attack(gc, id, ene_id)) continue;;
        Ptr<bc_Unit> enemy(bc_GameController_unit(gc, ene_id));
        bc_UnitType type = bc_Unit_unit_type(enemy);
        int priority = -1, health = bc_Unit_health(enemy);
        if(type == Worker) priority = 0;
        else if(type == Knight) priority = 1;
        else if(type == Mage) priority = 6;
        else if(type == Ranger) priority = 4;
        else if(type == Healer) priority = 2;
        else if(type == Factory) priority = 5;
        else if(type == Rocket) priority = 3;
        if(priority > max_priority) max_priority = priority, min_heallth = health, att_id = ene_id;
        else if(priority == max_priority && health < min_heallth) min_heallth = health, att_id = ene_id;
    }
    if(max_priority != -1)
    {
        bc_GameController_attack(gc, id, att_id);
        check_errors("Attacking");
        return 1;
    }
    check_errors("Attacking");
    return 0;
}

bool try_attack(int id, int loc, int dist)
{
    if(!bc_GameController_is_attack_ready(gc, id)) return 0;
    vector<Ptr<bc_Unit>> nearby_enemies;
    get_nearby_enemies(loc, dist, nearby_enemies);
    return try_attack(id, loc, dist, nearby_enemies);
}

bool try_javelin(int id, int loc, int dist)
{
    if(!bc_GameController_is_javelin_ready(gc, id)) return 0;
    vector<Ptr<bc_Unit>> nearby_enemies;
    get_nearby_enemies(loc, dist, nearby_enemies);
    vector<int> weight(nearby_enemies.size());
    for(int i = 0; i < nearby_enemies.size(); i++)
    {
        bc_UnitType type = bc_Unit_unit_type(nearby_enemies[i]);
        if(type == Worker) weight[i] = 1;
        else if(type == Knight) weight[i] = 100;
        else if(type == Mage) weight[i] = 10000;
        else if(type == Ranger) weight[i] = 2000;
        else if(type == Healer) weight[i] = 500;
        else if(type == Factory) weight[i] = 10000;
        else if(type == Rocket) weight[i] = 700;
    }
    vector<int> tmp(get_random_indices(weight));
    for(int i = 0; i < tmp.size(); i++)
    {
        if(bc_GameController_can_javelin(gc, id, bc_Unit_id(nearby_enemies[tmp[i]])))
        {
            bc_GameController_javelin(gc, id, bc_Unit_id(nearby_enemies[tmp[i]]));
            check_errors("Javelin");
            return 1;
        }
    }
    check_errors("Javelin");
    return 0;
}

void get_nearby_teammates(int loc, int dist, vector<Ptr<bc_Unit>>& nearby_teammates)
{
    nearby_teammates.clear();
    int sq = 0; while((sq+1)*(sq+1) <= dist) sq++;
    int now_x = loc%map_width[my_Planet], now_y = loc/map_width[my_Planet];
    auto first = find_by_x(teammates.begin(), teammates.end(), now_x-sq);
    auto last = find_by_x(teammates.begin(), teammates.end(), now_x+sq+1);
    if(last-first < 22*dist/7) for(auto it = first; it != last; it++)
    {
        int x = it->first.first, y = it->first.second;
        if((x-now_x)*(x-now_x)+(y-now_y)*(y-now_y) <= dist) nearby_teammates.push_back(it->second);
    }
    else for(int i = 0; i <= sq; i++) for(int neg_x = -1; neg_x < (i?3:1); neg_x+=2)
        {
            if(now_x+neg_x*i < 0 || now_x+neg_x*i >= map_width[my_Planet]) continue;
            for(int j = 0; j*j+i*i <= dist; j++) for(int neg_y = -1; neg_y < (j?3:1); neg_y += 2)
            {
                if(now_y+neg_y*j < 0 || now_y+neg_y*j >= map_height[my_Planet]) continue;
                int new_loc = now_x+neg_x*i + (now_y+neg_y*j)*map_width[my_Planet];
                if(units[new_loc] && bc_Unit_team(units[new_loc]) == my_Team) nearby_teammates.push_back(units[new_loc]);
            }
        }
}

bool try_heal(int id, int loc, int dist)
{
    if(!bc_GameController_is_heal_ready(gc, id)) return 0;
    vector<Ptr<bc_Unit>> nearby_teammates;
    get_nearby_teammates(loc, dist, nearby_teammates);
    vector<int> weight(nearby_teammates.size());
    for(int i = 0; i < nearby_teammates.size(); i++)
    {
        bc_UnitType type = bc_Unit_unit_type(nearby_teammates[i]);
        int mate_id = bc_Unit_id(nearby_teammates[i]);
        if(!is_robot(type)) continue;
        if(!bc_GameController_can_sense_unit(gc, mate_id)) continue;
        if(!bc_GameController_can_heal(gc, id, mate_id)) continue;
        Ptr<bc_Unit> tmp(bc_GameController_unit(gc, mate_id));
        int lost_health = bc_Unit_max_health(tmp) - bc_Unit_health(tmp);
        if(type == Worker) weight[i] = lost_health;
        else if(type == Knight) weight[i] = lost_health*20;
        else if(type == Ranger) weight[i] = lost_health*15;
        else if(type == Mage) weight[i] = lost_health*40;
        else if(type == Healer) weight[i] = lost_health*3;
    }
    vector<int> tmp(get_random_indices(weight));
    for(auto i:tmp)
    {
        if(bc_GameController_can_heal(gc, id, bc_Unit_id(nearby_teammates[i])))
        {
            bc_GameController_heal(gc, id, bc_Unit_id(nearby_teammates[i]));
            return 1;
        }
    }
    return 0;
}

void try_load(int id, int loc)
{
    for(int i = 0; i < 8; i++)
    {
        if(out_of_bound(loc, i) || !units[loc+go(i)]) continue;
        int id2 = bc_Unit_id(units[loc+go(i)]);
        bc_UnitType type = bc_Unit_unit_type(units[loc+go(i)]);
        if(!is_robot(type)) continue;
        if(bc_GameController_can_load(gc, id, id2))
        {
            bc_GameController_load(gc, id, id2);
            built_rocket[id].second++;
            if(type == Worker) built_rocket[id].first++;
            units[loc+go(i)] = Ptr<bc_Unit>();
            on_rocket.insert(id2);
        }
    }
}

void launch(int id, int loc)
{
    vector<int> weight(map_width[Mars]*map_height[Mars]);
    for(int i = 0; i < weight.size(); i++) weight[i] = min(connected_square_num[i], built_rocket[i].second+1);
    vector<int> tmp(get_random_indices(weight));
    int x = tmp[0]%map_width[Mars], y = tmp[0]/map_width[Mars];
    if(!bc_GameController_can_launch_rocket(gc, id, new_bc_MapLocation(Mars, x, y))) return;
    bc_GameController_launch_rocket(gc, id, new_bc_MapLocation(Mars, x, y));
    built_rocket.erase(id);
    units[loc] = Ptr<bc_Unit>();
    waiting_for_rocket = 0;
    check_errors("Launching");
}

bool random_walk(int id, int loc, vector<int>& weight)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    while(weight.size() < 9) weight.push_back(1);
    vector<int> tmp(get_random_indices(weight));
    for(auto dir:tmp)
    {
        if(bc_Direction(dir) == Center) return 0;
        if(bc_GameController_can_move(gc, id, bc_Direction(dir)))
        {
            bc_GameController_move_robot(gc, id, bc_Direction(dir));
            units[loc] = Ptr<bc_Unit>();
            units[loc+go(dir)] = bc_GameController_unit(gc, id);
            return 1;
        }
    }
    return 0;
}

bool random_walk(int id, int loc)
{
    vector<int> tmp(9,1);
    return random_walk(id, loc, tmp);
}

int get_walk_to_direction(int id, int now_loc, int dest_loc)
{
    if(shortest_distance[now_loc][dest_loc] == -1) return 8;
    if(now_loc == dest_loc) return 8;
    unsigned int shortest_dist = -1, nearest_dir;
    vector<int> tmp(8); for(int i = 0; i < 8; i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    for(auto i:tmp)
    {
        if(out_of_bound(now_loc, i)) continue;
        if(shortest_distance[now_loc+go(i)][dest_loc] < shortest_dist)//So this direction is the preferred one
            shortest_dist = shortest_distance[now_loc+go(i)][dest_loc], nearest_dir = i;
    }
    shortest_distance[now_loc][dest_loc] = shortest_distance[dest_loc][now_loc] = max(shortest_distance[now_loc][dest_loc],shortest_dist+1);
    vector<int> adjust({0,1,7});
    if(rand()%2) swap(adjust[0], adjust[1]);
    for(auto adjust_dir:adjust)
    {
        int dir = (nearest_dir+adjust_dir)%8;
        if(bc_GameController_can_move(gc, id, bc_Direction(dir)))
            return dir;
    }
    return 8;
}

bool walk_to(int id, int now_loc, int dest_loc)
{
    if(shortest_distance[now_loc][dest_loc] == -1) return 0;
    if(now_loc == dest_loc) return 0;
    unsigned int shortest_dist = -1, nearest_dir;
    vector<int> tmp(8); for(int i = 0; i < 8; i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    for(auto i:tmp)
    {
        if(out_of_bound(now_loc, i)) continue;
        if(shortest_distance[now_loc+go(i)][dest_loc] < shortest_dist)//So this direction is the preferred one
            shortest_dist = shortest_distance[now_loc+go(i)][dest_loc], nearest_dir = i;
    }
    shortest_distance[now_loc][dest_loc] = shortest_distance[dest_loc][now_loc] = max(shortest_distance[now_loc][dest_loc],shortest_dist+1);
    vector<int> adjust({0,1,7});
    if(rand()%2) swap(adjust[0], adjust[1]);
    for(auto adjust_dir:adjust)
    {
        int dir = (nearest_dir+adjust_dir)%8;
        if(bc_GameController_can_move(gc, id, bc_Direction(dir)))
        {
            bc_GameController_move_robot(gc, id, bc_Direction(dir));
            units[now_loc] = Ptr<bc_Unit>();
            units[now_loc+go(dir)] = bc_GameController_unit(gc, id);
            return 1;
        }
        Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, (now_loc+go(dir))%map_width[my_Planet], (now_loc+go(dir))/map_width[my_Planet]));
        if(!out_of_bound(now_loc, dir) && bc_GameController_has_unit_at_location(gc, tmp))
        {
            Ptr<bc_Unit> blocking_unit(bc_GameController_sense_unit_at_location(gc, tmp));
            if(bc_Unit_team(blocking_unit) == my_Team && is_robot(bc_Unit_unit_type(blocking_unit)))
                poked_direction[bc_Unit_id(blocking_unit)] = dir;
        }
    }
    return 0;
}

bool poked(int id, int now_loc)
{
    if(!bc_GameController_is_move_ready(gc, id) || poked_direction[id] == -1) return 0;
    int pdir = poked_direction[id];
    poked_direction[id] = -1;
    assert(0<=pdir && pdir<8);
    if(bc_GameController_can_move(gc, id, bc_Direction(pdir)))
    {
        if(get_total_damage(now_loc+go(pdir)) > get_total_damage(now_loc)) return 0;
        bc_GameController_move_robot(gc, id, bc_Direction(pdir));
        units[now_loc] = Ptr<bc_Unit>();
        units[now_loc+go(pdir)] = bc_GameController_unit(gc, id);
        return 1;
    }
    Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, (now_loc+go(pdir))%map_width[my_Planet], (now_loc+go(pdir))/map_width[my_Planet]));
    if(!out_of_bound(now_loc, pdir) && bc_GameController_has_unit_at_location(gc, tmp))
    {
        Ptr<bc_Unit> blocking_unit(bc_GameController_sense_unit_at_location(gc, tmp));
        if(bc_Unit_team(blocking_unit) == my_Team && is_robot(bc_Unit_unit_type(blocking_unit)))
            poked_direction[bc_Unit_id(blocking_unit)] = pdir;
        vector<int> adjust({1,7,2,6});
        if(rand()%2) swap(adjust[0], adjust[1]); if(rand()%2) swap(adjust[2], adjust[3]);
        for(auto adjust_dir:adjust)
        {
            int dir = (pdir+adjust_dir)%8;
            if(bc_GameController_can_move(gc, id, bc_Direction(dir)))
            {
                if(get_total_damage(now_loc+go(dir)) > get_total_damage(now_loc)) return 0;
                bc_GameController_move_robot(gc, id, bc_Direction(dir));
                units[now_loc] = Ptr<bc_Unit>();
                units[now_loc+go(dir)] = bc_GameController_unit(gc, id);
                check_errors("Poking");
                return 1;
            }
        }
        check_errors("Poking");
        return 0;
    }
    vector<int> try_dir({1,7,2,6,3,5,4});
    if(rand()%2) swap(try_dir[0], try_dir[1]); if(rand()%2) swap(try_dir[2], try_dir[3]); if(rand()%2) swap(try_dir[4], try_dir[5]);
    for(int i = 0; i < 7; i++)
    {
        int dir = (pdir+try_dir[i])%8;
        if(!out_of_bound(now_loc, dir) && passable[my_Planet][now_loc+go(dir)]) poked_direction[id] = dir;
        return poked(id, now_loc);
    }
    return poked(id, now_loc);
}


bool walk_to_enemy(int id, int now_loc, bc_UnitType type)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    unsigned int nearest_dist = -1, nearest_factory_dist = -1;
    int nearest_loc, nearest_factory_loc;
    vector<int> tmp(enemies.size()); for(int i = 0; i < tmp.size(); i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    int h = map_height[my_Planet], w = map_width[my_Planet];
    int nowx = now_loc%w, nowy = now_loc/w;
    bool ignore_worker = get_total_damage(now_loc) && type == Knight;
    for(int i = 0; i < enemies.size(); i++)
    {
        int x = enemies[i].first.first, y = enemies[i].first.second;
        if(ignore_worker && bc_Unit_unit_type(enemies[i].second) == Worker) continue;
        if(shortest_distance[y*w+x][nowy*w+nowx] < nearest_dist)
        {
            nearest_dist = shortest_distance[y*w+x][nowy*w+nowx];
            nearest_loc = y*w+x;
        }
        if(bc_Unit_unit_type(enemies[i].second) == Factory && shortest_distance[y*w+x][nowy*w+nowx] < nearest_factory_dist)
        {
            nearest_factory_dist = shortest_distance[y*w+x][nowy*w+nowx];
            nearest_factory_loc = y*w+x;
        }
    }
    if(nearest_dist == -1) return 0;
    if(type == Knight && nearest_factory_dist < 2*nearest_dist) return walk_to(id, now_loc, nearest_factory_loc);
    else if((type == Mage && nearest_dist <= 3) || (type == Ranger && nearest_dist <= 5))
    {
        for(int i = 0; i < 8; i++)
            if(!out_of_bound(now_loc, i) && shortest_distance[now_loc+go(i)][nearest_loc] > nearest_dist)
                if(walk_to(id, now_loc, now_loc+go(i))) return 1;
        return 0;
    }
    return walk_to(id, now_loc, nearest_loc);
}

bool walk_to_heal(int id, int now_loc, int round)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    int h = map_height[my_Planet], w = map_width[my_Planet];
    if(get_total_damage(now_loc))
    {
        int lowest_damage = 2147483647, best_dir;
        for(int i = 0; i < 8; i++) if(!out_of_bound(now_loc, i) && passable[my_Planet][now_loc+go(i)])
        {
            int new_loc = now_loc+go(i);
            int new_x = new_loc%w, new_y = new_loc/w;
            if(get_total_damage(new_loc) < lowest_damage) lowest_damage = get_total_damage(new_loc), best_dir = i;
        }
        return walk_to(id, now_loc, now_loc+go(best_dir));
    }
    int max_willingness = 0;
    int dest_loc;
    vector<int> tmp(teammates.size()); for(int i = 0; i < tmp.size(); i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    int nowx = now_loc%w, nowy = now_loc/w;
    for(int i = 0; i < teammates.size(); i++)
    {
        if(bc_Unit_unit_type(teammates[i].second) != Ranger) continue;
        int x = teammates[i].first.first, y = teammates[i].first.second;
        int lost_health = bc_Unit_max_health(teammates[i].second)-bc_Unit_health(teammates[i].second);
        int id = bc_Unit_id(teammates[i].second);
        if(shortest_distance[y*w+x][nowy*w+nowx] == -1) continue;
        int willingness = lost_health*1000/max(((int)shortest_distance[y*w+x][nowy*w+nowx])-4, 1);
        if(willingness > max_willingness)
        {
            max_willingness = willingness;
            dest_loc = y*w+x;
        }
    }
    if(max_willingness <= 15000) return walk_to_enemy(id, now_loc, Healer);
    int dest_x = dest_loc%w, dest_y = dest_loc/w;
    if((dest_x-nowx)*(dest_x-nowx)+(dest_y-nowy)*(dest_y-nowy) > 20)
        return walk_to(id, now_loc, dest_loc);
    return 0;
}

bool walk_to_opposite(int id, int loc, int start_loc)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    int dest_loc = corresponding_point(start_loc, symmetry);
    return walk_to(id, loc, dest_loc);
}

bool walk_to_build(int id, int now_loc, int dest_loc)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    unsigned int nearest_dist = -1;
    int nearest_loc;
    for(int i = 0; i < 8; i++)
    {
        if(out_of_bound(dest_loc, i)) continue;
        if(shortest_distance[now_loc][dest_loc+go(i)] < nearest_dist)
        {
            nearest_dist = shortest_distance[now_loc][dest_loc+go(i)];
            nearest_loc = dest_loc+go(i);
        }
    }
    if(nearest_dist == -1) return 0;
    return walk_to(id, now_loc, nearest_loc);
}


//This simulates the attraction and the repulsion by all other units with gravity force
pair<double, double> worker_gravity_force(int now_loc)
{
    double force_x = 0, force_y = 0;
    int h = map_height[my_Planet], w = map_width[my_Planet];
    int now_x = now_loc%w, now_y = now_loc/w;

    for(int i = 0; i < h*w; i++)
    {
//        Calculate the "mass": positive => attractive, negative => repel
        double mass = 0;
        int x = i%w, y = i/w;
        if(now_x == x && now_y == y) continue;
        mass += karbonite[i];
        if(!passable[my_Planet][i] || units[i]) mass -= 5;
        if(units[i])
        {
            bc_UnitType type = bc_Unit_unit_type(units[i]);
            if(bc_Unit_team(units[i]) == my_Team)
            {
                if(!is_robot(type))
                {
                    if(bc_Unit_structure_is_built(units[i]))
                    {
                        double lost_health = bc_Unit_max_health(units[i])-bc_Unit_health(units[i]);
                        double dmass = lost_health/bc_Unit_max_health(units[i])*(bc_UnitType_blueprint_cost(type)+50)-20;
                        mass += dmass;
                    }
                    else mass += bc_UnitType_blueprint_cost(type)-40;
                }
            }
            else if(is_robot(type) && type != Worker && type != Healer) mass -= 60;
        }
//        Calculate the force
        double difx = x-now_x, dify = y-now_y;
        double mag = sqrt(difx*difx + dify*dify);
        force_x += difx/mag/mag/mag*mass, force_y += dify/mag/mag/mag*mass;
    }
    double center_difx = map_width[my_Planet]/2 - now_x, center_dify = map_height[my_Planet]/2-now_y;
//    Prevent worker from sticking to the wall
    force_x += center_difx+rand()%7-3, force_y += center_dify+rand()%7-3; //Random term
    return make_pair(force_x, force_y);
}

pair<double, double> healer_gravity_force(int now_loc)
{
    double force_x = 0, force_y = 0;
    int h = map_height[my_Planet], w = map_width[my_Planet];
    int now_x = now_loc%w, now_y = now_loc/w;

    for(int i = 0; i < h*w; i++)
    {
//        Calculate the "mass": positive => attractive, negative => repel
        double mass = 0;
        int x = i%w, y = i/w;
        if(now_x == x && now_y == y) continue;
        if(!passable[my_Planet][i] || units[i]) mass -= 3;
        if(units[i])
        {
            bc_UnitType type = bc_Unit_unit_type(units[i]);
            if(bc_Unit_team(units[i]) == my_Team)
            {
                if(is_robot(type))
                {
                    if(type == Worker) mass++;
                    if(type == Knight) mass += 20;
                    if(type == Ranger) mass += 15;
                    if(type == Mage) mass += 40;
                    if(type == Healer) mass += 3;
                }
                else mass -= 30;
            }
            else if(is_robot(type) && type != Worker && type != Healer) mass -= 50;
        }
//        Calculate the force
        double difx = x-now_x, dify = y-now_y;
        double mag = sqrt(difx*difx + dify*dify);
        force_x += difx/mag/mag/mag*mass, force_y += dify/mag/mag/mag*mass;
    }
    double center_difx = map_width[my_Planet]/2 - now_x, center_dify = map_height[my_Planet]/2-now_y;
//    Prevent worker from sticking to the wall
    force_x += center_difx+rand()%7-3, force_y += center_dify+rand()%7-3; //Random term
    return make_pair(force_x, force_y);
}

bool can_move(int id)
{
    if(!bc_GameController_is_move_ready(gc,id)) return 0;
    for(int i = 0; i < 8; i++) if(bc_GameController_can_move(gc, id, bc_Direction(i))) return 1;
    return 0;
}

bool try_walk_to_rocket(int id, bc_UnitType type, int now_loc)
{
    if(!can_move(id)) return 0;
    if(target_rocket[id] != -1)
    {
        int tmp = target_rocket[id];
        target_rocket[id] = -1;
        return walk_to(id, now_loc, tmp);
    }
    int dest_loc = -1;
    unsigned int nearest_dist = -1;
    for(auto rocket:need_bot_rocket)
    {
        if(shortest_distance[now_loc][rocket.second] < nearest_dist)
        {
            nearest_dist = shortest_distance[now_loc][rocket.second];
            dest_loc = rocket.second;
        }
    }
    if(dest_loc == -1) return 0;
    return walk_to(id, now_loc, dest_loc);
}

int get_next_chunk(int label)
{
    if(chunk_next_chunk_to_karbonite[label] != -1) return chunk_next_chunk_to_karbonite[label];
    chunk_next_chunk_to_karbonite[label] = -2;
//    Dijkstra
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;
    vector<bool> visited(chunk_num);
    pq.push(make_pair(0, label));
    double max_willingness = 0;
    while(pq.size())
    {
        pair<int,int> tmp = pq.top(); pq.pop();
        int now_label = tmp.second, dist = tmp.first;
        if(visited[now_label]) continue;
        visited[now_label] = 1;
        if(now_label != label)
        {
            double willingness = ((double)chunk_karbonite[now_label])/chunk_size[now_label]/dist;
            if(willingness > max_willingness) max_willingness = willingness, chunk_next_chunk_to_karbonite[label] = now_label;
        }
        for(auto pr:chunk_edge[now_label]) if(!visited[pr.first])
        {
            assert(shortest_distance[chunk_rep[pr.first]][chunk_rep[now_label]] != -1);
            if(chunk_friend_fire[pr.first] < chunk_enemy_fire[pr.first]) continue;
            pq.push(make_pair(dist+shortest_distance[chunk_rep[pr.first]][chunk_rep[now_label]], pr.first));
        }
    }
    return chunk_next_chunk_to_karbonite[label];
}

bool walk_to_harvest(int id, int now_loc)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    if(karbonite[now_loc]) return 0;
    for(int i = 0; i < 8; i++) if(!out_of_bound(now_loc, i) && karbonite[now_loc+go(i)] && passable[my_Planet][now_loc+go(i)])
        return walk_to(id, now_loc, now_loc+go(i));
    unsigned int nearest_dist = -1;
    int dest_loc, now_x = now_loc%map_width[my_Planet], now_y = now_loc/map_width[my_Planet];
    for(int i = -3; i <= 3; i++)
    {
        if(now_x+i < 0 || now_x+i >= map_width[my_Planet]) continue;
        for(int j = -3; j <= 3; j++)
        {
            if(now_y+j < 0 || now_y+j >= map_height[my_Planet]) continue;
            int loc = (now_x+i) + (now_y+j)*map_width[my_Planet];
            if(karbonite[loc]) for(int k = 0; k < 9; k++)
                if(!out_of_bound(loc, k) && shortest_distance[now_loc][loc+go(k)] < nearest_dist)
                    dest_loc = loc+go(k), nearest_dist = shortest_distance[now_loc][loc+go(k)];
        }
    }
    if(nearest_dist != (unsigned int)-1) return walk_to(id, now_loc, dest_loc);
    int label = chunk_label[now_loc];
    int dest_label = get_next_chunk(label);
    if(dest_label != -2) return walk_to(id, now_loc, chunk_rep[dest_label]);
    return 0;
}

int get_replicate_direction(int id, int now_loc)
{
    if(karbonite[now_loc]) return 8;
    for(int i = 0; i < 8; i++) if(!out_of_bound(now_loc, i) && karbonite[now_loc+go(i)] && passable[my_Planet][now_loc+go(i)])
        return get_walk_to_direction(id, now_loc, now_loc+go(i));
    unsigned int nearest_dist = -1;
    int dest_loc, now_x = now_loc%map_width[my_Planet], now_y = now_loc/map_width[my_Planet];
    for(int i = -3; i <= 3; i++)
    {
        if(now_x+i < 0 || now_x+i >= map_width[my_Planet]) continue;
        for(int j = -3; j <= 3; j++)
        {
            if(now_y+j < 0 || now_y+j >= map_height[my_Planet]) continue;
            int loc = (now_x+i) + (now_y+j)*map_width[my_Planet];
            if(karbonite[loc]) for(int k = 0; k < 9; k++)
                if(!out_of_bound(loc, k) && shortest_distance[now_loc][loc+go(k)] < nearest_dist)
                    dest_loc = loc+go(k), nearest_dist = shortest_distance[now_loc][loc+go(k)];
        }
    }
    if(nearest_dist != (unsigned int)-1) return get_walk_to_direction(id, now_loc, dest_loc);
    int label = chunk_label[now_loc];
    int dest_label = get_next_chunk(label);
    if(dest_label != -2) return get_walk_to_direction(id, now_loc, chunk_rep[dest_label]);
    return 8;
}

pair<bool,int> replicate_to_harvest(int id, int loc)
{
    int dir = get_replicate_direction(id, loc);
    if(dir == 8 || !bc_GameController_can_replicate(gc, id, bc_Direction(dir))) return try_replicate(id, loc);
    bc_GameController_replicate(gc, id, bc_Direction(dir));
    int new_loc = loc+go(dir);
    Ptr<bc_MapLocation> tmpmloc(new_bc_MapLocation(my_Planet, new_loc%map_width[my_Planet], new_loc/map_width[my_Planet]));
    Ptr<bc_Unit> tmp(bc_GameController_sense_unit_at_location(gc, tmpmloc));
    units[new_loc] = tmp;
    if(my_Planet == Earth) left_worker_num[connected_comp_label[loc]]--;
    return make_pair(1, new_loc);
}

bool explore_unknown(int id, int now_loc)
{
    if(!bc_GameController_is_move_ready(gc, id) || invisible_loc == -1) return 0;
    return walk_to(id, now_loc, invisible_loc);
}

Ptr<bc_MapLocation> get_mlocation(Ptr<bc_Unit>& unit)
{
    Ptr<bc_Location> tmp(bc_Unit_location(unit));
    if(!bc_Location_is_on_map(tmp))
    {
        return Ptr<bc_MapLocation>();
    }
    return bc_Location_map_location(tmp);
}

void update_unit_location(int id, Ptr<bc_Unit>& unit, int& now_loc)
{
    unit = bc_GameController_unit(gc, id);
    Ptr<bc_MapLocation> now_mlocation = get_mlocation(unit);//Remember to Update
    now_loc = bc_MapLocation_x_get(now_mlocation)+bc_MapLocation_y_get(now_mlocation)*map_width[my_Planet];
}

bool should_replicate(int id, int now_loc, int karb, int round)
{
    if(my_Planet == Mars)
    {
        if(round >= 750) return 1;
        int karbo = 0;
        for(int i = 0; i < 9; i++)
            if(!out_of_bound(now_loc, i))
            {
                karbo += karbonite[now_loc+go(i)];
                if(karbo >= 90) return 1;
            }
        return 0;
    }
    if(left_worker_num[connected_comp_label[now_loc]] <= 0) return 0;
    if(should_build_rocket && karb-60 < ((teammates.size()+7)/8-building_rocket.size()-built_rocket.size()) * 150) return 0;
    return 1;
}

bool try_anihilate(int id)
{
    if(!bc_GameController_can_sense_unit(gc, id)) return 1;
    for(auto teammate:teammates)
    {
        bc_UnitType type = bc_Unit_unit_type(teammate.second);
        if(type == Knight || type == Ranger)
        {
            int attacker_id = bc_Unit_id(teammate.second);
            if(bc_GameController_is_attack_ready(gc, attacker_id) && bc_GameController_can_attack(gc, attacker_id, id))
            {
                bc_GameController_attack(gc, attacker_id, id);
                if(!bc_GameController_can_sense_unit(gc, id)) return 1;
            }
        }
    }
    return 0;
}

bool try_snipe(int id)
{
    if(!bc_GameController_is_begin_snipe_ready(gc, id)) return 0;
    int max_priority = -1, att_index;
    vector<int> tmp(enemies.size()); for(int i = 0; i < enemies.size(); i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    for(int i = 0; i < enemies.size(); i++)
    {
        int priority;
        int ene_id = bc_Unit_id(enemies[tmp[i]].second);
        if(!bc_GameController_can_sense_unit(gc, ene_id)) continue;
        bc_UnitType type = bc_Unit_unit_type(enemies[tmp[i]].second);
        if(type == Worker) priority = 0;
        else if(type == Knight) priority = 2;
        else if(type == Mage) priority = 3;
        else if(type == Ranger) priority = 4;
        else if(type == Healer) priority = 5;
        else if(type == Factory) priority = 6;
        else if(type == Rocket) priority = 1;
        if(priority > max_priority) max_priority = priority, att_index = tmp[i];
    }
    if(max_priority != -1)
    {
        Ptr<bc_MapLocation> maploc(new_bc_MapLocation(my_Planet, enemies[att_index].first.first, enemies[att_index].first.second));
        if(bc_GameController_can_begin_snipe(gc, id, maploc))
        {
            bc_GameController_begin_snipe(gc, id, maploc);
            return 1;
        }
    }
    return 0;
}

int dfs_for_karbonite(int loc, vector<bool>& visited)
{
    if(visited[loc]) return 0;
    int ans = karbonite[loc];
    visited[loc] = 1;
    connected_comp_label[loc] = left_worker_num.size();
    if(!passable[my_Planet][loc]) return ans;
    for(int i = 0; i < 8; i++) if(!out_of_bound(loc, i))
        ans += dfs_for_karbonite(loc+go(i), visited);
    return ans;
}

void to_be_harvest(int loc)
{
    for(int i = 0; i < 9; i++) if(!out_of_bound(loc, i) && !being_harvested[loc+go(i)])
    {
        being_harvested[loc+go(i)] = 1;
        if(passable[my_Planet][loc+go(i)]) chunk_karbonite[chunk_label[loc+go(i)]] -= karbonite[loc+go(i)];
        karbonite[loc+go(i)] = 0;
    }
}

bool should_knight_rush(int now_loc)
{
    unsigned int nearest_dist = -1, nearest_factory_dist = -1;
    int nearest_loc, nearest_factory_loc;
    vector<int> tmp(enemies.size()); for(int i = 0; i < tmp.size(); i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    int h = map_height[my_Planet], w = map_width[my_Planet];
    for(int i = 0; i < enemies.size(); i++)
    {
        int x = enemies[i].first.first, y = enemies[i].first.second;
        if(bc_Unit_unit_type(enemies[i].second) == Worker || bc_Unit_unit_type(enemies[i].second) == Healer) continue;
        for(int j = 0; j < 8; j++)
        {
            if(out_of_bound(now_loc, j)) continue;
            if(bc_Unit_unit_type(enemies[i].second) != Factory && bc_Unit_unit_type(enemies[i].second) != Knight) continue;
            if(shortest_distance[y*w+x][now_loc+go(j)] < nearest_dist)
            {
                nearest_dist = shortest_distance[y*w+x][now_loc+go(j)];
                nearest_loc = y*w+x;
            }
            if(bc_Unit_unit_type(enemies[i].second) == Factory && shortest_distance[y*w+x][now_loc+go(j)] < nearest_factory_dist)
            {
                nearest_factory_dist = shortest_distance[y*w+x][now_loc+go(j)];
                nearest_factory_loc = y*w+x;
            }
        }
    }
    if(nearest_dist <= 10) return 1;
    return 0;
}

bool try_overcharge(int id)
{
    for(auto healer_id:overcharge_healer) if(bc_GameController_can_overcharge(gc, healer_id, id))
    {
        bc_GameController_overcharge(gc, healer_id, id);
        overcharge_healer.erase(healer_id);
        return 1;
    }
    return 0;
}

int main() {
    printf("Meow Starting\n");

    srand(7122);
    gen.seed(7122);
    printf("Connecting to manager...\n");

    gc = new_bc_GameController();
    check_errors("Connecting");
    planetmap[0] = bc_GameController_starting_map(gc, Earth);
    planetmap[1] = bc_GameController_starting_map(gc, Mars);
    my_Planet = bc_GameController_planet(gc);
    my_Team = bc_GameController_team(gc);
    cout<<"I am team "<<my_Team<<','<<my_Planet<<endl;
    int print_round = 1001;

    organize_map_info();
    cout<<"HEY"<<endl;
    orbit_pattern = bc_GameController_orbit_pattern(gc);
    cout<<"I guess this map has a symmetry of "<<symmetry<<endl;

    int start_loc = -1;
    fill(poked_direction, poked_direction+65536, -1);
    fill(worker_build_target, worker_build_target+65536, -1);
    fill(target_rocket, target_rocket+65536, -1);

    if(my_Planet == Earth)
    {
        bc_GameController_queue_research(gc, Worker);
        bc_GameController_queue_research(gc, Knight);
        bc_GameController_queue_research(gc, Healer);
        bc_GameController_queue_research(gc, Healer);
        bc_GameController_queue_research(gc, Ranger);
        bc_GameController_queue_research(gc, Rocket);
        bc_GameController_queue_research(gc, Healer);
        bc_GameController_queue_research(gc, Ranger);
        bc_GameController_queue_research(gc, Ranger);
        for(int i = 0; i < 3; i++) bc_GameController_queue_research(gc, Mage);
        for(int i = 0; i < 3; i++) bc_GameController_queue_research(gc, Knight);
        for(int i = 0; i < 3; i++) bc_GameController_queue_research(gc, Worker);
    }

    while (true)
    {
        uint32_t round = bc_GameController_round(gc);
        printf("Round: %d\n", round);
        int karb = bc_GameController_karbonite(gc);
        cout<<"Karbonite: "<<karb<<endl;//For debug

        Ptr<bc_ResearchInfo> research_info(bc_GameController_research_info(gc));
        if(bc_ResearchInfo_get_level(research_info, Rocket) && my_Planet == Earth) can_build_rocket = 1;

        typecount.clear(); typecount.resize(7);
        fill(visible, visible+2500, 0);
        fill(units, units+2500, Ptr<bc_Unit>());
        fill(chunk_enemy_fire.begin(), chunk_enemy_fire.end(), 0);
        fill(chunk_friend_fire.begin(), chunk_friend_fire.end(), 0);
        fill(chunk_next_chunk_to_karbonite.begin(), chunk_next_chunk_to_karbonite.end(), -1);
        try_blueprint_rocket = success_bluesprint_rocket = 0;
        fill(enemies_max_total_damage, enemies_max_total_damage+2500, -1); //-1 = uncalculated
        teammates.clear();
        enemies.clear();
        need_bot_rocket.clear();
        invisible_loc = -1;
        int idle_num = 0;

        if(my_Planet == Mars) // Update karbonite levels from asteroid strikes
        {
            if(bc_AsteroidPattern_has_asteroid(asteroid_pattern, round))
            {
                Ptr<bc_AsteroidStrike> tmp = bc_AsteroidPattern_asteroid(asteroid_pattern, round);
                Ptr<bc_MapLocation> tmpmloc = bc_AsteroidStrike_location_get(tmp);
                int loc = bc_MapLocation_x_get(tmpmloc)+bc_MapLocation_y_get(tmpmloc)*map_width[my_Planet];
                int karbo = bc_AsteroidStrike_karbonite_get(tmp);
                if(passable[Mars][loc])
                    chunk_karbonite[chunk_label[loc]] += karbo;
                karbonite[loc] += karbo;
            }
        }
        overcharge_healer.clear();

         // For each square...
        for(int i = 0; i < map_width[my_Planet]; i++) for(int j = 0; j < map_height[my_Planet]; j++)
        {
            // Location at this square
            Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, i, j));

             // If we can see this square
            if(bc_GameController_can_sense_location(gc, tmp))
            {
                int loc = i+j*map_width[my_Planet]; // Encode loc
                visible[loc] = 1;
                 // SOMETHING WITH KARBONITE LEVELS
                if(!being_harvested[loc])
                {
                    if(passable[my_Planet][loc])
                        chunk_karbonite[chunk_label[loc]] += bc_GameController_karbonite_at(gc, tmp)-karbonite[loc];
                    karbonite[loc] = bc_GameController_karbonite_at(gc, tmp);
                }
                else
                {
                    if(passable[my_Planet][loc])
                        chunk_karbonite[chunk_label[loc]] -= karbonite[loc];
                    karbonite[loc] = 0;
                }
                being_harvested[loc] = 0;
                // If there is a unit on this square
                if(bc_GameController_has_unit_at_location(gc, tmp))
                {
                    units[loc] = bc_GameController_sense_unit_at_location(gc, tmp);
                    bc_UnitType type = bc_Unit_unit_type(units[loc]);
                    // If the unit is on my team
                    if(bc_Unit_team(units[loc]) == my_Team)
                    {
                        // We have a list of team units
                        teammates.push_back(make_pair(make_pair(i,j),units[loc]));
                        // and a counter for each type
                        typecount[bc_Unit_unit_type(units[loc])]++;
                        if(is_robot(type))
                        {
                            if(type == Healer)
                            {
                                 // Chunk_enemy_fire is used in some "willingness" decisions later
                                chunk_enemy_fire[chunk_label[loc]] -= 10;
                                // We have a list of healers ready to use overcharge
                                if(bc_GameController_is_overcharge_ready(gc, bc_Unit_id(units[loc])))
                                    overcharge_healer.insert(bc_Unit_id(units[loc]));
                            }
                            // Chunk_friend_fire is also used in these "willingness" decisions
                            else chunk_friend_fire[chunk_label[loc]] += bc_Unit_damage(units[loc]);
                        }
                        // We also have lists of live rockets and factories
                        else if(type == Rocket) alive_rockets.insert(bc_Unit_id(units[loc]));
                        else alive_factories.insert(bc_Unit_id(units[loc]));
                    }
                    // If the unit is an enemy
                    else
                    {
                        // We have a list of enemies
                        enemies.push_back(make_pair(make_pair(i,j),units[loc]));
                        if(is_robot(type))
                        {
                            // Enemy healers reduce our effective firepower in the chunk
                            if(type == Healer) chunk_friend_fire[chunk_label[loc]] -= 10;
                            // Other enemy robots add to their effective firepower in the chunk
                            else chunk_friend_fire[chunk_label[loc]] += bc_Unit_damage(units[loc]);
                        }
                    }
                }
            }
            // If the spot is not visible
            else if(passable[my_Planet][i+j*map_width[my_Planet]]) invisible_loc = i+j*map_width[my_Planet];
        }

        if(enemies.size()) have_enemy_round++;
        else have_enemy_round = 0;

        should_build_rocket = round >= 500;
        if(!should_build_rocket && can_build_rocket)
        {
            for(auto teammate:teammates)
            {
                bc_UnitType type = bc_Unit_unit_type(teammate.second);
                if(!is_robot(type) || type == Worker || type == Healer) continue;
                int now_loc = teammate.first.first+teammate.first.second*map_width[my_Planet];
                vector<Ptr<bc_Unit>> nearby_teammates;
                get_nearby_enemies(now_loc, 10, nearby_teammates);
                bool is_idle = 1;
                for(auto nearby_teammate:nearby_teammates)
                    if(round - last_attack_round[bc_Unit_id(nearby_teammate)] < 20)
                    {
                        is_idle = 0;
                        break;
                    }
                int id = bc_Unit_id(teammate.second);
                if(rocket_just_blueprinted) idle_round[id] = 0;
                if(is_idle) idle_round[id]++;
                else idle_round[id] = 0;
                if(idle_round[id] >= 75) idle_num++;
            }
            if(idle_num >= max((typecount[1]+typecount[2]+typecount[3])*3/5, 8)) should_build_rocket = 1;
        }
        rocket_just_blueprinted = 0;
        if(should_build_rocket) waiting_for_rocket++;

        if(enemies.size()) first_enemy = 1;

        vector<int> high_priority;
//        If cannot launch rockets for a long time, then kill teammates for spaces
        if(waiting_for_rocket >= 150 && my_Planet == Earth)
        {
            for(auto teammate:teammates) if(bc_Unit_unit_type(teammate.second) == Worker || bc_Unit_unit_type(teammate.second) == Factory)
            {
                high_priority.push_back(bc_Unit_id(teammate.second));
                int now_loc = teammate.first.first+teammate.first.second*map_width[my_Planet];
                for(int i = 0; i < 8; i++) if(!out_of_bound(now_loc, i) && units[now_loc+go(i)])
                {
                    Ptr<bc_Unit>& tmp = units[now_loc+go(i)];
                    bc_UnitType type = bc_Unit_unit_type(tmp);
                    if(!is_robot(type) || type == Worker) continue;
                    try_anihilate(bc_Unit_id(tmp));
                }
            }
        }
//        For the early stage strategy. Determine the starting point.
        if(round == 1 && my_Planet == Earth)
        {
            vector<int> willingness(teammates.size()); //willingness to take teammates[i] as start_loc
            vector<int> loc(teammates.size());
            vector<bool> visited(map_width[my_Planet]*map_height[my_Planet]);
            for(int i = 0; i < teammates.size(); i++)
            {
                loc[i] = teammates[i].first.first+teammates[i].first.second*map_width[my_Planet];
                if(!visited[loc[i]])
                {
                    int expected_karbonite = dfs_for_karbonite(loc[i], visited);
                    left_worker_num.push_back(expected_karbonite/(120+expected_karbonite/60)-1);
                }
                else left_worker_num[connected_comp_label[loc[i]]]--;
            }
            for(int i = 0; i < teammates.size(); i++)
            {
                for(int j = 0; j < teammates.size(); j++)
                    willingness[i] += 1000-min(1000u, shortest_distance[loc[i]][loc[j]]);
                willingness[i] += chunk_size[chunk_label[loc[i]]]/5;
                willingness[i] += 10-min(10u, shortest_distance[loc[i]][chunk_rep[chunk_label[loc[i]]]]);
            }
            int max_willingess = 0, index;
            for(int i = 0; i < teammates.size(); i++)
                if(max_willingess < willingness[i]) max_willingess = willingness[i], start_loc = loc[i], index = i;
//            high_priority.push_back(index);
        }
        else if(karb >= 200 && !should_build_rocket && my_Planet == Earth)
        {
            int max_willingness = -100000, index;
            for(int i = 0; i < teammates.size(); i++)
            {
                if(bc_Unit_unit_type(teammates[i].second) != Worker) continue;
                if(can_harvest(bc_Unit_id(teammates[i].second))) continue;
                int willingness = 0, loc = teammates[i].first.first + teammates[i].first.second*map_width[my_Planet];
                willingness -= 100*get_total_damage(loc);
                willingness += distance_to_wall[loc]*10;
                willingness += chunk_size[chunk_label[loc]];
                willingness += min(0, chunk_friend_fire[chunk_label[loc]]-chunk_enemy_fire[chunk_label[loc]]);
                for(int i = 0; i < 8; i++) if(out_of_bound(loc, i) || !passable[my_Planet][loc+go(i)]) willingness -= 10;
                if(willingness > max_willingness) max_willingness = willingness, index = i;
            }
            high_priority.push_back(index);
        }

        vector<int> to_erase;
        for(auto rocket:building_rocket) if(!alive_rockets.count(rocket.first)) to_erase.push_back(rocket.first);
        for(auto erase_id:to_erase) building_rocket.erase(erase_id);
        to_erase.clear();
        for(auto rocket:built_rocket) if(!alive_rockets.count(rocket.first)) to_erase.push_back(rocket.first);
        for(auto erase_id:to_erase) built_rocket.erase(erase_id), built_rocket_location.erase(erase_id);
        alive_rockets.clear(); to_erase.clear();
        for(auto factory:my_factories) if(!alive_factories.count(factory.first))
        {
            to_erase.push_back(factory.first);
            passable[my_Planet][factory.second] = 1;
            for(int i = 0; i < map_height[my_Planet]*map_width[my_Planet]; i++)
            {
                if(!passable[my_Planet][i]) continue;
                for(int k = 0; k < 8; k++)
                    if(!out_of_bound(factory.second, k) && shortest_distance[i][factory.second+go(k)] != -1)
                        shortest_distance[factory.second][i] = shortest_distance[i][factory.second]
                        = min(shortest_distance[i][factory.second+go(k)]+1, shortest_distance[i][factory.second]);
            }
        }
        for(auto erase_id:to_erase) my_factories.erase(erase_id);
        alive_factories.clear();

        for(auto rocket:building_rocket)
        {
            bool working = 0;
            int id = rocket.first, now_loc = rocket.second;
            for(int i = 0; i < 8; i++) if(!out_of_bound(now_loc, i) && units[now_loc+go(i)])
                if(bc_Unit_unit_type(units[now_loc+go(i)]) == Worker && bc_Unit_team(units[now_loc+go(i)]) == my_Team)
                    working = 1, should_stay[bc_Unit_id(units[now_loc+go(i)])] = 1;
            if(working) structure_wait_for_work_round[id] = 0;
            else
            {
                int wait_round = ++structure_wait_for_work_round[id];
                vector<Ptr<bc_Unit>> nearby_teammates;
                int x = now_loc%map_width[my_Planet], y = now_loc/map_width[my_Planet];
                for(int i = -wait_round-1; i <= wait_round+1; i++)
                {
                    if(x+i < 0 || x+i >= map_width[my_Planet]) continue;
                    for(int j = -wait_round-1; j <= wait_round+1; j++)
                    {
                        if(y+j < 0 || y+j >= map_height[my_Planet]) continue;
                        int new_loc = (x+i) + (y+j)*map_width[my_Planet];
                        if(units[new_loc] && bc_Unit_unit_type(units[new_loc]) == Worker && bc_Unit_team(units[new_loc]) == my_Team)
                            worker_build_target[bc_Unit_id(units[new_loc])] = now_loc;
                    }
                }
            }
        }

        for(auto rocket:built_rocket_location)
        {
            if(built_rocket[rocket.first].second == 8) continue;
            int cnt = 0, now_loc = rocket.second;
            int x = now_loc%map_width[my_Planet], y = now_loc/map_width[my_Planet];
            for(int i = -3; i <= 3; i++)
            {
                if(x+i < 0 || x+i >= map_width[my_Planet]) continue;
                for(int j = -3; j <= 3; j++)
                {
                    if(y+j < 0 || y+j >= map_height[my_Planet]) continue;
                    int new_loc = (x+i) + (y+j)*map_width[my_Planet];
                    if(units[new_loc] && bc_Unit_unit_type(units[new_loc]) != Worker && bc_Unit_team(units[new_loc]) == my_Team)
                        target_rocket[bc_Unit_id(units[new_loc])] = now_loc, cnt++;
                }
            }
            if(cnt + built_rocket[rocket.first].second < 8) need_bot_rocket[rocket.first] = rocket.second;
        }

        for(auto factory:my_factories)
        {
            if(!bc_Unit_structure_is_built(units[factory.second]))
            {
                bool working = 1;
                int id = factory.first, now_loc = factory.second;
                for(int i = 0; i < 8; i++) if(!out_of_bound(now_loc, i) && passable[my_Planet][now_loc+go(i)])
                    if(!units[now_loc+go(i)] || bc_Unit_unit_type(units[now_loc+go(i)]) != Worker
                       || bc_Unit_team(units[now_loc+go(i)]) != my_Team)
                        working = 0;
                if(working) structure_wait_for_work_round[id] = 0;
                else
                {
                    int wait_round = ++structure_wait_for_work_round[id];
                    vector<Ptr<bc_Unit>> nearby_teammates;
                    int x = now_loc%map_width[my_Planet], y = now_loc/map_width[my_Planet];
                    for(int i = max(-5,-wait_round-1); i <= min(5,wait_round+1); i++)
                    {
                        if(x+i < 0 || x+i >= map_width[my_Planet]) continue;
                        for(int j = max(-5,-wait_round-1); j <= min(5,wait_round+1); j++)
                        {
                            if(y+j < 0 || y+j >= map_height[my_Planet]) continue;
                            int new_loc = (x+i) + (y+j)*map_width[my_Planet];
                            if(units[new_loc] && bc_Unit_unit_type(units[new_loc]) == Worker && bc_Unit_team(units[new_loc]) == my_Team)
                                worker_build_target[bc_Unit_id(units[new_loc])] = now_loc;
                        }
                    }
                }
            }
            else
            {
                vector<Ptr<bc_Unit>> nearby_enemies;
                get_nearby_enemies(factory.second, 50, nearby_enemies);
                int id = factory.first, now_loc = factory.second;
                int x = now_loc%map_width[my_Planet], y = now_loc/map_width[my_Planet];
                if(nearby_enemies.size())
                {
                    high_priority.push_back(factory.first);
                    int ene_id = bc_Unit_id(nearby_enemies[0]);
                    int ene_loc;
                    update_unit_location(ene_id, nearby_enemies[0], ene_loc);
                    int ene_x = ene_loc%map_width[my_Planet], ene_y = ene_loc/map_width[my_Planet];
                    int index = 0;
                    if(ene_x > x) index += 6;
                    else if(ene_x == x) index += 3;
                    if(ene_y > y) index += 2;
                    else if(ene_y == y) index++;
                    int tmp[9] = {5, 6, 7, 4, 8, 0, 3, 2, 1};
                    enemy_direction[id] = tmp[index];
                }
                int lost_health = 300 - bc_Unit_health(units[factory.second]);
                for(int i = -lost_health/20; i <= lost_health/20; i++)
                {
                    if(x+i < 0 || x+i >= map_width[my_Planet]) continue;
                    for(int j = -lost_health/20; j <= lost_health/20; j++)
                    {
                        if(y+j < 0 || y+j >= map_height[my_Planet]) continue;
                        int new_loc = (x+i) + (y+j)*map_width[my_Planet];
                        if(units[new_loc] && bc_Unit_unit_type(units[new_loc]) == Worker && bc_Unit_team(units[new_loc]) == my_Team)
                            worker_build_target[bc_Unit_id(units[new_loc])] = now_loc;
                    }
                }
            }
        }

        vector<int> tmprandom(teammates.size()); for(int i = 0; i < teammates.size(); i++) tmprandom[i] = i;
        random_shuffle(tmprandom.begin(), tmprandom.end());
        for(int ii = 0; ii < high_priority.size(); ii++) for(int i = 0; i < teammates.size(); i++)
            if(high_priority[ii] == tmprandom[i])
                swap(tmprandom[i], tmprandom[ii]);
        for(int ii = 0; ii < teammates.size(); ii++)
        {
            if(round >= print_round) cout<<"INIT"<<endl;
            Ptr<bc_Unit> unit(teammates[tmprandom[ii]].second);
            int id = bc_Unit_id(unit);
            if(on_rocket.count(id)) continue;
            int now_loc = teammates[tmprandom[ii]].first.first + teammates[tmprandom[ii]].first.second*map_width[my_Planet];
            bc_UnitType type = bc_Unit_unit_type(unit);
            if(round >= print_round) cout<<"START POKED"<<endl;
            if(poked(id, now_loc)) update_unit_location(id, unit, now_loc);
            if(round >= print_round) cout<<"END POKED"<<endl;
            if(type != Worker && (round >= 425 || round-last_attack_round[id] >= 50))
                if(try_walk_to_rocket(id, type, now_loc)) update_unit_location(id, unit, now_loc);
            bool should_overcharge = 0;
//            Here is what a worker will do
            if(type == Worker)
            {
//                A worker tries stuff in this order:
//                1. Harvest 2. Build 3. Replicate 4. Blueprint rockets 5. Blueprint factories 6. Repair
//                TODO: Make every attempt into a function, and change the order based on some other numbers
                if(round >= print_round) cout<<"Worker"<<endl;
                bool done = 0, dontmove = 0;
                if(worker_build_target[id] != -1 && !should_stay[id])
                {
                    dontmove = 1;
                    if(bc_GameController_is_move_ready(gc, id))
                        if(walk_to_build(id, now_loc, worker_build_target[id])) update_unit_location(id, unit, now_loc);
                }
                should_stay[id] = 0;
                if(!done && try_build(id, now_loc)) done = dontmove = 1; //Stay still to continue building
                if(!done && try_repair(id, now_loc)) done = dontmove = 1; //Stay still to continue repairing
                if(can_harvest(id)) to_be_harvest(now_loc);
                else if(!dontmove && (my_Planet == Mars || round <= 150) && walk_to_harvest(id, now_loc)) update_unit_location(id, unit, now_loc);
                if(!done && should_replicate(id, now_loc, karb, round) && (my_Planet == Mars || typecount[0] <= passable_count[my_Planet]/20))
                {
                    auto tmprep = replicate_to_harvest(id, now_loc);
                    if(tmprep.first)
                    {
                        now_loc = tmprep.second, unit = units[now_loc], id = bc_Unit_id(unit);
                        if(walk_to_harvest(id, now_loc)) update_unit_location(id, unit, now_loc);
                        if(try_harvest(id)) done = dontmove = 1, to_be_harvest(now_loc);
                    }
                    else if(try_harvest(id)) done = dontmove = 1;
                }
                else if(!done && try_harvest(id)) done = dontmove = 1;
                if(!done && building_rocket.size()+built_rocket.size() < (teammates.size()+7)/8
                   && should_build_rocket && bc_GameController_karbonite(gc) >= 150 && worker_build_target[id] == -1)
                    if(try_blueprint(id, now_loc, Rocket)) done = dontmove = 1; //Stay still to build rocket
                if(!done && worker_build_target[id] == -1 && (!should_build_rocket || my_factories.size() <= 4))
                    if(try_blueprint(id, now_loc, Factory)) done = dontmove = 1; //Stay still to build factory
                if(!dontmove && can_move(id))
                {
                    if(!try_walk_to_rocket(id, type, now_loc))
                    {
                        auto f = worker_gravity_force(now_loc);
                        //cout<<"FORCE"<<f.first<<' '<<f.second<<endl;
                        vector<int> walk_weight;
                        for(int i = 0; i < 8; i++)
                        {
                            int difx = (map_width[my_Planet]+go(i)+1)%map_width[my_Planet] - 1;
                            int dify = (go(i)+1)/map_width[my_Planet];
                            walk_weight.push_back(max((int)(difx*f.first + dify*f.second), 0));
                        }
                        walk_weight.push_back(20);
                        random_walk(id, now_loc, walk_weight);
                    }
                }
                if(worker_build_target[id] == -1 || !units[worker_build_target[id]] || bc_Unit_structure_is_built(units[worker_build_target[id]]))
                    worker_build_target[id] = -1;
                worker_is_building[id] = 0;
                check_errors("Worker's turn");
            }
            else if(type == Factory)
            {
                if(round >= print_round) cout<<"Factory"<<endl;
                if(!bc_Unit_structure_is_built(unit)) continue;
                while(1)
                    if(!try_unload(id, now_loc).first) break;
                if(!bc_GameController_can_produce_robot(gc, id, Knight)) continue;
                if(should_build_rocket && !need_worker && karb-40 < ((teammates.size()+7)/8-building_rocket.size()-built_rocket.size())*200 && typecount[0]) continue;
                vector<int> weight({0,0,10,0,min(3, (int)round/20)});
                if(should_build_rocket) weight[0] = 1, weight[1] = 3, weight[3] = 3;
                if(round >= 500) for(int i = 0; i < 5; i++) weight[i] = ((i%2)?1:0);
                if(should_build_rocket && (typecount[0] <= building_rocket.size() || need_worker)) for(int i = 0; i < 5; i++) weight[i] = (i?0:1);
                if(should_knight_rush(now_loc)) for(int i = 0; i < 5; i++) weight[i] = (i==1?1:0);
                try_produce(id, weight);
                check_errors("Factory's turn");
            }
            else if(type == Knight)
            {
                if(round >= print_round) cout<<"Knight"<<endl;
                bool done = 0;
                if(try_javelin(id, now_loc, 10)) last_attack_round[id] = round, should_overcharge = 1;
                if(try_attack(id, now_loc, 2)) last_attack_round[id] = round, should_overcharge = 1;
                else
                {
                    if(my_Planet == Earth && !first_enemy) walk_to_opposite(id, now_loc, start_loc);
                    else walk_to_enemy(id, now_loc, type);
                    update_unit_location(id, unit, now_loc);
                    if(try_javelin(id, now_loc, 10)) last_attack_round[id] = round, should_overcharge = 1;
                    if(try_attack(id, now_loc, 2)) last_attack_round[id] = round, should_overcharge = 1;
                }
                check_errors("Knight's turn");
            }
            else if(type == Healer)
            {
                if(round >= print_round) cout<<"Healer"<<endl;
                check_errors("INIT");
                if(try_heal(id, now_loc, 30)) last_attack_round[id] = round;
                else if(walk_to_heal(id, now_loc, round))
                {
                    update_unit_location(id, unit, now_loc);
                    if(try_heal(id, now_loc, 30)) last_attack_round[id] = round;
                }
//                if(can_move(id))
//                {
//                    auto f = healer_gravity_force(now_loc);
//                    vector<int> walk_weight;
//                    for(int i = 0; i < 8; i++)
//                    {
//                        int difx = (map_width[my_Planet]+go(i)+1)%map_width[my_Planet] - 1;
//                        int dify = (go(i)+1)/map_width[my_Planet];
//                        walk_weight.push_back(max((int)(difx*f.first + dify*f.second), 0));
//                    }
//                    walk_weight.push_back(20);
//                    if(random_walk(id, now_loc, walk_weight))
//                    {
//                        update_unit_location(id, unit, now_loc);
//                        if(try_heal(id, now_loc, 30)) last_attack_round[id] = round;
//                    }
//                }
                check_errors("Healer's turn");
            }
            else if(type == Rocket)
            {
                if(round >= print_round) cout<<"Rocket"<<endl;
                if(my_Planet == Earth)
                {
                    if(bc_Unit_structure_is_built(unit))
                    {
                        building_rocket.erase(id);
                        if(!built_rocket.count(id)) built_rocket[id] = make_pair(0,0);
                        if(!built_rocket_location.count(id))built_rocket_location[id] = now_loc;
                        try_load(id, now_loc);
                        int garrison_num = built_rocket[id].second;
                        if((garrison_num && bc_Unit_health(unit) < 200)
                           || (garrison_num == min(8, max_connected_num-1) && bc_OrbitPattern_duration(orbit_pattern, round) <= bc_OrbitPattern_duration(orbit_pattern, round+1))
                           || round == 749) launch(id, now_loc);
                    }
                    else building_rocket[id] = now_loc;
                }
                else
                {
                    while(1)
                        if(!try_unload(id, now_loc).first) break;
                }
                check_errors("Rocket's turn");
            }
            else if(type == Mage)
            {
                if(round >= print_round) cout<<"Mage"<<endl;
                if(try_attack(id, now_loc, 30)) last_attack_round[id] = round, should_overcharge = 1;
                if(my_Planet == Earth && !first_enemy) walk_to_opposite(id, now_loc, start_loc);
                else if(enemies.size())
                {
                    if(!get_total_damage(now_loc) && walk_to_enemy(id, now_loc, type))
                    {
                        update_unit_location(id, unit, now_loc);
                        if(try_attack(id, now_loc, 30)) last_attack_round[id] = round, should_overcharge = 1;
                    }
                }
                else if(invisible_loc != -1 && (!(id%10) || my_Planet == Mars)) explore_unknown(id, now_loc);
                check_errors("Mage's turn");
            }
            else if(type == Ranger)
            {
                if(round >= print_round) cout<<"Ranger"<<endl;
                if(try_attack(id, now_loc, 50)) last_attack_round[id] = round, should_overcharge = 1;
                if(bc_GameController_is_begin_snipe_ready(gc, id) && round - last_attack_round[id] >= 10 && round - last_snipe_round[id] > 5 && have_enemy_round >= 10)
                    should_overcharge = 1;
                if(round - last_attack_round[id] >= 10 && round - last_snipe_round[id] > 5 && have_enemy_round >= 10 && try_snipe(id))
                    last_snipe_round[id] = round, should_overcharge = 0;
                else if(my_Planet == Earth && !first_enemy && typecount[2] >= 8) walk_to_opposite(id, now_loc, start_loc);
                else if(enemies.size())
                {
                    if(!get_total_damage(now_loc) && walk_to_enemy(id, now_loc, type))
                    {
                        update_unit_location(id, unit, now_loc);
                        if(try_attack(id, now_loc, 50)) last_attack_round[id] = round, should_overcharge = 1;
                    }
                }
                else if(invisible_loc != -1 && !(id%10)) explore_unknown(id, now_loc);
                check_errors("Ranger's turn");
            }
            if(should_overcharge && try_overcharge(id)) ii--;
        }
        need_worker = try_blueprint_rocket && !success_bluesprint_rocket;
        on_rocket.clear();
        bc_GameController_next_turn(gc);
    }
}
