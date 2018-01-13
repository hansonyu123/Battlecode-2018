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

/// Redundant for now.
//struct My_Unit
//{
//    Ptr<bc_Unit> unit;
//    uint32_t health;
//    bc_UnitType type;
//    Ptr<bc_Location> location;
//    My_Unit():unit(), health(), type(Worker), location(){};
//    My_Unit(uint16_t id):unit(bc_GameController_unit(gc, id))
//    {
//        health = bc_Unit_health(unit);
//        type = bc_Unit_unit_type(unit);
//        location = bc_Unit_location(unit);
//    }
//};

bc_Planet my_Planet; //Current planet
bc_Team my_Team; //Current team
int map_height[2], map_width[2]; //The size of the Earth map and the Mars map.
vector<bool> passable[2]; //The current map
int passable_count[2];
unsigned int shortest_distance[2500][2500]; //Shortest distance between squares. (x,y) corresponds to y*w+h.
Ptr<bc_PlanetMap> planetmap[2]; //The maps.
default_random_engine gen; //C++11 random engine

bool out_of_bound(int loc, int dir) //Test if going in the direction dir from (loc%w, loc/w) will go out the map
{
    int h = map_height[my_Planet], w = map_width[my_Planet];
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

int go(int dir) //The number added to go in the direction dir
{
    int w = map_width[my_Planet];
    int a[9] = {w, w+1, 1, -w+1, -w, -w-1, -1, w-1, 0};
    return a[dir];
}

void organize_map_info() //Organize all the map information
{
    cout<<"Start organizing"<<endl;
    for(int i = 0; i < 2; i++)
    {
        map_height[i] = bc_PlanetMap_height_get(planetmap[i]);
        map_width[i] = bc_PlanetMap_width_get(planetmap[i]);
        passable[i].resize(map_height[i]*map_width[i]);
        for(int x = 0; x < map_width[i]; x++) for(int y = 0; y < map_height[i]; y++) //Read in the map
        {
            Ptr<bc_MapLocation> tmp(new_bc_MapLocation(bc_Planet(i), x, y));
            passable_count[i] += passable[i][map_width[i]*y+x] = bc_PlanetMap_is_passable_terrain_at(planetmap[i], tmp);
        }
    }
    int h = map_height[my_Planet], w = map_width[my_Planet];
    vector<bool>& p = passable[my_Planet];
    for(int i = 0; i < h*w; i++) for(int j = 0; j < h*w; j++) shortest_distance[i][j] = (i==j)?0:-1;// -1 = Very large
    cout<<"Start BFS"<<endl;
    for(int i = 0; i < h*w; i++) // BFS.
    {
        if(!p[i]) continue;
        queue<int> q; q.push(i);
        while(q.size())
        {
            int j = q.front(); q.pop();
            for(int k = 0; k < 8; k++)
                if(!out_of_bound(j, k) && shortest_distance[i][j+go(k)] == -1 && p[j+go(k)])
                {
                    shortest_distance[i][j+go(k)] = shortest_distance[i][j]+1;
                    q.push(j+go(k));
                }
        }
    }
    for(int i = 0; i < h*w; i++) //For debug uses.
    {
        cout<<shortest_distance[2*w+2][i]<<' ';
        if(!(i+1)%w) cout<<endl;
    }
    cout<<"End organizing"<<endl;
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

bool try_harvest(int id, vector<Ptr<bc_MapLocation>>& v, vector<Ptr<bc_Unit>>& nearby_units, Ptr<bc_MapLocation> now_mlocation)
{
    for(int i = 0; i < v.size(); i++)
        if(bc_GameController_can_harvest(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i])))
        {
            bc_GameController_harvest(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i]));
            return 1;
        }
    return 0;
}

bool try_build(int id, vector<Ptr<bc_MapLocation>>& v, vector<Ptr<bc_Unit>>& nearby_units)
{
    for(int i = 0; i < v.size(); i++)
        if(nearby_units[i] && !is_robot(bc_Unit_unit_type(nearby_units[i]))
            && bc_GameController_can_build(gc, id, bc_Unit_id(nearby_units[i])))
        {
            bc_GameController_build(gc, id, bc_Unit_id(nearby_units[i]));
            return 1;
        }
    return 0;
}

bool try_replicate(int id, vector<Ptr<bc_MapLocation>>& v, Ptr<bc_MapLocation> now_mlocation)
{
    for(int i = 0; i < v.size(); i++)
        if(bc_GameController_can_replicate(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i])))
        {
            bc_GameController_replicate(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i]));
            return 1;
        }
    return 0;
}

bool try_blueprint(int id, vector<Ptr<bc_MapLocation>>& v, Ptr<bc_MapLocation> now_mlocation, bc_UnitType structure)
{
    for(int i = 0; i < v.size(); i++)
        if(bc_GameController_can_blueprint(gc, id, structure, bc_MapLocation_direction_to(now_mlocation, v[i])))
        {
            bc_GameController_blueprint(gc, id, structure, bc_MapLocation_direction_to(now_mlocation, v[i]));
            return 1;
        }
    return 0;
}

bool try_repair(int id, vector<Ptr<bc_MapLocation>>& v, vector<Ptr<bc_Unit>>& nearby_units)
{
    for(int i = 0; i < v.size(); i++)
        if(nearby_units[i] && !is_robot(bc_Unit_unit_type(nearby_units[i]))
           && bc_GameController_can_repair(gc, id, bc_Unit_id(nearby_units[i]))
           && bc_Unit_max_health(nearby_units[i]) != bc_Unit_health(nearby_units[i]))
       {
           bc_GameController_repair(gc, id, bc_Unit_id(nearby_units[i]));
           return 1;
       }
    return 0;
}

bool try_unload(int id)
{
    int tmp[8]; for(int i = 0; i < 8; i++) tmp[i] = i;
    random_shuffle(tmp, tmp+8);
    for(int i = 0; i < 8; i++) if(bc_GameController_can_unload(gc, id, bc_Direction(tmp[i])))
    {
        bc_GameController_unload(gc, id, bc_Direction(tmp[i]));
        return 1;
    }
    return 0;
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

bool try_attack(int id, vector<Ptr<bc_Unit>>& nearby_units)
{
    if(!bc_GameController_is_attack_ready(gc, id)) return 0;
    vector<int> tmp(nearby_units.size()); for(int i = 0; i < tmp.size(); i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    for(int i = 0; i < nearby_units.size(); i++)
    {
        if(nearby_units[tmp[i]] && bc_Unit_team(nearby_units[tmp[i]]) == my_Team) continue;
        if(bc_GameController_can_attack(gc, id, bc_Unit_id(nearby_units[tmp[i]])))
        {
            bc_GameController_attack(gc, id, bc_Unit_id(nearby_units[tmp[i]]));
            return 1;
        }
    }
    return 0;
}

//Not working now
//bool try_javelin(int id, vector<Ptr<bc_Unit>>& nearby_units)
//{
//    if(!bc_GameController_is_javelin_ready(gc, id)) return 0;
//    vector<int> tmp(nearby_units.size()); for(int i = 0; i < tmp.size(); i++) tmp[i] = i;
//    random_shuffle(tmp.begin(), tmp.end());
//    for(int i = 0; i < nearby_units.size(); i++)
//    {
//        if(nearby_units[tmp[i]] && bc_Unit_team(nearby_units[tmp[i]]) == my_Team) continue;
//        if(bc_GameController_can_javelin(gc, id, bc_Unit_id(nearby_units[tmp[i]])))
//        {
//            bc_GameController_javelin(gc, id, bc_Unit_id(nearby_units[tmp[i]]));
//            return 1;
//        }
//    }
//}

bool random_walk(int id, vector<int>& weight)
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
            return 1;
        }
    }
    return 0;
}

bool random_walk(int id)
{
    vector<int> tmp(9,1);
    return random_walk(id, tmp);
}

map<int, int> not_free; //For Bug pathing. not_free[id] = (destination_loc, last_direction);

bool walk_to(int id, Ptr<bc_MapLocation>& now_mlocation, Ptr<bc_MapLocation> destination)
{
    int nowx = bc_MapLocation_x_get(now_mlocation);
    int nowy = bc_MapLocation_y_get(now_mlocation);
    int h = map_height[my_Planet], w = map_width[my_Planet];
    int now_loc = nowx+nowy*w;
    int x = bc_MapLocation_x_get(destination);
    int y = bc_MapLocation_y_get(destination);
    int loc = x+y*w;
    if(not_free.count(id))
    {
        int clockwise = 2*((id+1)%2)-1; //counter-clockwise or clockwise, depends on id
        int last_dir = not_free[id];
        int j = (last_dir-2*clockwise+8)%8;
        while(1)
        {
            if(bc_GameController_can_move(gc, id, bc_Direction(j)))
            {
                if(shortest_distance[now_loc][loc] > shortest_distance[now_loc+go(j)][loc])
                    not_free.erase(id);
                else not_free[id] = j;
                bc_GameController_move_robot(gc, id, bc_Direction(j));
                return 1;
            }
            j = (j+8+clockwise)%8;
            if(j == (last_dir-2*clockwise+8)%8) return 0;
        }
    }
    for(int i = 0; i < 8; i++)
    {
        if(out_of_bound(now_loc, i)) continue;
        if(shortest_distance[now_loc+go(i)][loc] == shortest_distance[now_loc][loc]-1)//So this direction is the preferred one
        {
//            Bug Pathing
            int clockwise = 2*(id%2)-1; //clockwise or counter-clockwise, depends on id
            int j = i;
            while(1)
            {
                if(bc_GameController_can_move(gc, id, bc_Direction(j)))
                {
                    if(shortest_distance[now_loc][loc] <= shortest_distance[now_loc+go(j)][loc])
                        not_free[id] = j;
                    bc_GameController_move_robot(gc, id, bc_Direction(j));
                    return 1;
                }
                j = (j+8+clockwise)%8;
                if(j == i) return 0;
            }
        }
    }
    return 0; //Although it should never reach here
}

bool walk_to_enemy(int id, Ptr<bc_MapLocation>& now_mlocation, vector<Ptr<bc_MapLocation>>& enemy_location)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    unsigned int nearest_dist = -1;
    int nearest_x, nearest_y;
    vector<int> tmp(enemy_location.size()); for(int i = 0; i < tmp.size(); i++) tmp[i] = i;
    random_shuffle(tmp.begin(), tmp.end());
    int nowx = bc_MapLocation_x_get(now_mlocation);
    int nowy = bc_MapLocation_y_get(now_mlocation);
    int h = map_height[my_Planet], w = map_width[my_Planet];
    int now_loc = nowx+nowy*w;
    for(int i = 0; i < enemy_location.size(); i++)
    {
        int x = bc_MapLocation_x_get(enemy_location[tmp[i]]);
        int y = bc_MapLocation_y_get(enemy_location[tmp[i]]);
        if(shortest_distance[y*w+x][nowy*w+nowx] < nearest_dist)
        {
            nearest_dist = shortest_distance[y*w+x][nowy*w+nowx];
            nearest_x = x, nearest_y = y;
        }
    }
    if(nearest_dist == -1) return 0;
    return walk_to(id, now_mlocation, new_bc_MapLocation(my_Planet, nearest_x, nearest_y));
}

//This simulates the attraction and the repulsion by all other units with gravity force
pair<double, double> worker_gravity_force(Ptr<bc_MapLocation> now_mlocation, vector<Ptr<bc_MapLocation>>& whole_map, vector<Ptr<bc_Unit>>& all_units)
{
    double force_x = 0, force_y = 0;
    int now_x = bc_MapLocation_x_get(now_mlocation), now_y = bc_MapLocation_y_get(now_mlocation);

    for(int i = 0; i < whole_map.size(); i++)
    {
//        Calculate the "mass": positive => attractive, negative => repel
        double mass = 0;
        int x = bc_MapLocation_x_get(whole_map[i]), y = bc_MapLocation_y_get(whole_map[i]);
        if(now_x == x && now_y == y) continue;
        mass += bc_GameController_karbonite_at(gc, whole_map[i]);
        if(!bc_GameController_is_occupiable(gc, whole_map[i])) mass -= 5;
        if(all_units[i])
        {
            bc_UnitType type = bc_Unit_unit_type(all_units[i]);
            if(bc_Unit_team(all_units[i]) == my_Team)
            {
                if(!is_robot(type))
                {
                    if(bc_Unit_structure_is_built(all_units[i]))
                    {
                        double lost_health = bc_Unit_max_health(all_units[i])-bc_Unit_health(all_units[i]);
                        double dmass = lost_health/bc_Unit_max_health(all_units[i])*(bc_UnitType_blueprint_cost(type)+50)-10;
                        mass += dmass;
                    }
                    else mass += bc_UnitType_blueprint_cost(type)-20;
                }
            }
            else if(is_robot(type) && type != Worker && type != Healer) mass -= 30;
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

int main() {
    printf("Meow Starting\n");

    srand(7122);
    printf("Connecting to manager...\n");

    gc = new_bc_GameController();
    check_errors("Connecting");
    planetmap[0] = bc_GameController_starting_map(gc, Earth);
    planetmap[1] = bc_GameController_starting_map(gc, Mars);
    my_Planet = bc_GameController_planet(gc);
    my_Team = bc_GameController_team(gc);
    cout<<"I am team "<<my_Team<<endl;
    organize_map_info();
    for(int i = 0; i < 3; i++) bc_GameController_queue_research(gc, Knight);
    bc_GameController_queue_research(gc, Rocket);
    while (true)
    {
        uint32_t round = bc_GameController_round(gc);
        printf("Round: %d\n", round);
        cout<<"Karbonite:"<<bc_GameController_karbonite(gc)<<endl;//For debug

        Ptr<bc_VecUnit> units(bc_GameController_my_units(gc));
        int len = bc_VecUnit_len(units);
        vector<int> typecount(7); //Count the number of robots of a certain type.
        for(int i = 0; i < len; i++)
        {
            Ptr<bc_Unit> unit(bc_VecUnit_index(units, i));
            typecount[bc_Unit_unit_type(unit)]++;
        }

        vector<Ptr<bc_MapLocation>> whole_map; // Save the whole visible map
        for(int i = 0; i < map_width[my_Planet]; i++) for(int j = 0; j < map_width[my_Planet]; j++)
        {
            Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, i, j));
            if(bc_GameController_can_sense_location(gc, tmp)) whole_map.push_back(tmp);
        }

        vector<Ptr<bc_Unit>> all_units(whole_map.size()); // Save the whole visible units
        for(int i = 0; i < whole_map.size(); i++) if(bc_GameController_has_unit_at_location(gc, whole_map[i]))
            all_units[i] = bc_GameController_sense_unit_at_location(gc, whole_map[i]);

        vector<Ptr<bc_MapLocation>> enemy_location; // Save the locations of enemies
        for(int i = 0; i < whole_map.size(); i++) if(all_units[i] && bc_Unit_team(all_units[i]) != my_Team)
            enemy_location.push_back(whole_map[i]);

        vector<int> tmp(len); for(int i = 0; i < len; i++) tmp[i] = i;
        random_shuffle(tmp.begin(), tmp.end());
        for(int ii = 0; ii < len; ii++)
        {
            int i = tmp[ii];
            Ptr<bc_Unit> unit(bc_VecUnit_index(units, i));
            int id = bc_Unit_id(unit);
            bc_UnitType type = bc_Unit_unit_type(unit);
            Ptr<bc_Location> tmp(bc_Unit_location(unit));
            if(!bc_Location_is_on_map(tmp)) continue; //It should always be true though
            Ptr<bc_MapLocation> now_mlocation(bc_Location_map_location(tmp));
//            Ranger's active ability is infinity, so we have to be careful by taking min with 50.
            unsigned int dist = 2;
            if(is_robot(type)) dist = min(50u, max(bc_Unit_ability_range(unit), bc_Unit_attack_range(unit)));
            Ptr<bc_VecMapLocation> vmap(bc_GameController_all_locations_within(gc, now_mlocation, dist));
            int vmap_len = bc_VecMapLocation_len(vmap);
            vector<Ptr<bc_MapLocation>> v; //Save the squares within attack/ability range
            for(int i = 0; i < vmap_len; i++) v.push_back(bc_VecMapLocation_index(vmap, i));
            vector<Ptr<bc_Unit>> nearby_units(vmap_len); //Save the units within attack/ability range
            for(int i = 0; i < vmap_len; i++) if(bc_GameController_has_unit_at_location(gc, v[i]))
                nearby_units[i] = bc_GameController_sense_unit_at_location(gc, v[i]);
//            Here is what a worker will do
            if(type == Worker)
            {
//                A worker tries stuff in this order:
//                1. Harvest 2. Build 3. Replicate 4. Blueprint rockets 5. Blueprint factories 6. Repair
//                TODO: Make every attempt into a function, and change the order based on some other numbers
                bool done = 0, dontmove = 0;
                if(!done) done = dontmove = try_harvest(id, v, nearby_units, now_mlocation); //Stay still to continue harvesting
                if(!done) done = dontmove = try_build(id, v, nearby_units); //Stay still to continue building
                if(!done && typecount[Worker] < passable_count[my_Planet]/20) //Don't want too many workers
                    done = try_replicate(id, v, now_mlocation);
                if(!done) done = dontmove = try_blueprint(id, v, now_mlocation, Rocket); //Stay still to build rocket
                if(!done) done = dontmove = try_blueprint(id, v, now_mlocation, Factory); //Stay still to build factory
                if(!done) done = dontmove = try_repair(id, v, nearby_units); //Stay still to continue repairing
                if(!dontmove)
                {
                    auto f = worker_gravity_force(now_mlocation, whole_map, all_units);
                    //cout<<"FORCE"<<f.first<<' '<<f.second<<endl;
                    vector<int> walk_weight;
                    for(int i = 0; i < 8; i++)
                    {
                        int difx = (map_width[my_Planet]+go(i)+1)%map_width[my_Planet] - 1;
                        int dify = (go(i)+1)/map_width[my_Planet];
                        walk_weight.push_back(max((int)(difx*f.first + dify*f.second), 0));
                    }
                    walk_weight.push_back(20);
                    random_walk(id, walk_weight);
                }
            }
            else if(type == Factory)
            {
                if(!bc_Unit_structure_is_built(unit)) continue;
                try_unload(id);
                vector<int> weight({0,0,0,0,0});
                try_produce(id, weight);
            }
            else if(type == Knight)
            {
                bool done = 0;
                //try_javelin(id, nearby_units);
                if(try_attack(id, nearby_units)) random_walk(id), not_free.erase(id);
                else
                {
                    walk_to_enemy(id, now_mlocation, enemy_location);
                    //Reload map
                    vmap = bc_GameController_all_locations_within(gc, now_mlocation, 10);
                    vmap_len = bc_VecMapLocation_len(vmap);
                    v.clear(); //Save the squares within attack/ability range
                    for(int i = 0; i < vmap_len; i++) v.push_back(bc_VecMapLocation_index(vmap, i));
                    nearby_units.clear(); nearby_units.resize(vmap_len);//Save the units within attack/ability range
                    for(int i = 0; i < vmap_len; i++) if(bc_GameController_has_unit_at_location(gc, v[i]))
                        nearby_units[i] = bc_GameController_sense_unit_at_location(gc, v[i]);
                    try_attack(id, nearby_units);
                }
            }
            else
            {
                try_attack(id, nearby_units);
                random_walk(id);
            }
        }
        bc_GameController_next_turn(gc);
    }
}
