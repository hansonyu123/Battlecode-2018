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
unsigned int shortest_distance[2500][2500]; //Shortest distance between squares. (x,y) corresponds to y*w+h.
Ptr<bc_PlanetMap> planetmap[2]; //The maps.

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
            passable[i][map_width[i]*y+x] = bc_PlanetMap_is_passable_terrain_at(planetmap[i], tmp);
        }
    }
    int h = map_height[my_Planet], w = map_width[my_Planet];
    vector<bool>& p = passable[my_Planet];
    for(int i = 0; i < h*w; i++) for(int j = 0; j < h*w; j++) shortest_distance[i][j] = (i==j)?0:-1;// -1 = Very large
    for(int i = 0; i < h*w; i++) // BFS. Apology for the ugly code.
    {
        if(!p[i]) continue;
        queue<int> q; q.push(i);
        while(q.size())
        {
            int j = q.front(); q.pop();
            bool up = j/w != h-1, down = j/w, lft = j%w, rght = j%w != w-1;
            if(up && shortest_distance[i][j+w] == -1 && p[j+w]) shortest_distance[i][j+w] = shortest_distance[i][j]+1, q.push(j+w);
            if(up && rght && shortest_distance[i][j+w+1] == -1 && p[j+w+1]) shortest_distance[i][j+w+1] = shortest_distance[i][j]+1, q.push(j+w+1);
            if(rght && shortest_distance[i][j+1] == -1 && p[j+1]) shortest_distance[i][j+1] = shortest_distance[i][j]+1, q.push(j+1);
            if(rght && down && shortest_distance[i][j-w+1] == -1 && p[j-w+1]) shortest_distance[i][j-w+1] = shortest_distance[i][j]+1, q.push(j-w+1);
            if(down && shortest_distance[i][j-w] == -1 && p[j-w]) shortest_distance[i][j-w] = shortest_distance[i][j]+1, q.push(j-w);
            if(lft && down && shortest_distance[i][j-w-1] == -1 && p[j-w-1]) shortest_distance[i][j-w-1] = shortest_distance[i][j]+1, q.push(j-w-1);
            if(lft && shortest_distance[i][j-1] == -1 && p[j-1]) shortest_distance[i][j-1] = shortest_distance[i][j]+1, q.push(j-1);
            if(up && lft && shortest_distance[i][j+w-1] == -1 && p[j+w-1]) shortest_distance[i][j+w-1] = shortest_distance[i][j]+1, q.push(j+w-1);
        }
    }
    for(int i = 0; i < h*w; i++) //For debug uses.
    {
        cout<<shortest_distance[2*w+2][i]<<' ';
        if(!(i+1)%w) cout<<endl;
    }
}

bool is_robot(bc_UnitType s) //Check if a unit is a robot.
{
    return (s != Rocket) && (s != Factory);
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
           && bc_GameController_can_repair(gc, id, bc_Unit_id(nearby_units[i])))
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
    int sum = 0; for(auto w:weight) sum += w;
    while(sum)
    {
        int random_number = rand()%sum;
        for(int i = 0; i < 5; i++)
        {
            if(random_number >= weight[i])
            {
                random_number -= weight[i];
                continue;
            }
//            Decide to try to produce a bot of type i
            if(bc_GameController_can_produce_robot(gc, id, bc_UnitType(i)))
            {
                bc_GameController_produce_robot(gc, id, bc_UnitType(i));
                return 1;
            }
//            Failed. Set the weight of type i to zero and retry.
            sum -= weight[i];
            weight[i] = 0;
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
        if(bc_Unit_team(nearby_units[tmp[i]]) == my_Team) continue;
        if(bc_GameController_can_attack(gc, id, bc_Unit_id(nearby_units[tmp[i]])))
        {
            bc_GameController_attack(gc, id, bc_Unit_id(nearby_units[tmp[i]]));
            return 1;
        }
    }
    return 0;
}

bool random_walk(int id)
{
    if(!bc_GameController_is_move_ready(gc, id)) return 0;
    while(1)
    {
        int random_number = rand()%9; //a random direction
        if(bc_Direction(random_number) == Center) return 0;
        if(bc_GameController_can_move(gc, id, bc_Direction(random_number)))
        {
            bc_GameController_move_robot(gc, id, bc_Direction(random_number));
            return 1;
        }
    }
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
    organize_map_info();
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
            all_units.push_back(bc_GameController_sense_unit_at_location(gc, whole_map[i]));

        for(int i = 0; i < len; i++)
        {
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
                if(!done && typecount[Worker] < map_height[my_Planet]*map_width[my_Planet]/10) //Don't want too many workers
                    done = try_replicate(id, v, now_mlocation);
                if(!done) done = dontmove = try_blueprint(id, v, now_mlocation, Rocket); //Stay still to build rocket
                if(!done) done = dontmove = try_blueprint(id, v, now_mlocation, Factory); //Stay still to build factory
                if(!done) done = dontmove = try_repair(id, v, nearby_units); //Stay still to continue repairing
                if(!dontmove) random_walk(id);
            }
            else if(type == Factory)
            {
                if(!bc_Unit_structure_is_built(unit)) continue;
                try_unload(id);
                vector<int> weight({0,1,1,7,1});
                try_produce(id, weight);
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
