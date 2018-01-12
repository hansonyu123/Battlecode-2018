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
using namespace std;

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

inline void release(bc_PlanetMap *planetmap){delete_bc_PlanetMap(planetmap);}
inline void release(bc_GameController *gc){delete_bc_GameController(gc);}
inline void release(bc_Unit *unit){delete_bc_Unit(unit);}
inline void release(bc_Location *location){delete_bc_Location(location);}
inline void release(bc_MapLocation *maplocation){delete_bc_MapLocation(maplocation);}
inline void release(bc_VecUnit *units){delete_bc_VecUnit(units);}
inline void release(bc_VecMapLocation *vecmaplocation){delete_bc_VecMapLocation(vecmaplocation);}

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
        //cout<<"READY TO DELETE!"<<endl;
        if(!cntr->cnt && cntr->ptr) release(cntr->ptr), delete cntr;
        //cout<<"DELETED!"<<endl;
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

Ptr<bc_GameController> gc;

struct My_Unit
{
    Ptr<bc_Unit> unit;
    uint32_t health;
    bc_UnitType type;
    Ptr<bc_Location> location;
    My_Unit():unit(), health(), type(Worker), location(){};
    My_Unit(uint16_t id):unit(bc_GameController_unit(gc, id))
    {
        health = bc_Unit_health(unit);
        type = bc_Unit_unit_type(unit);
        location = bc_Unit_location(unit);
    }
};

bc_Planet my_Planet;
int map_height[2], map_width[2];
vector<bool> passable[2];
unsigned int shortest_distance[2500][2500];
Ptr<bc_PlanetMap> planetmap[2];

void organize_map_info()
{
    for(int i = 0; i < 2; i++)
    {
        map_height[i] = bc_PlanetMap_height_get(planetmap[i]);
        map_width[i] = bc_PlanetMap_width_get(planetmap[i]);
        passable[i].resize(map_height[i]*map_width[i]);
        for(int x = 0; x < map_width[i]; x++) for(int y = 0; y < map_height[i]; y++)
        {
            Ptr<bc_MapLocation> tmp(new_bc_MapLocation(bc_Planet(i), x, y));
            passable[i][map_width[i]*y+x] = bc_PlanetMap_is_passable_terrain_at(planetmap[i], tmp);
        }
    }
    int h = map_height[my_Planet], w = map_width[my_Planet];
    vector<bool>& p = passable[my_Planet];
    for(int i = 0; i < h*w; i++) for(int j = 0; j < h*w; j++) shortest_distance[i][j] = (i==j)?0:-1;
    for(int i = 0; i < h*w; i++)
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
    for(int i = 0; i < h*w; i++)
    {
        cout<<shortest_distance[2*w+2][i]<<' ';
        if(!(i+1)%w) cout<<endl;
    }
}

bool is_robot(bc_UnitType s)
{
    return (s != Rocket) && (s != Factory);
}

int main() {
    printf("Meow Starting\n");

    srand(7122);
    printf("Connecting to manager...\n");

    // Most methods return pointers; methods returning integers or enums are the only exception.
    gc = new_bc_GameController();
    check_errors("Connecting");
    planetmap[0] = bc_GameController_starting_map(gc, Earth);
    planetmap[1] = bc_GameController_starting_map(gc, Mars);
    organize_map_info();
    bc_GameController_queue_research(gc, Rocket);
    while (true)
    {
        uint32_t round = bc_GameController_round(gc);
        printf("Round: %d\n", round);
        cout<<"Karbonite:"<<bc_GameController_karbonite(gc)<<endl;

        Ptr<bc_VecUnit> units(bc_GameController_my_units(gc));
        int len = bc_VecUnit_len(units);
        vector<int> typecount(7);
        for(int i = 0; i < len; i++)
        {
            Ptr<bc_Unit> unit(bc_VecUnit_index(units, i));
            typecount[bc_Unit_unit_type(unit)]++;
        }

        vector<Ptr<bc_MapLocation>> whole_map;
        for(int i = 0; i < map_width[my_Planet]; i++) for(int j = 0; j < map_width[my_Planet]; j++)
        {
            Ptr<bc_MapLocation> tmp(new_bc_MapLocation(my_Planet, i, j));
            if(bc_GameController_can_sense_location(gc, tmp)) whole_map.push_back(tmp);
        }

        vector<Ptr<bc_Unit>> all_units(whole_map.size());
        for(int i = 0; i < whole_map.size(); i++) if(bc_GameController_has_unit_at_location(gc, whole_map[i]))
            all_units.push_back(bc_GameController_sense_unit_at_location(gc, whole_map[i]));

        for(int i = 0; i < len; i++)
        {
            Ptr<bc_Unit> unit(bc_VecUnit_index(units, i));
            int id = bc_Unit_id(unit);
            bc_UnitType type = bc_Unit_unit_type(unit);
            Ptr<bc_Location> tmp(bc_Unit_location(unit));
            if(!bc_Location_is_on_map(tmp)) continue;
            Ptr<bc_MapLocation> now_mlocation(bc_Location_map_location(tmp));
            unsigned int dist = 2;
            if(is_robot(type)) dist = min(50u, max(bc_Unit_ability_range(unit), bc_Unit_attack_range(unit)));
            Ptr<bc_VecMapLocation> vmap(bc_GameController_all_locations_within(gc, now_mlocation, dist));
            int vmap_len = bc_VecMapLocation_len(vmap);
            vector<Ptr<bc_MapLocation>> v;
            for(int i = 0; i < vmap_len; i++) v.push_back(bc_VecMapLocation_index(vmap, i));
            vector<Ptr<bc_Unit>> nearby_units(vmap_len);
            for(int i = 0; i < vmap_len; i++) if(bc_GameController_has_unit_at_location(gc, v[i]))
                nearby_units[i] = bc_GameController_sense_unit_at_location(gc, v[i]);
            if(type == Worker)
            {
                bool done = 0, dontmove = 0;
                if(!done) for(int i = 0; i < vmap_len; i++)
                    if(bc_GameController_can_harvest(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i])))
                    {
                        bc_GameController_harvest(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i]));
                        done = 1;
                        dontmove = 1;
                        break;
                    }
                if(!done) for(int i = 0; i < vmap_len; i++)
                    if(nearby_units[i] && (bc_Unit_unit_type(nearby_units[i]) == Factory || bc_Unit_unit_type(nearby_units[i]) == Rocket)
                       && bc_GameController_can_build(gc, id, bc_Unit_id(nearby_units[i])))
                    {
                        bc_GameController_build(gc, id, bc_Unit_id(nearby_units[i]));
                        done =1 ;
                        dontmove = 1;
                        break;
                    }
                if(!done && typecount[Worker] < map_height[my_Planet]*map_width[my_Planet]/10) for(int i = 0; i < vmap_len; i++)
                    if(bc_GameController_can_replicate(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i])))
                    {
                        bc_GameController_replicate(gc, id, bc_MapLocation_direction_to(now_mlocation, v[i]));
                        done = 1;
                        break;
                    }
                if(!done) for(int i = 0; i < vmap_len; i++)
                    if(bc_GameController_can_blueprint(gc, id, Rocket, bc_MapLocation_direction_to(now_mlocation, v[i])))
                    {
                        bc_GameController_blueprint(gc, id, Rocket, bc_MapLocation_direction_to(now_mlocation, v[i]));
                        done = 1;
                        dontmove = 1;
                        break;
                    }
                if(!done) for(int i = 0; i < vmap_len; i++)
                    if(bc_GameController_can_blueprint(gc, id, Factory, bc_MapLocation_direction_to(now_mlocation, v[i])))
                    {
                        bc_GameController_blueprint(gc, id, Factory, bc_MapLocation_direction_to(now_mlocation, v[i]));
                        done = 1;
                        dontmove = 1;
                        break;
                    }
                if(!done) for(int i = 0; i < vmap_len; i++)
                    if(nearby_units[i] && (bc_Unit_unit_type(nearby_units[i]) == Factory || bc_Unit_unit_type(nearby_units[i]) == Rocket)
                       && bc_GameController_can_repair(gc, id, bc_Unit_id(nearby_units[i])))
                   {
                       bc_GameController_repair(gc, id, bc_Unit_id(nearby_units[i]));
                       done = 1;
                       dontmove = 1;
                       break;
                   }
                if(bc_GameController_is_move_ready(gc, id)) while(!dontmove)
                {
                    int random_number = rand()%9;
                    if(bc_Direction(random_number) == Center) dontmove = 1;
                    if(bc_GameController_can_move(gc, id, bc_Direction(random_number)))
                    {
                        bc_GameController_move_robot(gc, id, bc_Direction(random_number));
                        dontmove = 1;
                    }
                }
            }
        }
        bc_GameController_next_turn(gc);
    }
}
