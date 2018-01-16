#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#ifdef MEOW
#include "bc.h"
#else
#include <bc.h>
#endif

/// See bc.h at xxx for the API you have access too.
/// Note: the API is not thread safe; don't use pthreads.
/// It's not very pretty, sorry. If you want a prettier low-level language, maybe consider rust?

// Any method in the API may set an error.
// Call check_errors() to get the most recent error.
bool check_errors() {
    /// Check if we have an error...
    if (bc_has_err()) {
        char *err;
        /// Note: this clears the current global error.
        int8_t code = bc_get_last_err(&err);
        printf("Engine error code %d: %s\n", code, err);
        bc_free_string(err);
        return true;
    } else {
        return false;
    }
}

int main() {
    printf("Meow Starting\n");

    srand(0);
    printf("Connecting to manager...\n");

    // Most methods return pointers; methods returning integers or enums are the only exception.
    bc_GameController *gc = new_bc_GameController();

    if (check_errors()) {
        // If there was an error creating gc, just die.
        printf("Failed, dying.\n");
        exit(1);
    }
    printf("Connected!\n");


    // loop through the whole game.
    while (true)
    {
        // The API is "object-oriented" - most methods take pointers to some object.
        uint32_t round = bc_GameController_round(gc);
        printf("Round: %d\n", round);

        // Note that all operations perform copies out of their data structures, returning new objects.
        // You're responsible for freeing objects.
        bc_VecUnit *units = bc_GameController_my_units(gc);
        // it's good to cache things like this. Calls into the API have around a 20ns overhead, plus the cost of
        // copying the data out of the engine. not horrible, but good to avoid more than necessary.
        int len = bc_VecUnit_len(units);
        bc_Planet my_Planet = Earth;
        if(round == 1)
        {
            if(len)
                printf("I'm on EARTH!\n");
            else
                printf("I'm on MARS!\n"), my_Planet = Mars;
        }
        fflush(stdout);
        for (int i = 0; i < len; i++) {
            // Get the current unit. This also copies.
            bc_Unit *unit = bc_VecUnit_index(units, i);

            // Calls on the controller take unit IDs for ownership reasons.
            uint16_t id = bc_Unit_id(unit);
            if (bc_GameController_can_move(gc, id, North) && bc_GameController_is_move_ready(gc, id)) {
                bc_Location *now_location = bc_Unit_location(unit);
                printf("%s\n", bc_Location_debug(now_location));
                bc_GameController_move_robot(gc, id, North);
                delete_bc_Location(now_location);
                delete_bc_Unit(unit);
                unit = bc_GameController_unit(gc, id);
                now_location = bc_Unit_location(unit);
                printf("%s\n", bc_Location_debug(now_location));
            }
            // don't want memory leaks!
            delete_bc_Unit(unit);

        }
        delete_bc_VecUnit(units);

        // this line helps the output logs make more sense by forcing output to be sent
        // to the manager.
        // it's not strictly necessary, but it helps.
        fflush(stdout);

        // pause and wait for the next turn.
        bc_GameController_next_turn(gc);
    }
}
