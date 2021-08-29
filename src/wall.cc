// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file wall.cc
 * \brief Function implementations for the base wall classes and wall lists. */

#include "config.hh"
#include "wall.hh"

namespace voro {

/** The wall_list constructor sets up an array of pointers to wall classes. */
wall_list_2d::wall_list_2d() : walls(new wall_2d*[init_wall_size]), wep(walls), wel(walls+init_wall_size),
    current_wall_size(init_wall_size) {}

/** The wall_list destructor frees the array of pointers to the wall classes.
 */
wall_list_2d::~wall_list_2d() {
    delete [] walls;
}

/** Adds all of the walls on another wall_list to this class.
 * \param[in] wl a reference to the wall class. */
void wall_list_2d::add_wall(wall_list_2d &wl) {
    for(wall_2d **wp=wl.walls;wp<wl.wep;wp++) add_wall(*wp);
}

/** Deallocates all of the wall classes pointed to by the wall_list. */
void wall_list_2d::deallocate() {
    for(wall_2d **wp=walls;wp<wep;wp++) delete *wp;
}

/** Increases the memory allocation for the walls array. */
void wall_list_2d::increase_wall_memory() {
    current_wall_size<<=1;
    if(current_wall_size>max_wall_size)
        voro_fatal_error("Wall memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
    wall_2d **nwalls=new wall_2d*[current_wall_size],**nwp=nwalls,**wp=walls;
    while(wp<wep) *(nwp++)=*(wp++);
    delete [] walls;
    walls=nwalls;wel=walls+current_wall_size;wep=nwp;
}

/** The wall_list constructor sets up an array of pointers to wall classes. */
wall_list_3d::wall_list_3d() : walls(new wall_3d*[init_wall_size]), wep(walls), wel(walls+init_wall_size),
    current_wall_size(init_wall_size) {}

/** The wall_list destructor frees the array of pointers to the wall classes.
 */
wall_list_3d::~wall_list_3d() {
    delete [] walls;
}

/** Adds all of the walls on another wall_list to this class.
 * \param[in] wl a reference to the wall class. */
void wall_list_3d::add_wall(wall_list_3d &wl) {
    for(wall_3d **wp=wl.walls;wp<wl.wep;wp++) add_wall(*wp);
}

/** Deallocates all of the wall classes pointed to by the wall_list. */
void wall_list_3d::deallocate() {
    for(wall_3d **wp=walls;wp<wep;wp++) delete *wp;
}

/** Increases the memory allocation for the walls array. */
void wall_list_3d::increase_wall_memory() {
    current_wall_size<<=1;
    if(current_wall_size>max_wall_size)
        voro_fatal_error("Wall memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
    wall_3d **nwalls=new wall_3d*[current_wall_size],**nwp=nwalls,**wp=walls;
    while(wp<wep) *(nwp++)=*(wp++);
    delete [] walls;
    walls=nwalls;wel=walls+current_wall_size;wep=nwp;
}

}
