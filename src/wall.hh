// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file wall.hh
 * \brief Header file for the base wall classes and wall lists. */

#ifndef VOROPP_WALL_HH
#define VOROPP_WALL_HH

#include "cell_2d.hh"
#include "cell_3d.hh"

namespace voro {

/** \brief Pure virtual class from which wall objects are derived.
 *
 * This is a pure virtual class for a generic wall object. A wall object can be
 * specified by deriving a new class from this and specifying the functions. */
class wall_2d {
    public:
        virtual ~wall_2d() {};
        /** A pure virtual function for testing whether a point is
         * inside the wall object. */
        virtual bool point_inside(double x,double y) = 0;
        /** A pure virtual function for cutting a cell without
         * neighbor-tracking with a wall. */
        virtual bool cut_cell(voronoicell_2d &c,double x,double y) = 0;
        /** A pure virtual function for cutting a cell with
         * neighbor-tracking enabled with a wall. */
        virtual bool cut_cell(voronoicell_neighbor_2d &c,double x,double y) = 0;
};

/** \brief A class for storing a list of pointers to walls.
 *
 * This class stores a list of pointers to wall classes. It contains several
 * simple routines that make use of the wall classes (such as telling whether a
 * given position is inside all of the walls or not). It can be used by itself,
 * but also forms part of container_base, for associating walls with this
 * class. */
class wall_list_2d {
    public:
        /** An array holding pointers to wall objects. */
        wall_2d **walls;
        /** A pointer to the next free position to add a wall pointer. */
        wall_2d **wep;
        wall_list_2d();
        ~wall_list_2d();
        /** Adds a wall to the list.
         * \param[in] w the wall to add. */
        inline void add_wall(wall_2d *w) {
            if(wep==wel) increase_wall_memory();
            *(wep++)=w;
        }
        /** Adds a wall to the list.
         * \param[in] w a reference to the wall to add. */
        inline void add_wall(wall_2d &w) {add_wall(&w);}
        void add_wall(wall_list_2d &wl);
        /** Determines whether a given position is inside all of the walls on
         * the list.
         * \param[in] (x,y) the position to test.
         * \return True if it is inside, false if it is outside. */
        inline bool point_inside_walls(double x,double y) {
            for(wall_2d **wp=walls;wp<wep;wp++) if(!((*wp)->point_inside(x,y))) return false;
            return true;
        }
        /** Cuts a Voronoi cell by all of the walls currently on the list.
         * \param[in] c a reference to the Voronoi cell class.
         * \param[in] (x,y) the position of the cell.
         * \return True if the cell still exists, false if the cell is deleted.
         */
        template<class c_class_2d>
        bool apply_walls(c_class_2d &c,double x,double y) {
            for(wall_2d **wp=walls;wp<wep;wp++) if(!((*wp)->cut_cell(c,x,y))) return false;
            return true;
        }
        void deallocate();
    protected:
        void increase_wall_memory();
        /** A pointer to the limit of the walls array, used to determine when
         * array is full. */
        wall_2d **wel;
        /** The current amount of memory allocated for walls. */
        int current_wall_size;
};

/** \brief Pure virtual class from which wall objects are derived.
 *
 * This is a pure virtual class for a generic wall object. A wall object can be
 * specified by deriving a new class from this and specifying the functions. */
class wall_3d {
    public:
        virtual ~wall_3d() {}
        /** A pure virtual function for testing whether a point is
         * inside the wall object. */
        virtual bool point_inside(double x,double y,double z) = 0;
        /** A pure virtual function for cutting a cell without
         * neighbor-tracking with a wall. */
        virtual bool cut_cell(voronoicell_3d &c,double x,double y,double z) = 0;
        /** A pure virtual function for cutting a cell with
         * neighbor-tracking enabled with a wall. */
        virtual bool cut_cell(voronoicell_neighbor_3d &c,double x,double y,double z) = 0;
};

/** \brief A class for storing a list of pointers to walls.
 *
 * This class stores a list of pointers to wall classes. It contains several
 * simple routines that make use of the wall classes (such as telling whether a
 * given position is inside all of the walls or not). It can be used by itself,
 * but also forms part of container_base, for associating walls with this
 * class. */
class wall_list_3d {
    public:
        /** An array holding pointers to wall objects. */
        wall_3d **walls;
        /** A pointer to the next free position to add a wall pointer.
         */
        wall_3d **wep;
        wall_list_3d();
        ~wall_list_3d();
        /** Adds a wall to the list.
         * \param[in] w the wall to add. */
        inline void add_wall(wall_3d *w) {
            if(wep==wel) increase_wall_memory();
            *(wep++)=w;
        }
        /** Adds a wall to the list.
         * \param[in] w a reference to the wall to add. */
        inline void add_wall(wall_3d &w) {add_wall(&w);}
        void add_wall(wall_list_3d &wl);
        /** Determines whether a given position is inside all of the walls on
         * the list.
         * \param[in] (x,y,z) the position to test.
         * \return True if it is inside, false if it is outside. */
        inline bool point_inside_walls(double x,double y,double z) {
            for(wall_3d **wp=walls;wp<wep;wp++) if(!((*wp)->point_inside(x,y,z))) return false;
            return true;
        }
        /** Cuts a Voronoi cell by all of the walls currently on the list.
         * \param[in] c a reference to the Voronoi cell class.
         * \param[in] (x,y,z) the position of the cell.
         * \return True if the cell still exists, false if the cell is deleted.
         */
        template<class c_class>
        bool apply_walls(c_class &c,double x,double y,double z) {
            for(wall_3d **wp=walls;wp<wep;wp++) if(!((*wp)->cut_cell(c,x,y,z))) return false;
            return true;
        }
        void deallocate();
    protected:
        void increase_wall_memory();
        /** A pointer to the limit of the walls array, used to determine when
         * array is full. */
        wall_3d **wel;
        /** The current amount of memory allocated for walls. */
        int current_wall_size;
};

}

#endif
