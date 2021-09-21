// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file iter_2d.hh
 * \brief Header file for the 2D container iterators and related classes. */

#ifndef VOROPP_ITER_2D_HH
#define VOROPP_ITER_2D_HH

#include "particle_order.hh"
#include "container_2d.hh"
#include "c_info.hh"

namespace voro {

class container_base_2d::iterator : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        c_info ptr;
        int* co;
        int nxy;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        iterator(){}
        iterator(int* co_,int _nxy);
        /** Initializes the iterator.
         * \param[in] co_ a pointer to the particle count array.
         * \param[in] ptr_ information on the particle to point to. */
        iterator(int* co_,c_info ptr_,int _nxy) : ptr(ptr_), co(co_), nxy(_nxy) {}
        /** Initializes the iterator as a copy of another.
         * \param[in] ci a reference to an existing iterator. */
        iterator(const iterator& ci) : ptr(ci.ptr), co(ci.co), nxy(ci.nxy) {}
        /** Sets the iterator to equal another.
         * \param[in] other the iterator to copy. */
        iterator& operator=(iterator other){
            co=other.co;
            ptr=other.ptr;
            nxy=other.nxy;
            return *this;
        }
        /** Evaluates if this iterator is equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they are equal, false otherwise. */
        inline bool operator==(const iterator& rhs) const {
            return ptr.ijk==rhs.ptr.ijk&&ptr.q==rhs.ptr.q;
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        bool operator!=(const iterator& rhs) const {
            return ptr.ijk!=rhs.ptr.ijk||ptr.q!=rhs.ptr.q;
        }
        /** Dereferences the iterator as an rvalue. */
        inline c_info& operator*() {return ptr;}
        /** Dereferences the iterator as an rvalue. */
        inline c_info* operator->() {return &ptr;}
        iterator& operator++();
        iterator operator++(int);
        iterator& operator--();
        iterator operator--(int);
        difference_type operator-(const iterator& rhs) const;
        /** Calculates a new iterator by adding elements.
         * \param[in] incre the number of elements to increment by. */
        inline iterator operator+(const difference_type& incre) const {
            iterator tmp(*this);
            return tmp+=incre;
        }
        /** Calculates a new iterator by subtracting elements.
         * \param[in] decre the number of elements to decrement by. */
        inline iterator operator-(const difference_type& decre) const {
            iterator tmp(*this);
            return tmp-=decre;
        }
        /** Evaluates if this iterator is greater than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater, false otherwise. */
        inline bool operator>(const iterator& rhs) const {
            return ptr.ijk>rhs.ptr.ijk||(ptr.ijk==rhs.ptr.ijk&&ptr.q>rhs.ptr.q);
        }
        /** Evaluates if this iterator is less than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less, false otherwise. */
        inline bool operator<(const iterator& rhs) const {
            return ptr.ijk<rhs.ptr.ijk||(ptr.ijk==rhs.ptr.ijk&&ptr.q<rhs.ptr.q);
        }
        /** Evaluates if this iterator is greater than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater or equal, false otherwise. */
        inline bool operator>=(const iterator& rhs) const {
            return ptr.ijk>rhs.ptr.ijk||(ptr.ijk==rhs.ptr.ijk&&ptr.q>=rhs.ptr.q);
        }
        /** Evaluates if this iterator is less than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less or equal, false otherwise. */
        inline bool operator<=(const iterator& rhs) const {
            return ptr.ijk<rhs.ptr.ijk||(ptr.ijk==rhs.ptr.ijk&&ptr.q<=rhs.ptr.q);
        }
        iterator& operator+=(const difference_type& incre);
        iterator& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
        friend class container_base_2d;
        friend void swap(iterator& a, iterator& b){
            std::swap(a.ptr.ijk, b.ptr.ijk);std::swap(a.ptr.q, b.ptr.q);
        }
};

class subset_info_2d {
    public:
        friend class container_base_2d::iterator_subset;
        friend class container_base_2d;
        template<class c_class_2d>
        subset_info_2d(c_class_2d& con): nx(con.nx), ny(con.ny), nxy(con.nxy),
            ps(con.ps), p(con.p), id(con.id), co(con.co), ax(con.ax),
            ay(con.ay), sx(con.bx-ax), sy(con.by-ay), xsp(con.xsp),
            ysp(con.ysp), x_prd(con.x_prd), y_prd(con.y_prd) {}
        ~subset_info_2d() {}
        subset_mode mode;
        int nx;
        int ny;
        int nxy;
        int ps;
        double **p;
        int **id;
        int *co;
        double apx,apy;
        void setup_circle(double vx,double vy,double r,bool bounds_test=true);
        void setup_box(double xmin,double xmax,double ymin,double ymax,bool bounds_test=true);
        /** Sets up the class constants to loop over all particles in a
         * rectangular subgrid of blocks.
         * \param[in] (ai_,bi_) the subgrid range in the x direction, inclusive
         *                      of both ends.
         * \param[in] (aj_,bj_) the subgrid range in the y direction, inclusive
         *                      of both ends. */
        inline void setup_intbox(int ai_,int bi_,int aj_,int bj_) {
            ai=ai_;bi=bi_;aj=aj_;bj=bj_;
            mode=no_check;
            setup_common();
        }
    private:
        double ax,ay,sx,sy,xsp,ysp;
        bool x_prd,y_prd;
        double v0,v1,v2,v3;
        int ai,bi,aj,bj;
        int di,dj,inc1;
        int ddi,ddj;
        double aapx,aapy; // XXX CHR - shouldn't ddi and ddj be integers?
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        void setup_common();
};

class container_base_2d::iterator_subset : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        c_info ptr;
        subset_info_2d* cl_iter;
        int i,j,ci,cj;
        double px,py;
        /** Computes whether the current point is out of bounds, relative to
         * the current loop setup.
         * \return True if the point is out of bounds, false otherwise. */
        bool out_of_bounds();
        bool next_block();
        bool previous_block();
        iterator_subset(){}
        iterator_subset(subset_info_2d* si_);
        iterator_subset(subset_info_2d* si_,c_info ptr_,int i_,int j_);
        /** Initializes the iterator to be a copy of another.
         * \param[in] ci_ a reference to an existing iterator. */
        iterator_subset(const iterator_subset& ci_) : ptr(ci_.ptr),
            cl_iter(ci_.cl_iter), i(ci_.i), j(ci_.j), ci(ci_.ci), cj(ci_.cj),
            px(ci_.px), py(ci_.py) {}
        iterator_subset& operator=(iterator_subset other);
         /** Evaluates if this iterator is equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they are equal, false otherwise. */
        inline bool operator==(const iterator_subset& rhs) const {
            return ptr.q==rhs.ptr.q&&i==rhs.i&&j==rhs.j;
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        inline bool operator!=(const iterator_subset& rhs) const {
            return ptr.q!=rhs.ptr.q||i!=rhs.i||j!=rhs.j;
        }
        /** Dereferences the iterator as an rvalue. */
        inline c_info& operator*() {return ptr;}
        /** Dereferences the iterator as an rvalue. */
        inline c_info* operator->() {return &ptr;}
        iterator_subset& operator++();
        iterator_subset operator++(int);
        iterator_subset& operator--();
        iterator_subset operator--(int);
        difference_type operator-(const iterator_subset& rhs) const;
        iterator_subset& operator+=(const difference_type& incre);
        iterator_subset& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
        /** Calculates a new iterator by adding elements.
         * \param[in] incre the number of elements to increment by. */
        inline iterator_subset operator+(const difference_type& incre) const {
            iterator_subset tmp(*this);
            return tmp+=incre;
        }
        /** Calculates a new iterator by subtracting elements.
         * \param[in] decre the number of elements to decrement by. */
        inline iterator_subset operator-(const difference_type& decre) const {
            iterator_subset tmp(*this);
            return tmp-=decre;
        }
        /** Evaluates if this iterator is greater than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater, false otherwise. */
        inline bool operator>(const iterator_subset& rhs) const {
            return j>rhs.j||(j==rhs.j&&(i>rhs.i||(i==rhs.i&&ptr.q>rhs.ptr.q)));
        }
        /** Evaluates if this iterator is less than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less, false otherwise. */
        inline bool operator<(const iterator_subset& rhs) const {
            return j<rhs.j||(j==rhs.j&&(i<rhs.i||(i==rhs.i&&ptr.q<rhs.ptr.q)));
        }
        /** Evaluates if this iterator is greater than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater or equal, false otherwise. */
        inline bool operator>=(const iterator_subset& rhs) const {
            return j>rhs.j||(j==rhs.j&&(i>rhs.i||(i==rhs.i&&ptr.q>=rhs.ptr.q)));
        }
        /** Evaluates if this iterator is less than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less or equal, false otherwise. */
        inline bool operator<=(const iterator_subset& rhs) const {
            return j<rhs.j||(j==rhs.j&&(i<rhs.i||(i==rhs.i&&ptr.q<=rhs.ptr.q)));
        }
        friend class container_base_2d;
        friend void swap(iterator_subset& a, iterator_subset& b){
            std::swap(a.ptr.ijk, b.ptr.ijk);std::swap(a.ptr.q, b.ptr.q);
            std::swap(a.i,b.i);std::swap(a.j,b.j);
            std::swap(a.ci,b.ci);std::swap(a.cj,b.cj);
            std::swap(a.px,b.px);std::swap(a.py,b.py);
        }
};

class container_base_2d::iterator_order : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        c_info ptr;
        int* cp_iter;
        int* op_iter;
        int ptr_n;
        int pn_upper_bound; //(op_iter-cp_iter)/2, ie. number of particles; ptr_n< than this number to be in range
        int nxy;
        iterator_order(){}
        /** Initializes the iterator.
         * \param[in] vo_ a reference to the particle_order class to follow. */
        iterator_order(particle_order& vo_, int nxy_) : cp_iter(vo_.o), op_iter(vo_.op), ptr_n(0), pn_upper_bound((op_iter-cp_iter)/2),
        nxy(nxy_)
        {
            if(pn_upper_bound!=0){ptr.set(cp_iter[0],cp_iter[1]);} else{ptr.set(nxy,0);} //if empty particle_order, set to one over the end 
        }
        /** Initializes the iterator, and sets it to point at a given particle.
         * \param[in] vo_ a reference to the particle_order class to follow.
         * \param[in] ptr_ the particle to point to.
         * \param[in] ptr_n_ the position in the particle_order class of the
         *                   particle. */
        iterator_order(particle_order& vo_,int ptr_n_,int _nxy)
            : cp_iter(vo_.o), op_iter(vo_.op), ptr_n(ptr_n_),pn_upper_bound((op_iter-cp_iter)/2),nxy(_nxy)
        {
            if(ptr_n>=0 && ptr_n<pn_upper_bound){ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);}
            else if(ptr_n<0){ptr_n=-1;ptr.set(0,-1);} //out of range, return 1-before-the-start
            else{ptr_n=pn_upper_bound; ptr.set(nxy,0);} //out of range, return 1-over-the-last
        }
        /** Initializes the iterator as a copy of another.
         * \param[in] ci a reference to an existing iterator. */
        iterator_order(const iterator_order& ci) : ptr(ci.ptr),
            cp_iter(ci.cp_iter), op_iter(ci.op_iter), ptr_n(ci.ptr_n), pn_upper_bound(ci.pn_upper_bound), nxy(ci.nxy)
        {}
        /** Initializes the iterator as a copy of another.
         * \param[in] other an existing iterator. */
        iterator_order& operator=(iterator_order other) {
            cp_iter=other.cp_iter;
            op_iter=other.op_iter;
            ptr_n=other.ptr_n;
            ptr=other.ptr;
            pn_upper_bound=other.pn_upper_bound;
            nxy=other.nxy;
            return *this;
        }
        /** Evaluates if this iterator is equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they are equal, false otherwise. */
        bool operator==(const iterator_order& rhs) const {
            return ptr_n==rhs.ptr_n;
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        bool operator!=(const iterator_order& rhs) const {
            return ptr_n!=rhs.ptr_n;
        }
        /** Dereferences the iterator as an rvalue. */
        inline c_info& operator*() {return ptr;}
        /** Dereferences the iterator as an rvalue. */
        inline c_info* operator->() {return &ptr;}
        iterator_order& operator++();
        iterator_order operator++(int);
        iterator_order& operator--();
        iterator_order operator--(int);
        difference_type operator-(const iterator_order& rhs) const;
        /** Calculates a new iterator by adding elements.
         * \param[in] incre the number of elements to increment by. */
        inline iterator_order operator+(const difference_type& incre) const {
            iterator_order tmp(*this);
            return tmp+=incre;
        }
        /** Calculates a new iterator by subtracting elements.
         * \param[in] decre the number of elements to decrement by. */
        inline iterator_order operator-(const difference_type& decre) const {
            iterator_order tmp(*this);
            return tmp-=decre;
        }
        /** Evaluates if this iterator is greater than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater, false otherwise. */
        inline bool operator>(const iterator_order& rhs) const {
            return ptr_n>rhs.ptr_n;
        }
        /** Evaluates if this iterator is less than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less, false otherwise. */
        inline bool operator<(const iterator_order& rhs) const {
            return ptr_n<rhs.ptr_n;
        }
        /** Evaluates if this iterator is greater than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater or equal, false otherwise. */
        inline bool operator>=(const iterator_order& rhs) const {
            return ptr_n>=rhs.ptr_n;
        }
        /** Evaluates if this iterator is less than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less or equal, false otherwise. */
        inline bool operator<=(const iterator_order& rhs) const {
            return ptr_n<=rhs.ptr_n;
        }
        iterator_order& operator+=(const difference_type& incre);
        iterator_order& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
        friend class container_base_2d;
        friend void swap(iterator_order& a, iterator_order& b){
            std::swap(a.ptr.ijk, b.ptr.ijk);std::swap(a.ptr.q, b.ptr.q);std::swap(a.ptr_n,b.ptr_n);
        }
};

}

#endif
