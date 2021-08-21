// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file iter_2d.hh
 * \brief Header file for the 2D container iterators and related classes. */

#ifndef VOROPP_ITER_2D_HH
#define VOROPP_ITER_2D_HH

#include "container_2d.hh"
#include "c_info.hh"

class container_base_2d::iterator : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        c_info ptr;
        int* co_iter;
        iterator() : co_iter(0) {}
        iterator(int* _co);
        /** Initializes the iterator.
         * \param[in] co_ a pointer to the particle count array.
         * \param[in] ptr_ information on the particle to point to. */
        iterator(int* co_,c_info ptr_)
            : ptr(_ptr), co_iter(_co) {}
        /** Initializes the iterator as a copy of another.
         * \param[in] ci a reference to an existing iterator. */
        iterator(const iterator& ci)
            : ptr(ci.ptr), co(ci.co) {}
        iterator& operator=(iterator other);
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
};

class subset_info_2d {
    public:
        /** The types of geometrical region to loop over. */
        enum subset_mode {
            circle,
            box,
            no_check
        };
        template<class c_class_2d>
        subset_info(c_class_2d& con): nx(con.nx), ny(con.ny), nxy(con.nxy),
            ps(con.ps), p(con.p), id(con.id), co(con.co), ax(con.ax),
            ay(con.ay), sx(con.bx-ax), sy(con.by-ay), xsp(con.xsp),
            ysp(con.ysp), x_prd(con.x_prd), y_prd(con.y_prd) {}
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
        void setup_intbox(int ai_,int bi_,int aj_,int bj_) {
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
        double ddi,ddj,aapx,aapy; // XXX CHR - shouldn't ddi and ddj be integers?
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        void setup_common();
        //to be used in con.begin(cli) and con.end(cli)
        void previous_block_iter(int &ijk_,int &i_,int &j_,int &ci_,int &cj_,double &px_,double &py_);
        bool out_of_bounds(int ijk_,int q_,double px_,double py_);
        friend class container_base_2d;
        friend class container_base_2d::iterator_subset;
};

class container_base_2d::iterator_subset : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        c_info ptr;
        subset_info* cl_iter;
        double px,py;
        int ci,cj,i,j;
        /** Computes whether the current point is out of bounds, relative to
         * the current loop setup.
         * \param[in] ijk_ the current block.
         * \param[in] q_ the index of the particle in the current block.
         * \return True if the point is out of bounds, false otherwise. */
        inline bool out_of_bounds_iter(int ijk_,int q_) const {
            return cl_iter->out_of_bounds(ijk_,q_,px,py);
        }
        void next_block_iter(int &ijk_);
        void previous_block_iter(int &ijk_);
        iterator_subset() : cl_iter(0) {};
        iterator_subset(subset_info* si_);
        iterator_subset(subset_info* si_,c_info ptr_,int i_,int j_);
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
            return ptr.ijk==rhs.ptr.ijk&&ptr.q==rhs.ptr.q&&i==rhs.i&&j==rhs.j;
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        inline bool operator!=(const iterator_subset& rhs) const {
            return ptr.ijk!=rhs.ptr.ijk||ptr.q!=rhs.ptr.q||i!=rhs.i||j!=rhs.j;
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
        /** Calculates a new iterator by adding elements.
         * \param[in] incre the number of elements to increment by. */
        inline iterator_subset operator+(const difference_type& incre) const {
            iterator tmp(*this);
            return tmp+=incre;
        }
        /** Calculates a new iterator by subtracting elements.
         * \param[in] decre the number of elements to decrement by. */
        inline iterator_subset operator-(const difference_type& decre) const {
            iterator tmp(*this);
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
        iterator_subset& operator+=(const difference_type& incre);
        iterator_subset& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
        friend class container_base_2d;
};

class container_base_2d::iterator_order : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        c_info ptr;
        int* cp_iter;
        int ptr_n;
        iterator_order() : cp_iter(0) {}
        /** Initializes the iterator.
         * \param[in] vo_ a reference to the particle_order class to follow. */
        iterator_order(particle_order& vo_) : cp_iter(vo_.o), ptr_n(0) {
            ptr.set(cp_iter[0],cp_iter[1]);
        }
        /** Initializes the iterator, and sets it to point at a given particle.
         * \param[in] vo_ a reference to the particle_order class to follow.
         * \param[in] ptr_ the particle to point to.
         * \param[in] ptr_n_ the position in the particle_order class of the
         *                   particle. */
        iterator_order(particle_order& vo_,c_info ptr_,int ptr_n_)
            : ptr(ptr_), cp_iter(vo_.o), ptr_n(ptr_n_) {}
        /** Initializes the iterator as a copy of another.
         * \param[in] ci a reference to an existing iterator. */
        iterator_order(const iterator_order& ci) : ptr(ci.ptr),
            cp_iter(ci.cp_iter), ptr_n(ci.ptr_n) {}
        /** Initializes the iterator as a copy of another.
         * \param[in] other an existing iterator. */
        iterator_order& operator=(iterator_order other) {
            cp_iter=other.cp_iter;
            ptr_n=other.ptr_n;
            ptr=other.ptr;
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
};

#endif
