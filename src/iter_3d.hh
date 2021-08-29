// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file iter_3d.hh
 * \brief Header file for the 3D container iterators and related classes. */

#ifndef VOROPP_ITER_3D_HH
#define VOROPP_ITER_3D_HH

#include "particle_order.hh"
#include "container_3d.hh"
#include "c_info.hh"

namespace voro {

// XXX CHR - for small functions (i.e. those that take one or two lines) I
// moved them to be defined inline here. That can improve efficiency since the
// compiler doesn't have to make a separate function call.
//
// Note that constructors defined in the class can be inlined (even without the
// inline keyword).
//
// For functions that are declared here but not defined, I removed the
// comments. The comments appear in the .cc file to avoid repetition.

class container_base_3d::iterator : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        c_info ptr;
        int* co;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        iterator() : co(0) {}
        iterator(int* _co);
        /** Initializes the iterator.
         * \param[in] co_ a pointer to the particle count array.
         * \param[in] ptr_ information on the particle to point to. */
        iterator(int* co_,c_info ptr_) : ptr(ptr_), co(co_) {}
        /** Initializes the iterator as a copy of another.
         * \param[in] ci a reference to an existing iterator. */
        iterator(const iterator& ci) : ptr(ci.ptr), co(ci.co) {}
        /** Sets the iterator to equal another.
         * \param[in] other the iterator to copy. */
        inline iterator& operator=(iterator other) {
            co=other.co;
            ptr=other.ptr;
            return *this;
        }
        /** Evaluates if this iterator is equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they are equal, false otherwise. */
        inline bool operator==(const iterator& rhs) const {
            return ptr.ijk==rhs.ptr.ijk&&ptr.q==rhs.ptr.q;
            // XXX CHR - I guess we don't check the co pointers are equal?
            // That's probably fine, and is likely more efficient. I presume
            // that for these comparisons, you can assume the iterators are
            // declared for the same container.
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        inline bool operator!=(const iterator& rhs) const {
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
        friend class container_base_3d;
};

class subset_info_3d {
    public:
        friend class container_base_3d::iterator_subset;
        friend class container_base_3d;
        template<class c_class>
        subset_info_3d(c_class& con) : nx(con.nx), ny(con.ny), nz(con.nz),
            nxy(con.nxy), nxyz(con.nxyz), ps(con.ps), p(con.p), id(con.id),
            co(con.co), ax(con.ax), ay(con.ay), az(con.az), sx(con.bx-ax),
            sy(con.by-ay), sz(con.bz-az), xsp(con.xsp), ysp(con.ysp),
            zsp(con.zsp), x_prd(con.x_prd), y_prd(con.y_prd), z_prd(con.z_prd) {}
        ~subset_info_3d() {}
        subset_mode mode;
        int nx;
        int ny;
        int nz;
        int nxy;
        int nxyz;
        int ps;
        double **p;
        int **id;
        int *co;
        void setup_sphere(double vx,double vy,double vz,double r,bool bounds_test=true);
        void setup_box(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,bool bounds_test=true);
        /** Sets up the class constants to loop over all particles in a
         * rectangular subgrid of blocks.
         * \param[in] (ai_,bi_) the subgrid range in the x direction, inclusive
         *                      of both ends.
         * \param[in] (aj_,bj_) the subgrid range in the y direction, inclusive
         *                      of both ends.
         * \param[in] (ak_,bk_) the subgrid range in the z direction, inclusive
         *                      of both ends. */
        void setup_intbox(int ai_,int bi_,int aj_,int bj_,int ak_,int bk_) {
            ai=ai_;bi=bi_;aj=aj_;bj=bj_;ak=ak_;bk=bk_;
            mode=no_check;
            setup_common();
        }
        double apx,apy,apz;
    private:
        double ax,ay,az,sx,sy,sz,xsp,ysp,zsp;
        bool x_prd,y_prd,z_prd;
        double v0,v1,v2,v3,v4,v5;
        int ai,bi,aj,bj,ak,bk;
        int di,dj,dk,inc1,inc2;
        double ddi,ddj,ddk,aapx,aapy,aapz; // XXX CHR - shouldn't ddi, ddj, and ddk be integers?
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        void setup_common();
        void previous_block_iter(int &ijk_,int &i_,int &j_,int &k_,int &ci_,int &cj_,int &ck_,double &px_,double &py_,double &pz_);
        bool out_of_bounds(int ijk_,int q_,double px_,double py_,double pz_);
};

class container_base_3d::iterator_subset : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        c_info ptr;
        subset_info_3d* cl_iter;
        int i,j,k,ci,cj,ck;
        double px,py,pz;
        /** Computes whether the current point is out of bounds, relative to the
         * current loop setup.
         * \param[in] ijk_ the current block.
         * \param[in] q_ the index of the particle in the current block.
         * \return True if the point is out of bounds, false otherwise. */
        inline bool out_of_bounds_iter(int ijk_,int q_) const {
            return cl_iter->out_of_bounds(ijk_,q_,px,py,pz);
        }
        void next_block_iter(int &ijk_);
        void previous_block_iter(int &ijk_);
        // XXX CHR - is there any reason to initialize cl_iter here? A blank
        // iterator is never going to be used - you'd have to later copy-assign
        // it, or set up the variables another way. In those cases cl_iter will
        // be initialized then.
        iterator_subset() : cl_iter(0) {}
        iterator_subset(subset_info_3d* si_);
        iterator_subset(subset_info_3d* si_,c_info ptr_,int i_,int j_,int k_);
        /** Initializes the iterator as a copy of another.
         * \param[in] ci_ a reference to an existing iterator. */
        iterator_subset(const iterator_subset& ci_) : ptr(ci_.ptr),
            cl_iter(ci_.cl_iter), i(ci_.i), j(ci_.j), k(ci_.k), ci(ci_.ci),
            cj(ci_.cj), ck(ci_.ck), px(ci_.px), py(ci_.py), pz(ci_.pz) {}
        iterator_subset& operator=(iterator_subset other);
        /** Evaluates if this iterator is equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they are equal, false otherwise. */
        inline bool operator==(const iterator_subset& rhs) const {
            return ptr.ijk==rhs.ptr.ijk&&ptr.q==rhs.ptr.q&&i==rhs.i&&j==rhs.j&&k==rhs.k;
            // XXX CHR - Is it necessary to check ijk? Also, in the >= and <=
            // comparisons later, ijk is not checked.
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        inline bool operator!=(const iterator_subset& rhs) const {
            return ptr.ijk!=rhs.ptr.ijk||ptr.q!=rhs.ptr.q||i!=rhs.i||j!=rhs.j||k!=rhs.k;
        }
        /** Dereferences the iterator as an rvalue. */
        c_info& operator*() {return ptr;}
        /** Dereferences the iterator as an rvalue. */
        c_info* operator->() {return &ptr;}
        iterator_subset& operator++();
        iterator_subset operator++(int);
        iterator_subset& operator--();
        iterator_subset operator--(int);
        difference_type operator-(const iterator_subset& rhs) const;
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
            return k>rhs.k||(k==rhs.k&&(j>rhs.j||(j==rhs.j&&(i>rhs.i||(i==rhs.i&&ptr.q>rhs.ptr.q)))));
        }
        /** Evaluates if this iterator is less than another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less, false otherwise. */
        inline bool operator<(const iterator_subset& rhs) const {
            return k<rhs.k||(k==rhs.k&&(j<rhs.j||(j==rhs.j&&(i<rhs.i||(i==rhs.i&&ptr.q<rhs.ptr.q)))));
        }
        /** Evaluates if this iterator is greater than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is greater or equal, false otherwise. */
        inline bool operator>=(const iterator_subset& rhs) const {
            return k>rhs.k||(k==rhs.k&&(j>rhs.j||(j==rhs.j&&(i>rhs.i||(i==rhs.i&&ptr.q>=rhs.ptr.q)))));
        }
        /** Evaluates if this iterator is less than or equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if it is less or equal, false otherwise. */
        inline bool operator<=(const iterator_subset& rhs) const {
            return k<rhs.k||(k==rhs.k&&(j<rhs.j||(j==rhs.j&&(i<rhs.i||(i==rhs.i&&ptr.q<=rhs.ptr.q)))));
        }
        iterator_subset& operator+=(const difference_type& incre);
        iterator_subset& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
        friend class container_base_3d;
};

class container_base_3d::iterator_order : public std::iterator<std::random_access_iterator_tag,c_info,int> {
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
        // XXX CHR - Do we need to pass in both a c_info and ptr_n? If we know
        // ptr_n, then we can set c_info. I guess it depends on the situations
        // where this function is called?
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
        inline bool operator==(const iterator_order& rhs) const {
            return ptr_n==rhs.ptr_n;
        }
        /** Evaluates if this iterator is not equal to another.
         * \param[in] rhs a reference to another iterator.
         * \return True if they aren't equal, false if they are. */
        inline bool operator!=(const iterator_order& rhs) const {
            return ptr_n!=rhs.ptr_n;
        }
        /** Dereferences the iterator as an rvalue. */
        inline c_info& operator*(){return ptr;}
        /** Dereferences the iterator as an rvalue. */
        inline c_info* operator->(){return &ptr;}
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
        friend class container_base_3d;
};

}

#endif
