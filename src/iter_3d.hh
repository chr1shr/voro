// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file iter_3d.hh
 * \brief Header file for the 3D container iterators and related classes. */

#ifndef VOROPP_ITER_3D_HH
#define VOROPP_ITER_3D_HH

#include "container_3d.hh"

namespace voro {

struct c_info {
    int ijk;
    int q;
    c_info() {}
    // XXX CHR - I've added a constructor here, which should simplify some of
    // the code you wrote. Rather than writing "c_info ci;ci.set(ijk,q);return
    // ci;" you can write "return c_info(ijk,q);".
    c_info(int ijk_,int q_) ijk(ijk_), q(q_) {}
    inline void set(int ijk_,int q_){
        ijk=ijk_; q=q_;
    }
};

class container_base_3d::iterator : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        c_info ptr;
        int* co_iter;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        iterator() : co_iter(0){}
        iterator(int* _co);
        iterator(int* _co,c_info _ptr);
        // XXX CHR - note that I removed comments on functions here, to avoid redundancy with
        // the comments in the .cc file
        iterator(const iterator& ci);
        iterator& operator=(iterator other);
        bool operator==(const iterator& rhs) const;
        bool operator!=(const iterator& rhs) const;
        c_info& operator*();
        c_info* operator->();
        iterator& operator++();
        iterator operator++(int);
        iterator& operator--();
        iterator operator--(int);
        difference_type operator-(const iterator& rhs) const;
        iterator operator+(const difference_type& incre) const;
        iterator operator-(const difference_type& decre) const;
        bool operator>(const iterator& rhs) const;
        bool operator<(const iterator& rhs) const;
        bool operator>=(const iterator& rhs) const;
        bool operator<=(const iterator& rhs) const;
        iterator& operator+=(const difference_type& incre);
        iterator& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
        friend class container_base_3d;
};

class subset_info {
    public:
        friend class container_base_3d::iterator_subset;
        friend class container_base_3d;
        template<class c_class>
        subset_info(c_class& con) : nx(con.nx), ny(con.ny), nz(con.nz),
            nxy(con.nxy), nxyz(con.nxyz), ps(con.ps), p(con.p), id(con.id),
            co(con.co), ax(con.ax), ay(con.ay), az(con.az), sx(con.bx-ax),
            sy(con.by-ay), sz(con.bz-az), xsp(con.xsp), ysp(con.ysp),
            zsp(con.zsp), x_prd(con.x_prd), y_prd(con.y_prd), z_prd(con.z_prd)
        ~subset_info(){}
        c_loop_subset_mode mode;
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
        void setup_intbox(int ai_,int bi_,int aj_,int bj_,int ak_,int bk_);
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
        bool out_of_bounds(int ijk_,int q_,double px_,double py_,double pz_);
};

class container_base_3d::iterator_subset : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        c_info ptr;
        subset_info* cl_iter;
        double px,py,pz;
        int ci,cj,ck,i,j,k;
        bool out_of_bounds_iter(int ijk_,int q_) const;
        void next_block_iter(int &ijk_);
        void previous_block_iter(int &ijk_);
        friend class container_base_3d;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        iterator_subset() : cl_iter(0) {} // XXX CHR - is there any reason to initialize cl_iter here? A blank iterator is never going to be used - you'd have to later copy-assign it, or set up the variables another way. In those cases cl_iter will be initialized then.
        iterator_subset(subset_info* _si);
        iterator_subset(subset_info* _si,c_info _ptr,int _i,int _j,int _k);
        iterator_subset(const iterator_subset& _ci) : ptr(_ci.ptr),
            cl_iter(_ci.cl_iter), i(_ci.i), j(_ci.j), k(_ci.k), ci(_ci.ci),
            cj(_ci.cj), ck(_ci.ck), px(_ci.px), py(_ci.py), pz(_ci.pz) {}
        iterator_subset& operator=(iterator_subset other);
        bool operator==(const iterator_subset& rhs) const;
        bool operator!=(const iterator_subset& rhs) const;
        c_info& operator*() {return ptr;}
        c_info* operator->() {return &ptr;}
        iterator_subset& operator++();
        iterator_subset operator++(int);
        iterator_subset& operator--();
        iterator_subset operator--(int);
        difference_type operator-(const iterator_subset& rhs) const;
        iterator_subset operator+(const difference_type& incre) const;
        iterator_subset operator-(const difference_type& decre) const;
        bool operator>(const iterator_subset& rhs) const;
        bool operator<(const iterator_subset& rhs) const;
        bool operator>=(const iterator_subset& rhs) const;
        bool operator<=(const iterator_subset& rhs) const;
        iterator_subset& operator+=(const difference_type& incre);
        iterator_subset& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
};

class container_base_3d::iterator_order : public std::iterator<std::random_access_iterator_tag,c_info,int> {
    public:
        c_info ptr;
        int* cp_iter;
        int ptr_n;
        friend class container_base_3d;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;
        iterator_order() : cp_iter(0) {}
        iterator_order(particle_order& vo_);
        iterator_order(particle_order& vo_,c_info _ptr,int ptr_n_)
            : ptr(_ptr), cp_iter(vo_.o), ptr_n(ptr_n_) {}
        iterator_order(const iterator_order& ci) : ptr(ci.ptr),
            cp_iter(ci.cp_iter), ptr_n(ci.ptr_n) {}
        iterator_order& operator=(iterator_order other);
        bool operator==(const iterator_order& rhs) const;
        bool operator!=(const iterator_order& rhs) const;
        c_info& operator*(){return ptr;}
        c_info* operator->(){return &ptr;}
        iterator_order& operator++();
        iterator_order operator++(int);
        iterator_order& operator--();
        iterator_order operator--(int);
        difference_type operator-(const iterator_order& rhs) const;
        iterator_order operator+(const difference_type& incre) const;
        iterator_order operator-(const difference_type& decre) const;
        bool operator>(const iterator_order& rhs) const;
        bool operator<(const iterator_order& rhs) const;
        bool operator>=(const iterator_order& rhs) const;
        bool operator<=(const iterator_order& rhs) const;
        iterator_order& operator+=(const difference_type& incre);
        iterator_order& operator-=(const difference_type& decre);
        c_info& operator[](const difference_type& incre) const;
};

}

#endif
