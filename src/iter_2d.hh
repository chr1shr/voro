// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#ifndef VOROPP_ITER_2D_HH
#define VOROPP_ITER_2D_HH

#include "container_2d.hh"

struct c_info {
    int ijk;
    int q;
    inline void set(int ijk_,int q_){
        ijk=ijk_; q=q_;
    }
};

//forward declaration,so it can be used
//inside container_base_2d for iterator_subset,before complete definition.
class subset_info;

class container_base_2d::iterator : public std::iterator<std::random_access_iterator_tag,c_info,int>
{
    public:
        c_info ptr;
        int* co_iter;

        friend class container_base_2d;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;

        //default-constructible??
        iterator();
        //constructor: point to first particle in container
        iterator(int* _co);
        //constructor: point to a particle in container
        iterator(int* _co,c_info _ptr);
        //copy-constructible
        iterator(const iterator& ci);
        //copy-assignable
        iterator& operator=(iterator other);
        //destructible??
        ~iterator();

        //Can be compared for equivalence using the equality/inequality operators
        bool operator==(const iterator& rhs) const;
        bool operator!=(const iterator& rhs) const;

        //Can be dereferenced as an rvalue (if in a dereferenceable state).
        c_info& operator*();
        c_info* operator->();

        //Can be incremented (if in a dereferenceable state).
        //The result is either also dereferenceable or a past-the-end iterator.
        //Two iterators that compare equal,keep comparing equal after being both increased.
        iterator& operator++();
        iterator operator++(int);

        //Can be decremented (if a dereferenceable iterator value precedes it).
        iterator& operator--();
        iterator operator--(int);

        //Supports the arithmetic operators + and - between an iterator and an integer value,or subtracting an iterator from another.
        difference_type operator-(const iterator& rhs) const;
        iterator operator+(const difference_type& incre) const;
        iterator operator-(const difference_type& decre) const;

        //Can be compared with inequality relational operators (<,>,<= and >=).
        bool operator>(const iterator& rhs) const;
        bool operator<(const iterator& rhs) const;
        bool operator>=(const iterator& rhs) const;
        bool operator<=(const iterator& rhs) const;

        //Supports compound assignment operations += and -=
        iterator& operator+=(const difference_type& incre);
        iterator& operator-=(const difference_type& decre);

        //Supports the offset dereference operator ([])
        c_info& operator[](const difference_type& incre) const;

        //swappable??

};

class subset_info {
    public:
        friend class container_base_2d::iterator_subset;
        friend class container_base_2d;
        template<class c_class_2d>
        subset_info(c_class_2d& con): nx(con.nx),ny(con.ny),nxy(con.nxy),
                                ps(con.ps),p(con.p),id(con.id),
                                co(con.co),ax(con.ax),ay(con.ay),
                                sx(con.bx-ax),sy(con.by-ay),xsp(con.xsp),ysp(con.ysp),
                                xperiodic(con.xperiodic),yperiodic(con.yperiodic) {}
        ~subset_info(){}

        c_loop_subset_mode_2d mode;
        int nx;
        int ny;
        int nxy;
        int ps;
        double **p;
        int **id;
        int *co;
        void setup_circle(double vx,double vy,double r,bool bounds_test=true);
        void setup_box(double xmin,double xmax,double ymin,double ymax,bool bounds_test=true);
        void setup_intbox(int ai_,int bi_,int aj_,int bj_);
        double apx,apy;
    private:
        double ax,ay,sx,sy,xsp,ysp;
        bool xperiodic,yperiodic;
        double v0,v1,v2,v3;
        int ai,bi,aj,bj;
        int di,dj,inc1;
        double ddi,ddj,aapx,aapy;
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        void setup_common();
        //to be used in con.begin(cli) and con.end(cli)
        void previous_block_iter2(int &ijk_,int &i_,int &j_,int &ci_,int &cj_,double &px_,double &py_);
        bool out_of_bounds_2(int ijk_,int q_,double px_,double py_);

};

class container_base_2d::iterator_subset : public std::iterator<std::random_access_iterator_tag,c_info,int>
{
    public:
        c_info ptr;
        subset_info* cl_iter;
        double px,py;
        int ci,cj,i,j;

        bool out_of_bounds_iter(int ijk_,int q_) const;
        void next_block_iter(int &ijk_);
        void previous_block_iter(int &ijk_);

        friend class container_base_2d;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;

        //default-constructible
        iterator_subset() : cl_iter(0) {};
        //constructor: point to first particle in container based on subset_info
        iterator_subset(subset_info* _si);

        //constructor: point to a specific particle in container
        iterator_subset(subset_info* _si,c_info _ptr,int _i,int _j);

        //copy-constructible
        iterator_subset(const iterator_subset& _ci);
        //copy-assignable
        iterator_subset& operator=(iterator_subset other);
        //destructible
        ~iterator_subset(){}

        //Can be compared for equivalence using the equality/inequality operators
        bool operator==(const iterator_subset& rhs) const;
        bool operator!=(const iterator_subset& rhs) const;

        //Can be dereferenced as an rvalue (if in a dereferenceable state).
        c_info& operator*();
        c_info* operator->();

        //Can be incremented (if in a dereferenceable state).
        //The result is either also dereferenceable or a past-the-end iterator.
        //Two iterators that compare equal,keep comparing equal after being both increased.
        iterator_subset& operator++();
        iterator_subset operator++(int);
        //Can be decremented (if a dereferenceable iterator value precedes it).
        iterator_subset& operator--();
        iterator_subset operator--(int);
        //Supports the arithmetic operators + and - between an iterator and an integer value,or subtracting an iterator from another.
        difference_type operator-(const iterator_subset& rhs) const;
        iterator_subset operator+(const difference_type& incre) const;
        iterator_subset operator-(const difference_type& decre) const;

        //Can be compared with inequality relational operators (<,>,<= and >=).
        bool operator>(const iterator_subset& rhs) const;
        bool operator<(const iterator_subset& rhs) const;
        bool operator>=(const iterator_subset& rhs) const;
        bool operator<=(const iterator_subset& rhs) const;

        //Supports compound assignment operations += and -=
        iterator_subset& operator+=(const difference_type& incre);
        iterator_subset& operator-=(const difference_type& decre);
        //Supports the offset dereference operator ([])
        c_info& operator[](const difference_type& incre) const;
        //swappable??
};

class container_base_2d::iterator_order : public std::iterator<std::random_access_iterator_tag,c_info,int>
{
    public:
        c_info ptr;
        int* cp_iter;
        int ptr_n;

        friend class container_base_2d;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::pointer pointer;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::reference reference;
        typedef typename std::iterator<std::random_access_iterator_tag,c_info,int>::difference_type difference_type;

        //default-constructible??
        iterator_order() : cp_iter(0){}
        //constructor: point to first particle in container
        iterator_order(particle_order& vo_);
        //constructor: point to a particle in container
        iterator_order(particle_order& vo_,c_info _ptr,int ptr_n_) : ptr(_ptr),cp_iter(vo_.o),ptr_n(ptr_n_) {}
        //copy-constructible
        iterator_order(const iterator_order& ci) : ptr(ci.ptr),cp_iter(ci.cp_iter),ptr_n(ci.ptr_n) {}
        //copy-assignable
        iterator_order& operator=(iterator_order other);
        //destructible??
        ~iterator_order(){}

        //Can be compared for equivalence using the equality/inequality operators
        bool operator==(const iterator_order& rhs) const;
        bool operator!=(const iterator_order& rhs) const;

        //Can be dereferenced as an rvalue (if in a dereferenceable state).
        c_info& operator*();
        c_info* operator->();

        //Can be incremented (if in a dereferenceable state).
        //The result is either also dereferenceable or a past-the-end iterator.
        //Two iterators that compare equal,keep comparing equal after being both increased.
        iterator_order& operator++();
        iterator_order operator++(int);

        //Can be decremented (if a dereferenceable iterator value precedes it).
        iterator_order& operator--();
        iterator_order operator--(int);

        //Supports the arithmetic operators + and - between an iterator and an integer value,or subtracting an iterator from another.
        difference_type operator-(const iterator_order& rhs) const;
        iterator_order operator+(const difference_type& incre) const;
        iterator_order operator-(const difference_type& decre) const;
        //Can be compared with inequality relational operators (<,>,<= and >=).
        bool operator>(const iterator_order& rhs) const;
        bool operator<(const iterator_order& rhs) const;
        bool operator>=(const iterator_order& rhs) const;
        bool operator<=(const iterator_order& rhs) const;
        //Supports compound assignment operations += and -=
        iterator_order& operator+=(const difference_type& incre);
        iterator_order& operator-=(const difference_type& decre);
        //Supports the offset dereference operator ([])
        c_info& operator[](const difference_type& incre) const;
        //swappable??

};

#endif
