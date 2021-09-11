// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#include "iter_2d.hh"

namespace voro {

/** Initializes the iterator, setting it to point at the first particle in the
 * container.
 * \param[in] co_ a pointer to the particle count array. */
container_base_2d::iterator::iterator(int* co_) : co(co_) {
    int ij=0;
    while(co[ij]==0 && ij<nxy) ij++;
    ptr.set(ij,0);
}

/** Sets the iterator to equal another.
 * \param[in] other the iterator to copy. */
container_base_2d::iterator& container_base_2d::iterator::operator=(iterator other) {
    co=other.co;
    ptr=other.ptr;
    return *this;
}

/** Increments the iterator by one element. */
container_base_2d::iterator& container_base_2d::iterator::operator++() {
    int q_=ptr.q,ijk_=ptr.ijk,n=1,diff=q_+n-co[ijk_];
    while(diff>=0) {
        n=n-co[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co[ijk_];
    }
    ptr.set(ijk_,q_+n);
    return *this;
}

/** Increments the iterator by one element. */
container_base_2d::iterator container_base_2d::iterator::operator++(int) {
    iterator tmp(*this);
    int q_=ptr.q,ijk_=ptr.ijk,n=1,diff=q_+n-co[ijk_];
    while(diff>=0) {
        n=n-co[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co[ijk_];
    }
    ptr.set(ijk_,q_+n);
    return tmp;
}

/** Decrements the iterator by one element. */
container_base_2d::iterator& container_base_2d::iterator::operator--() {
    int q_=ptr.q,ijk_=ptr.ijk,n=1,diff=q_-n;
    while(diff<0) {
        n=n-q_-1;
        ijk_--;
        q_=co[ijk_]-1;
        diff=q_-n;
    }
    ptr.set(ijk_,q_-n);
    return *this;
}

/** Decrements the iterator by one element. */
container_base_2d::iterator container_base_2d::iterator::operator--(int) {
    iterator tmp(*this);
    int q_=ptr.q,ijk_=ptr.ijk,n=1,diff=q_-n;
    while(diff<0) {
        n=n-q_-1;
        ijk_--;
        q_=co[ijk_]-1;
        diff=q_-n;
    }
    ptr.set(ijk_,q_-n);
    return tmp;
}

/** Calculates the number of elements between this iterator and another.
 * \param[in] rhs a reference to another iterator. */
container_base_2d::iterator::difference_type container_base_2d::iterator::operator-(const iterator& rhs) const {
    difference_type diff=0;
    if(ptr.ijk==rhs.ptr.ijk) {
        diff=ptr.q-rhs.ptr.q;
    } else {
        int ijk_small=rhs.ptr.ijk,q_small=rhs.ptr.q,
            ijk_big=ptr.ijk,q_big=ptr.q;
        bool negative=false;
        if(ptr.ijk<rhs.ptr.ijk) {
            negative=true;
            ijk_small=ptr.ijk;q_small=ptr.q;
            ijk_big=rhs.ptr.ijk;q_big=rhs.ptr.q;
        }
        for(int ijk_diff=ijk_small+1;ijk_diff<ijk_big;ijk_diff++) {
            diff+=co[ijk_diff];
        }
        diff=diff+q_big+(co[ijk_small]-q_small);
        if(negative) {diff=-diff;}
    }
    return diff;
}

/** Increments the iterator.
 * \param[in] incre the number of elements to increment by. */
container_base_2d::iterator& container_base_2d::iterator::operator+=(const difference_type& incre) {
    int q_=ptr.q,ijk_=ptr.ijk,
        n=incre,diff=q_+n-co[ijk_];
    while(diff>=0) {
        n=n-co[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co[ijk_];
    }
    ptr.set(ijk_,q_+n);
    return *this;
}

/** Decrements the iterator.
 * \param[in] decre the number of elements to decrement by. */
container_base_2d::iterator& container_base_2d::iterator::operator-=(const difference_type& decre) {
    int q_=ptr.q,ijk_=ptr.ijk,n=decre,diff=q_-n;
    while(diff<0) {
        n=n-q_-1;
        ijk_--;
        q_=co[ijk_]-1;
        diff=q_-n;
    }
    ptr.set(ijk_,q_-n);
    return *this;
}

/* Dereferences the iterator.
 * \param[in] incre the number of elements to offset by. */
c_info& container_base_2d::iterator::operator[](const difference_type& incre) const {
    c_info ci;
    int q_=ptr.q,ijk_=ptr.ijk,n=incre,diff=q_+n-co[ijk_];
    while(diff>=0) {
        n=n-co[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co[ijk_];
    }
    ci.set(ijk_,q_+n);
    return ci;
}

/** Returns an iterator pointing to the first particle in the container.
 * \return The iterator. */
container_base_2d::iterator container_base_2d::begin() {
    return iterator(co);
}

/** Returns an iterator pointing past the last particle in the container.
 * \return The iterator. */
container_base_2d::iterator container_base_2d::end() {
    c_info ci(nxy,0);
    return iterator(co,ci);
}

/** Sets up the class constants to loop over all particles inside a circle.
 * \param[in] (vx,vy) the center of the circle.
 * \param[in] r the radius of the sphere.
 * \param[in] bounds_test whether to do detailed bounds checking. If this is
 *                        false then the class will loop over all particles in
 *                        blocks that overlap the given circle. If it is true,
 *                        the particle will only loop over the particles which
 *                        actually lie within the circle. */
void subset_info_2d::setup_circle(double vx,double vy,double r,bool bounds_test) {
    if(bounds_test) {mode=circle;v0=vx;v1=vy;v2=r*r;} else mode=no_check;
    ai=step_int((vx-ax-r)*xsp);
    bi=step_int((vx-ax+r)*xsp);
    aj=step_int((vy-ay-r)*ysp);
    bj=step_int((vy-ay+r)*ysp);
    setup_common();
}

/** Initializes the class to loop over all particles in a rectangular box.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] bounds_test whether to do detailed bounds checking. If this is
 *                        false then the class will loop over all particles in
 *                        blocks that overlap the given box. If it is true, the
 *                        particle will only loop over the particles which
 *                        actually lie within the box. */
void subset_info_2d::setup_box(double xmin,double xmax,double ymin,double ymax,bool bounds_test) {
    if(bounds_test) {mode=box;v0=xmin;v1=xmax;v2=ymin;v3=ymax;} else mode=no_check;
    ai=step_int((xmin-ax)*xsp);
    bi=step_int((xmax-ax)*xsp);
    aj=step_int((ymin-ay)*ysp);
    bj=step_int((ymax-ay)*ysp);
    setup_common();
}

/** Sets up all of the common constants used for the loop. */
void subset_info_2d::setup_common() {

    // For any non-periodic directions, truncate the block ranges so that they
    // lie within the container grid
    if(!x_prd) {
        if(ai<0) {ai=0;if(bi<0) bi=0;}
        if(bi>=nx) {bi=nx-1;if(ai>=nx) ai=nx-1;}
    }
    if(!y_prd) {
        if(aj<0) {aj=0;if(bj<0) bj=0;}
        if(bj>=ny) {bj=ny-1;if(aj>=ny) aj=ny-1;}
    }

    // Set (di,dj) to be the initial block in the main container grid to
    // consider (taking into account wrapping from periodicity). Set (apx,apy)
    // to be the periodic displacement vector to apply for this initial block.
    di=step_mod(ai,nx);
    apx=step_div(ai,nx)*sx;
    dj=step_mod(aj,ny);
    apy=step_div(aj,ny)*sy;

    // Compute block increment that is frequently used in the loop
    inc1=nx+di-step_mod(bi,nx);

    // Set (ddi,ddj) to be the final block in the main container grid to
    // consider, and set (aapx,aapy) to be the periodic displacement vector to
    // apply for this final block
    ddi=step_mod(bi,nx);
    ddj=step_mod(bj,ny);
    aapx=step_div(bi,nx)*sx;
    aapy=step_div(bj,ny)*sy;
}

/** Moves to the next block, updating all of the required vectors and indices.
 * \param[in,out] ijk_ the index of the block. */
void container_base_2d::iterator_subset::next_block_iter(int &ijk_) {
    if(i<cl_iter->bi) {
        i++;
        if(ci<cl_iter->nx-1) {ci++;ijk_++;} else {ci=0;ijk_+=1-cl_iter->nx;px+=cl_iter->sx;}
    } else if(j<cl_iter->bj) {
        i=cl_iter->ai;ci=cl_iter->di;px=cl_iter->apx;j++;
        if(cj<cl_iter->ny-1) {cj++;ijk_+=cl_iter->inc1;} else {cj=0;ijk_+=cl_iter->inc1-cl_iter->nxy;py+=cl_iter->sy;}
    }
}

/** Moves to the next block, updating all of the required vectors and indices.
 * \param[in,out] ijk_ the index of the block. */
void container_base_2d::iterator_subset::previous_block_iter(int &ijk_) {
    cl_iter->previous_block_iter(ijk_,i,j,ci,cj,px,py);
}

/** Moves to the previous block, updating all of the required vectors and
 * indices.
 * \param[in,out] ijk_ the index of the block.
 * \param[in,out] (i_,j_) the block coordinates.
 * \param[in,out] (ci_,cj_) the block coordinates in the primary grid.
 * \param[in,out] (px_,py_) the periodicity vector. */
void subset_info_2d::previous_block_iter(int &ijk_,int &i_,int &j_,int &ci_,int &cj_,double &px_,double &py_) {
    if(i_>ai) {
        i_--;
        if(ci_>0) {ci_--;ijk_--;} else {ci_=nx-1;ijk_+=nx-1;px_-=sx;}
    } else if(j_>aj) {
        i_=bi;ci_=ddi;px_=aapx;j_--;
        if(cj_>0) {cj_--;ijk_-=inc1;} else {cj_=ny-1;ijk_+=nxy-inc1;py_-=sy;}
    }
}
/** Computes whether the current point is out of bounds, relative to the
 * current loop setup.
 * \param[in] ijk_ the current block.
 * \param[in] q_ the index of the particle in the current block.
 * \param[in] (px_,py_,) the periodicity vector.
 * \return True if the point is out of bounds, false otherwise. */
bool subset_info_2d::out_of_bounds(int ijk_,int q_,double px_,double py_) {
    double *pp=p[ijk_]+ps*q_;
    if(mode==circle) {
        double fx=*pp+px_-v0,fy=pp[1]+py_-v1;
        return fx*fx+fy*fy>v2;
    }
    double f=*pp+px_;if(f<v0||f>v1) return true;
    f=pp[1]+py_;return f<v2||f>v3;
}

/** Initializes the iterator, setting it to point at the first particle in the
 * container.
 * \param[in] si_ a pointer to the information about the particle subset to
 *                consider. */
container_base_2d::iterator_subset::iterator_subset(subset_info_2d* _si) :
    cl_iter(_si), i(cl_iter->ai), j(cl_iter->aj) {

    ci=cl_iter->step_mod(i,cl_iter->nx);
    cj=cl_iter->step_mod(j,cl_iter->ny);

    px=cl_iter->step_div(i,cl_iter->nx)*cl_iter->sx;
    py=cl_iter->step_div(j,cl_iter->ny)*cl_iter->sy;

    int ijk_=ci+cl_iter->nx*cj;
    int q_=0;

    while(cl_iter->co[ijk_]==0) {
        next_block_iter(ijk_);
    }
    while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
        q_++;
        while(q_>=cl_iter->co[ijk_]) {
            q_=0;
            next_block_iter(ijk_);
        }
    }

    ptr.set(ijk_,q_);
}

/** Initializes the iterator.
 * \param[in] si_ a pointer to the information about the particle subset to
 *                consider.
 * \param[in] (ptr_,i_,j_) information about the particle for the iterator
 *                         to point to. */
container_base_2d::iterator_subset::iterator_subset(subset_info_2d* si_,c_info ptr_,int i_,int j_)
  : ptr(ptr_), cl_iter(si_), i(i_), j(j_) {
    ci=cl_iter->step_mod(i,cl_iter->nx);
    cj=cl_iter->step_mod(j,cl_iter->ny);
    px=cl_iter->step_div(i,cl_iter->nx)*cl_iter->sx;
    py=cl_iter->step_div(j,cl_iter->ny)*cl_iter->sy;
}

/** Sets the iterator to equal another.
 * \param[in] other the iterator to copy. */
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator=(iterator_subset other)
{
    cl_iter=other.cl_iter;
    ptr=other.ptr;
    i=other.i;
    j=other.j;
    ci=other.ci;
    cj=other.cj;
    px=other.px;
    py=other.py;
    return *this;
}

/** Increments the iterator by one element. */
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator++() {
    int q_=ptr.q,ijk_=ptr.ijk,n=1;
    while(n>0) {
        q_++;
        while(q_>=cl_iter->co[ijk_]) {
            q_=0;
            next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_++;
            while(q_>=cl_iter->co[ijk_]) {
                q_=0;
                next_block_iter(ijk_);
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}

/** Increments the iterator by one element. */
container_base_2d::iterator_subset container_base_2d::iterator_subset::operator++(int) {
    iterator_subset tmp(*this);
    int q_=ptr.q,ijk_=ptr.ijk,n=1;
    while(n>0) {
        q_++;
        while(q_>=cl_iter->co[ijk_]) {
            q_=0;
            next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_++;
            while(q_>=cl_iter->co[ijk_]) {
                q_=0;
                next_block_iter(ijk_);
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return tmp;
}

/** Decrements the iterator by one element. */
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator--() {
    int q_=ptr.q,ijk_=ptr.ijk,n=1;
    while(n>0) {
        q_--;
        while(q_<0) {
            previous_block_iter(ijk_);
            q_=cl_iter->co[ijk_]-1;
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_--;
            while(q_<0) {
                previous_block_iter(ijk_);
                q_=cl_iter->co[ijk_]-1;
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}

/** Decrements the iterator by one element. */
container_base_2d::iterator_subset container_base_2d::iterator_subset::operator--(int) {
    iterator_subset tmp(*this);
    int q_=ptr.q,ijk_=ptr.ijk,n=1;
    while(n>0) {
        q_--;
        while(q_<0) {
            previous_block_iter(ijk_);
            q_=cl_iter->co[ijk_]-1;
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_--;
            while(q_<0) {
                previous_block_iter(ijk_);
                q_=cl_iter->co[ijk_]-1;
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return tmp;
}

/** Calculates the difference in the number of elements between two iterators.
 * \param[in] rhs the other iterator to compare to.
 * \return The difference. */
container_base_2d::iterator_subset::difference_type container_base_2d::iterator_subset::operator-(const iterator_subset& rhs) const {
    difference_type diff=0;
    if(*this==rhs) {
        diff=0;
    } else if(*this<rhs) {
        iterator_subset tmp(*this);
        while(tmp!=rhs) {
            tmp++;
            diff++;
        }
        diff=-diff;
    } else { //*this>rhs
        iterator_subset tmp(*this);
        while(tmp!=rhs) {
            tmp--;
            diff++;
        }
    }
    return diff;
}

/** Increments the iterator.
 * \param[in] incre the number of elements to increment by. */
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator+=(const difference_type& incre) {
    int q_=ptr.q,ijk_=ptr.ijk,n=incre;
    while(n>0) {
        q_++;
        while(q_>=cl_iter->co[ijk_]) {
            q_=0;
            next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_++;
            while(q_>=cl_iter->co[ijk_]) {
                q_=0;
                next_block_iter(ijk_);
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}

/** Decrements the iterator.
 * \param[in] decre the number of elements to decrement by. */
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator-=(const difference_type& decre) {
    int q_=ptr.q,ijk_=ptr.ijk,n=decre;
    while(n>0) {
        q_--;
        while(q_<0) {
            previous_block_iter(ijk_);
            q_=cl_iter->co[ijk_]-1;
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_--;
            while(q_<0) {
                previous_block_iter(ijk_);
                q_=cl_iter->co[ijk_]-1;
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}

/* Dereferences the iterator.
 * \param[in] incre the number of elements to offset by. */
c_info& container_base_2d::iterator_subset::operator[](const difference_type& incre) const {
    c_info ci;
    iterator_subset tmp(*this);
    int q_=ptr.q,ijk_=ptr.ijk,n=incre;
    while(n>0) {
        q_++;
        while(q_>=cl_iter->co[ijk_]) {
            q_=0;
            tmp.next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)) {
            q_++;
            while(q_>=cl_iter->co[ijk_]) {
                q_=0;
                tmp.next_block_iter(ijk_);
            }
        }
        n--;
    }
    ci.set(ijk_,q_);
    return ci;
}

/** Returns an iterator pointing to the first particle in the container.
 * \return The iterator. */
container_base_2d::iterator_subset container_base_2d::begin(subset_info_2d& si) {
    return iterator_subset(&si);
}

/** Returns an iterator pointing past the last particle in the container.
 * \return The iterator. */
container_base_2d::iterator_subset container_base_2d::end(subset_info_2d& si) {
    c_info cinfo;
    //find the last particle to point to
    int i_=si.bi,j_=si.bj,
        ci_=si.ddi,cj_=si.ddj,
        ijk_=si.ddi+si.nx*si.ddj,
        q_=si.co[ijk_]-1;
    double px_=si.aapx,py_=si.aapy;

    while(q_<0) {
        si.previous_block_iter(ijk_,i_,j_,ci_,cj_,px_,py_);
        q_=si.co[ijk_]-1;
    }
    while(si.mode!=no_check&&si.out_of_bounds(ijk_,q_,px_,py_)) {
        q_--;
        while(q_<0) {
            si.previous_block_iter(ijk_,i_,j_,ci_,cj_,px_,py_);
            q_=si.co[ijk_]-1;
        }
    }
    cinfo.set(ijk_,q_+1);
    return iterator_subset(&si,cinfo,i_,j_);
}

/** Increments the iterator by one element. */
container_base_2d::iterator_order& container_base_2d::iterator_order::operator++() {
    ptr_n++;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/** Increments the iterator by one element. */
container_base_2d::iterator_order container_base_2d::iterator_order::operator++(int) {
    iterator_order tmp(*this);
    ptr_n++;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return tmp;
}

/** Decrements the iterator by one element. */
container_base_2d::iterator_order& container_base_2d::iterator_order::operator--() {
    ptr_n--;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/** Decrements the iterator by one element. */
container_base_2d::iterator_order container_base_2d::iterator_order::operator--(int) {
    iterator_order tmp(*this);
    ptr_n--;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return tmp;
}

//Supports the arithmetic operators + and - between an iterator and an integer value, or subtracting an iterator from another.
container_base_2d::iterator_order::difference_type container_base_2d::iterator_order::operator-(const iterator_order& rhs) const {
    difference_type diff=ptr_n-rhs.ptr_n;
    return diff;
}

/** Increments the iterator.
 * \param[in] incre the number of elements to increment by. */
container_base_2d::iterator_order& container_base_2d::iterator_order::operator+=(const difference_type& incre) {
    ptr_n+=incre;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/** Decrements the iterator.
 * \param[in] decre the number of elements to decrement by. */
container_base_2d::iterator_order& container_base_2d::iterator_order::operator-=(const difference_type& decre) {
    ptr_n-=decre;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/* Dereferences the iterator.
 * \param[in] incre the number of elements to offset by. */
c_info& container_base_2d::iterator_order::operator[](const difference_type& incre) const {
    int ci_n=ptr_n+incre;
    c_info ci(cp_iter[2*ci_n],cp_iter[2*ci_n+1]);
    return ci;
}

/** Returns an iterator pointing to the first particle in the container.
 * \return The iterator. */
container_base_2d::iterator_order container_base_2d::begin(particle_order &vo) {return iterator_order(vo);}

/** Returns an iterator pointing past the last particle in the container.
 * \return The iterator. */
container_base_2d::iterator_order container_base_2d::end(particle_order &vo) {    //vo, ptr, n
    int ptr_n_=0.5*(vo.op-vo.o);//1-over-the-last-particle, eg. if 0,1,2,3,4 particle, here, ptr_n=5
    return iterator_order(vo, c_info(nxy,0), ptr_n_);
}

}
