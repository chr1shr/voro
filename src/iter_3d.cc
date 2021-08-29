// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#include "iter_3d.hh"

namespace voro {

/** Initializes the iterator, setting it to point at the first particle in the
 * container.
 * \param[in] co_ a pointer to the particle count array. */
container_base_3d::iterator::iterator(int* co_) : co(co_) {

    // Find the first empty block
    int ijk=0;
    while(co[ijk]==0) ijk++;
    // XXX CHR - what about an empty container? Do we need to check for ijk out of range?
    ptr.set(ijk,0);
}

/** Increments the iterator by one element. */
container_base_3d::iterator& container_base_3d::iterator::operator++() {
    int q_=ptr.q,ijk_=ptr.ijk,n=1,diff=q_+n-co[ijk_];
    // XXX CHR - I think the number of arithmetic operations in the loop below
    // can be reduced. q_ is being set multiple times to zero, and both diff
    // and n are recalculated.
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
container_base_3d::iterator container_base_3d::iterator::operator++(int) {
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
container_base_3d::iterator& container_base_3d::iterator::operator--() {
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    int diff=q_-n;
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
container_base_3d::iterator container_base_3d::iterator::operator--(int) {
    iterator tmp(*this);
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    int diff=q_-n;
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
container_base_3d::iterator::difference_type container_base_3d::iterator::operator-(const iterator& rhs) const {
    difference_type diff=0;
    if(ptr.ijk==rhs.ptr.ijk) {
        if(ptr.q==rhs.ptr.q) { // XXX CHR - simplify these six lines to just diff=pts.q-rhs.ptr.q ?
            diff=0;
        }
        else{
            diff=ptr.q-rhs.ptr.q;
        }
    }
    else{
        int ijk_small=rhs.ptr.ijk;int q_small=rhs.ptr.q;
        int ijk_big=ptr.ijk;int q_big=ptr.q;
        bool negative=false;
        if(ptr.ijk < rhs.ptr.ijk) {
            negative=true;
            ijk_small=ptr.ijk;q_small=ptr.q;
            ijk_big=rhs.ptr.ijk;q_big=rhs.ptr.q;
        }
        for(int ijk_diff=ijk_small+1;ijk_diff<ijk_big; ijk_diff++) {
            diff+=co[ijk_diff];
        }
        diff=diff+q_big+(co[ijk_small]-q_small);
        if(negative) {diff=-diff;}
    }
    return diff;
}

/** Increments the iterator.
 * \param[in] incre the number of elements to increment by. */
container_base_3d::iterator& container_base_3d::iterator::operator+=(const difference_type& incre) {
    int q_=ptr.q,ijk_=ptr.ijk,n=incre,diff=q_+n-co[ijk_];
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
container_base_3d::iterator& container_base_3d::iterator::operator-=(const difference_type& decre) {
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
c_info& container_base_3d::iterator::operator[](const difference_type& incre) const {
    c_info ci;
    int q_=ptr.q,ijk_=ptr.ijk,n=incre,diff=q_+n-co[ijk_];
    // XXX CHR - the following code assumes that incre is positive. For normal
    // array lookups that it not necessary: writing something like a[-2] is
    // fine, and means *(a-2). Do we need to take into account if incre is
    // negative?
    while(diff>=0) {
        n=n-co[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co[ijk_];
    }
    ci.set(ijk_,q_+n);
    return ci;
    // XXX CHR - I think there is a conceptual issue here with this function.
    // We are returning a reference to a local copy of ci, which will disappear
    // once the function finishes.
}

/** Returns an iterator pointing to the first particle in the container.
 * \return The iterator. */
container_base_3d::iterator container_base_3d::begin() {return iterator(co);}

/** Returns an iterator pointing past the last particle in the container.
 * \return The iterator. */
container_base_3d::iterator container_base_3d::end() {
    c_info ci;
    //find the last particle to point to
    int ijk_=nxyz-1;
    while(co[ijk_]==0) {
        ijk_--;
    }
    int q_=co[ijk_];  //1 over the end of the particles
    ci.set(ijk_,q_);
    return iterator(co, ci);

    // XXX CHR - here you are scanning backward through the blocks to put
    // "end()" after the last particle, at (ijk,q+1). But if there are n
    // blocks, then couldn't we just put "end()" at (n,0)? That wouldn't
    // require a scan. (Same issue for iterator_subset below.)
}

//swappable??
// XXX CHR - the library is designed to compile with the original ANSI C++
// standard, and the swappable property is not needed. However some people may
// wish to compile the library with C++ 11, and the swappable property is
// required. To ensure interoperability with C++ 11 we should check that the
// iterators can be swapped (using swap(a,b)).

/** Sets up the class constants to loop over all particles inside a sphere.
 * \param[in] (vx,vy,vz) the center of the sphere.
 * \param[in] r the radius of the sphere.
 * \param[in] bounds_test whether to do detailed bounds checking. If this is
 *                        false then the class will loop over all particles in
 *                        blocks that overlap the given sphere. If it is true,
 *                        the particle will only loop over the particles which
 *                        actually lie within the sphere. */
void subset_info_3d::setup_sphere(double vx,double vy,double vz,double r,bool bounds_test) {
    if(bounds_test) {mode=sphere;v0=vx;v1=vy;v2=vz;v3=r*r;} else mode=no_check;
    ai=step_int((vx-ax-r)*xsp);
    bi=step_int((vx-ax+r)*xsp);
    aj=step_int((vy-ay-r)*ysp);
    bj=step_int((vy-ay+r)*ysp);
    ak=step_int((vz-az-r)*zsp);
    bk=step_int((vz-az+r)*zsp);
    setup_common();
}

/** Initializes the class to loop over all particles in a rectangular box.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box.
 * \param[in] bounds_test whether to do detailed bounds checking. If this is
 *                        false then the class will loop over all particles in
 *                        blocks that overlap the given box. If it is true, the
 *                        particle will only loop over the particles which
 *                        actually lie within the box. */
void subset_info_3d::setup_box(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,bool bounds_test) {
    if(bounds_test) {mode=box;v0=xmin;v1=xmax;v2=ymin;v3=ymax;v4=zmin;v5=zmax;} else mode=no_check;
    ai=step_int((xmin-ax)*xsp);
    bi=step_int((xmax-ax)*xsp);
    aj=step_int((ymin-ay)*ysp);
    bj=step_int((ymax-ay)*ysp);
    ak=step_int((zmin-az)*zsp);
    bk=step_int((zmax-az)*zsp);
    setup_common();
}

/** Sets up all of the common constants used for the loop. */
void subset_info_3d::setup_common() {

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
    if(!z_prd) {
            if(ak<0) {ak=0;if(bk<0) bk=0;}
            if(bk>=nz) {bk=nz-1;if(ak>=nz) ak=nz-1;}
    }

    // Set (di,dj,dk) to be the initial block in the main container grid to
    // consider (taking into account wrapping from periodicity). Set
    // (apx,apy,apz) to be the periodic displacement vector to apply for this
    // initial block.
    di=step_mod(ai,nx);
    apx=step_div(ai,nx)*sx;
    dj=step_mod(aj,ny);
    apy=step_div(aj,ny)*sy;
    dk=step_mod(ak,nz);
    apz=step_div(ak,nz)*sz;

    // Compute block increments that are frequently used during the loop
    inc1=di-step_mod(bi,nx);
    inc2=nx*(ny+dj-step_mod(bj,ny))+inc1;
    inc1+=nx;

    // Set (ddi,ddj,ddk) to be the final block in the main container grid to
    // consider, and set (aapx,aapy,aapz) to be the periodic displacement
    // vector to apply for this final block
    ddi=step_mod(bi,nx);
    ddj=step_mod(bj,ny);
    ddk=step_mod(bk,nz);
    aapx=step_div(bi,nx)*sx;
    aapy=step_div(bj,ny)*sy;
    aapz=step_div(bk,nz)*sz;
}

// XXX CHR - I don't think it is necessary to append "_iter" to the end of
// function/variable names. If functions/variables in two different classes
// represent the same thing, then it is fine (and likely preferable) for them
// to have the same name. C++ name scoping avoids potential clashes.

/** Moves to the previous block, updating all of the required vectors and
 * indices.
 * \param[in,out] ijk_ the index of the block.
 * \param[in,out] (i_,j_,k_) the block coordinates.
 * \param[in,out] (ci_,cj_,ck_) the block coordinates in the primary grid.
 * \param[in,out] (px_,py_,pz_) the periodicity vector. */
void subset_info_3d::previous_block_iter(int &ijk_,int &i_,int &j_,int &k_,int &ci_,int &cj_,int &ck_,double &px_,double &py_,double &pz_) {
    if(i_>ai) {
        i_--;
        if(ci_>0) {ci_--;ijk_--;} else {ci_=nx-1;ijk_+=nx-1;px_-=sx;}
    } else if(j_>aj) {
        i_=bi;ci_=ddi;px_=aapx;j_--;
        if(cj_>0) {cj_--;ijk_-=inc1;} else {cj_=ny-1;ijk_+=nxy-inc1;py_-=sy;}
    } else if(k_>ak) {
        i_=bi;ci_=ddi;px_=aapx;
        j_=bj;cj_=ddj;py_=aapy;k_--;
        if(ck_>0) {ck_--;ijk_-=inc2;} else {ck_=nz-1;ijk_+=nxyz-inc2;pz_-=sz;}
    }
}

/** Computes whether the current point is out of bounds, relative to the
 * current loop setup.
 * \param[in] ijk_ the current block.
 * \param[in] q_ the index of the particle in the current block.
 * \param[in] (px_,py_,pz_) the periodicity vector.
 * \return True if the point is out of bounds, false otherwise. */
bool subset_info_3d::out_of_bounds(int ijk_,int q_,double px_,double py_,double pz_) {
    double *pp=p[ijk_]+ps*q_;
    if(mode==sphere) {
            double fx=*pp+px_-v0,fy=pp[1]+py_-v1,fz=pp[2]+pz_-v2;
            return fx*fx+fy*fy+fz*fz>v3;
    }
    double f=*pp+px_;if(f<v0||f>v1) return true;
    f=pp[1]+py_;if(f<v2||f>v3) return true;
    f=pp[2]+pz_;return f<v4||f>v5;
}

/** Moves to the next block, updating all of the required vectors and indices.
 * \param[in,out] ijk_ the index of the block. */
void container_base_3d::iterator_subset::next_block_iter(int &ijk_) {
    // XXX CHR - it's not clear to me that this function needs to take ijk_ as
    // an argument. It already has access to ijk - it is ptr.ijk. You could
    // even write "int &ijk=ptr.ijk;" at the start to make a shorthand to this
    // variable.
    if(i<cl_iter->bi) {
        i++;
        if(ci<cl_iter->nx-1) {ci++;ijk_++;}
        else {ci=0;ijk_+=1-cl_iter->nx;px+=cl_iter->sx;}
    } else if(j<cl_iter->bj) {
        i=cl_iter->ai;ci=cl_iter->di;px=cl_iter->apx;j++;
        if(cj<cl_iter->ny-1) {cj++;ijk_+=cl_iter->inc1;}
        else {cj=0;ijk_+=cl_iter->inc1-cl_iter->nxy;py+=cl_iter->sy;}
    } else if(k<cl_iter->bk) {
        i=cl_iter->ai;ci=cl_iter->di;j=cl_iter->aj;cj=cl_iter->dj;px=cl_iter->apx;py=cl_iter->apy;k++;
        if(ck<cl_iter->nz-1) {ck++;ijk_+=cl_iter->inc2;}
        else {ck=0;ijk_+=cl_iter->inc2-cl_iter->nxyz;pz+=cl_iter->sz;}
    }
}

/** Moves to the next block, updating all of the required vectors and indices.
 * \param[in,out] ijk_ the index of the block. */
void container_base_3d::iterator_subset::previous_block_iter(int &ijk_) {
    cl_iter->previous_block_iter(ijk_,i,j,k,ci,cj,ck,px,py,pz);
    // XXX CHR - why is the previous_block_iter routine done within
    // subset_info_3d, but the next_block_iter routine is done within the iterator
    // itself? I think the next_block_iter approach is probably better, since
    // it can operate on its own data (i, j, k, etc.) rather than having to pass
    // all of them by reference to subset_info_3d.
}

/** Initializes the iterator, setting it to point at the first particle in the
 * container.
 * \param[in] si_ a pointer to the information about the particle subset to
 *                consider. */
container_base_3d::iterator_subset::iterator_subset(subset_info_3d* si_)
    : cl_iter(si_), i(cl_iter->ai), j(cl_iter->aj), k(cl_iter->ak) {

    ci=cl_iter->step_mod(i,cl_iter->nx);
    cj=cl_iter->step_mod(j,cl_iter->ny);
    ck=cl_iter->step_mod(k,cl_iter->nz);

    px=cl_iter->step_div(i,cl_iter->nx)*cl_iter->sx;
    py=cl_iter->step_div(j,cl_iter->ny)*cl_iter->sy;
    pz=cl_iter->step_div(k,cl_iter->nz)*cl_iter->sz;

    int ijk_=ci+cl_iter->nx*(cj+cl_iter->ny*ck);
    int q_=0;

    while(cl_iter->co[ijk_]==0) {
        // XXX CHR - need to catch case when container is empty and ijk goes
        // out of range? (Same issue as for previous class.)
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
 * \param[in] (ptr_,i_,j_,k_) information about the particle for the iterator
 *                            to point to. */
container_base_3d::iterator_subset::iterator_subset(subset_info_3d* si_,c_info ptr_,int i_,int j_,int k_)
    : ptr(ptr_), cl_iter(si_), i(i_), j(j_), k(k_) {
    ci=cl_iter->step_mod(i,cl_iter->nx);
    cj=cl_iter->step_mod(j,cl_iter->ny);
    ck=cl_iter->step_mod(k,cl_iter->nz);
    px=cl_iter->step_div(i,cl_iter->nx)*cl_iter->sx;
    py=cl_iter->step_div(j,cl_iter->ny)*cl_iter->sy;
    pz=cl_iter->step_div(k,cl_iter->nz)*cl_iter->sz;
}

/** Sets the iterator to equal another.
 * \param[in] other the iterator to copy. */
container_base_3d::iterator_subset& container_base_3d::iterator_subset::operator=(iterator_subset other) {
    cl_iter=other.cl_iter;
    ptr=other.ptr;
    i=other.i;
    j=other.j;
    k=other.k;
    ci=other.ci;
    cj=other.cj;
    ck=other.ck;
    px=other.px;
    py=other.py;
    pz=other.pz;
    return *this;
}

/** Increments the iterator by one element. */
container_base_3d::iterator_subset& container_base_3d::iterator_subset::operator++() {
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
container_base_3d::iterator_subset container_base_3d::iterator_subset::operator++(int) {
    iterator_subset tmp=*this;
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
container_base_3d::iterator_subset& container_base_3d::iterator_subset::operator--() {
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
    // XXX CHR - What happens if "--" is applied when you are at the first
    // element? Is the iterator meant to handle that case?
    ptr.set(ijk_,q_);
    return *this;
}

/** Decrements the iterator by one element. */
container_base_3d::iterator_subset container_base_3d::iterator_subset::operator--(int) {
    iterator_subset tmp=*this;
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
container_base_3d::iterator_subset::difference_type container_base_3d::iterator_subset::operator-(const iterator_subset& rhs) const {
    difference_type diff=0;
    if(*this==rhs) {
        diff=0;
    }
    else if(*this<rhs) {
        // XXX CHR - This is a fairly expensive way to compute the difference.
        // You make a copy of the current iterator, and then you step
        // individually through each particle. I guess, though, that with the
        // out_of_bounds routine, it is difficult to avoid this.
        iterator_subset tmp(*this);
        while(tmp!=rhs) {
            tmp++;
            diff++;  // XXX CHR - change to diff--
        }
        diff=-diff;  // XXX CHR - and delete this line
    }
    else {
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
container_base_3d::iterator_subset& container_base_3d::iterator_subset::operator+=(const difference_type& incre) {
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
container_base_3d::iterator_subset& container_base_3d::iterator_subset::operator-=(const difference_type& decre) {
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
c_info& container_base_3d::iterator_subset::operator[](const difference_type& incre) const {
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
container_base_3d::iterator_subset container_base_3d::begin(subset_info_3d& si) {
    return iterator_subset(&si);
}

/** Returns an iterator pointing past the last particle in the container.
 * \return The iterator. */
container_base_3d::iterator_subset container_base_3d::end(subset_info_3d& si) {
    c_info cinfo;
    //find the last particle to point to
    int i_=si.bi,j_=si.bj,k_=si.bk,
        ci_=si.ddi,cj_=si.ddj,ck_=si.ddk,
        ijk_=si.ddi+si.nx*(si.ddj+si.ny*si.ddk),
        q_=si.co[ijk_]-1;
    double px_=si.aapx,py_=si.aapy,pz_=si.aapz;

    while(q_<0) {
        si.previous_block_iter(ijk_,i_,j_,k_,ci_,cj_,ck_,px_,py_,pz_);
        q_=si.co[ijk_]-1;
    }
    while(si.mode!=no_check&&si.out_of_bounds(ijk_,q_,px_,py_,pz_)) {
        q_--;
        while(q_<0) {
            si.previous_block_iter(ijk_,i_,j_,k_,ci_,cj_,ck_,px_,py_,pz_);
            q_=si.co[ijk_]-1;
        }
    }
    cinfo.set(ijk_,q_+1);
    return iterator_subset(&si,cinfo,i_,j_,k_);
}

/** Increments the iterator by one element. */
container_base_3d::iterator_order& container_base_3d::iterator_order::operator++() {
    ptr_n++;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/** Increments the iterator by one element. */
container_base_3d::iterator_order container_base_3d::iterator_order::operator++(int) {
    iterator_order tmp(*this);
    ptr_n++;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return tmp;
}

/** Decrements the iterator by one element. */
container_base_3d::iterator_order& container_base_3d::iterator_order::operator--() {
    ptr_n--;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/** Decrements the iterator by one element. */
container_base_3d::iterator_order container_base_3d::iterator_order::operator--(int) {
    iterator_order tmp(*this);
    ptr_n--;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return tmp;
}

/** Calculates the number of elements between this iterator and another.
 * \param[in] rhs a reference to another iterator. */
container_base_3d::iterator_order::difference_type container_base_3d::iterator_order::operator-(const iterator_order& rhs) const {
    difference_type diff=ptr_n-rhs.ptr_n;
    return diff;
}

/** Increments the iterator.
 * \param[in] incre the number of elements to increment by. */
container_base_3d::iterator_order& container_base_3d::iterator_order::operator+=(const difference_type& incre) {
    ptr_n+=incre;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/** Decrements the iterator.
 * \param[in] decre the number of elements to decrement by. */
container_base_3d::iterator_order& container_base_3d::iterator_order::operator-=(const difference_type& decre) {
    ptr_n-=decre;
    ptr.set(cp_iter[2*ptr_n],cp_iter[2*ptr_n+1]);
    return *this;
}

/* Dereferences the iterator.
 * \param[in] incre the number of elements to offset by. */
c_info& container_base_3d::iterator_order::operator[](const difference_type& incre) const {
    int ci_n=ptr_n+incre;
    c_info ci(cp_iter[2*ci_n],cp_iter[2*ci_n+1]);
    return ci;
}

/** Returns an iterator pointing to the first particle in the container.
 * \return The iterator. */
container_base_3d::iterator_order container_base_3d::begin(particle_order &vo) {
    return iterator_order(vo);
}

/** Returns an iterator pointing past the last particle in the container.
 * \return The iterator. */
container_base_3d::iterator_order container_base_3d::end(particle_order &vo) {
    int ptr_n_=(vo.op-vo.o)/2; //1-over-the-last-particle, eg. if 0,1,2,3,4 particle, here, ptr_n=5
    c_info ci;
    int ijk_=-1; //dummy
    int q_=-1; //dummy
    ci.set(ijk_,q_);
    // XXX CHR - Do we need to set dummy values here? [Also, if we do need to set them, then presumably
    // you can write ci.set(-1,-1) .]
    return iterator_order(vo,ci,ptr_n_);
}

}
