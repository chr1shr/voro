// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#include "iter_2d.hh"

//default-constructible??
container_base_2d::iterator::iterator() : co_iter(0){}
//constructor: point to first particle in container
container_base_2d::iterator::iterator(int* _co) : co_iter(_co) {
    //find the first particle to point to
    int ijk_=0; int q_=0;
    while(co_iter[ijk_]==0){
        ijk_++;
    }
    ptr.set(ijk_,q_);
}
//constructor: point to a particle in container
container_base_2d::iterator::iterator(int* _co, c_info _ptr) : ptr(_ptr), co_iter(_co) {}
//copy-constructible
container_base_2d::iterator::iterator(const iterator& ci) : ptr(ci.ptr), co_iter(ci.co_iter) {}
//copy-assignable
container_base_2d::iterator& container_base_2d::iterator::operator=(iterator other){
    co_iter=other.co_iter;
    ptr=other.ptr;
    return *this;
}
//destructible??
container_base_2d::iterator::~iterator(){}

//Can be compared for equivalence using the equality/inequality operators
bool container_base_2d::iterator::operator==(const iterator& rhs) const {
    if(ptr.ijk==rhs.ptr.ijk && ptr.q==rhs.ptr.q){return true;}
    else{return false;}
}
bool container_base_2d::iterator::operator!=(const iterator& rhs) const {
    if(ptr.ijk==rhs.ptr.ijk && ptr.q==rhs.ptr.q){return false;}
    else{return true;}
}

//Can be dereferenced as an rvalue (if in a dereferenceable state).
c_info& container_base_2d::iterator::operator*() {
    return ptr;
}
c_info* container_base_2d::iterator::operator->() {
    return &ptr;
}

//Can be incremented (if in a dereferenceable state).
//The result is either also dereferenceable or a past-the-end iterator.
//Two iterators that compare equal, keep comparing equal after being both increased.
container_base_2d::iterator& container_base_2d::iterator::operator++() {
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    int diff=q_+n-co_iter[ijk_];
    while(diff>=0){
        n=n-co_iter[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co_iter[ijk_];
    }
    ptr.set(ijk_,q_+n);
    return *this;
}
container_base_2d::iterator container_base_2d::iterator::operator++(int) {
    iterator tmp(*this);
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    int diff=q_+n-co_iter[ijk_];
    while(diff>=0){
        n=n-co_iter[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co_iter[ijk_];
    }
    ptr.set(ijk_,q_+n);
    return tmp;
}

//Can be decremented (if a dereferenceable iterator value precedes it).
container_base_2d::iterator& container_base_2d::iterator::operator--() {
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    int diff=q_-n;
    while(diff<0){
        n=n-q_-1;
        ijk_--;
        q_=co_iter[ijk_]-1;
        diff=q_-n;
    }
    ptr.set(ijk_,q_-n);
    return *this;
}
container_base_2d::iterator container_base_2d::iterator::operator--(int) {
    iterator tmp(*this);
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    int diff=q_-n;
    while(diff<0){
        n=n-q_-1;
        ijk_--;
        q_=co_iter[ijk_]-1;
        diff=q_-n;
    }
    ptr.set(ijk_,q_-n);
    return tmp;
}

//Supports the arithmetic operators + and - between an iterator and an integer value, or subtracting an iterator from another.
container_base_2d::iterator::difference_type container_base_2d::iterator::operator-(const iterator& rhs) const {
    difference_type diff=0;
    if(ptr.ijk==rhs.ptr.ijk){
        if(ptr.q==rhs.ptr.q){
            diff=0;
        }
        else{
            diff=ptr.q-rhs.ptr.q;
        }
    }
    else{
        int ijk_small=rhs.ptr.ijk; int q_small=rhs.ptr.q;
        int ijk_big=ptr.ijk; int q_big=ptr.q;
        bool negative=false;
        if(ptr.ijk < rhs.ptr.ijk){
            negative=true;
            ijk_small=ptr.ijk; q_small=ptr.q;
            ijk_big=rhs.ptr.ijk; q_big=rhs.ptr.q;
        }
        for(int ijk_diff=ijk_small+1; ijk_diff<ijk_big; ijk_diff++){
            diff+=co_iter[ijk_diff];
        }
        diff=diff+q_big+(co_iter[ijk_small]-q_small);
        if(negative==true){diff=-diff;}
    }
    return diff;
}

container_base_2d::iterator container_base_2d::iterator::operator+(const difference_type& incre) const {
    iterator tmp(*this);
    tmp+=incre;
    return tmp;
}
container_base_2d::iterator container_base_2d::iterator::operator-(const difference_type& decre) const {
    iterator tmp(*this);
    tmp-=decre;
    return tmp;
}

//Can be compared with inequality relational operators (<, >, <= and >=).
bool container_base_2d::iterator::operator>(const iterator& rhs) const {
    if(ptr.ijk > rhs.ptr.ijk){return true;}
    else if(ptr.ijk==rhs.ptr.ijk && ptr.q > rhs.ptr.q){return true;}
    else{return false;}
}
bool container_base_2d::iterator::operator<(const iterator& rhs) const {
    if(ptr.ijk < rhs.ptr.ijk){return true;}
    else if(ptr.ijk==rhs.ptr.ijk && ptr.q < rhs.ptr.q){return true;}
    else{return false;}
}
bool container_base_2d::iterator::operator>=(const iterator& rhs) const {
    if(ptr.ijk > rhs.ptr.ijk){return true;}
    else if(ptr.ijk==rhs.ptr.ijk && ptr.q >= rhs.ptr.q){return true;}
    else{return false;}
}
bool container_base_2d::iterator::operator<=(const iterator& rhs) const {
    if(ptr.ijk < rhs.ptr.ijk){return true;}
    else if(ptr.ijk==rhs.ptr.ijk && ptr.q <= rhs.ptr.q){return true;}
    else{return false;}
}

//Supports compound assignment operations += and -=
container_base_2d::iterator& container_base_2d::iterator::operator+=(const difference_type& incre) {
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=incre;
    int diff=q_+n-co_iter[ijk_];
    while(diff>=0){
        n=n-co_iter[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co_iter[ijk_];
    }
    ptr.set(ijk_,q_+n);
    return *this;
}
container_base_2d::iterator& container_base_2d::iterator::operator-=(const difference_type& decre) {
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=decre;
    int diff=q_-n;
    while(diff<0){
        n=n-q_-1;
        ijk_--;
        q_=co_iter[ijk_]-1;
        diff=q_-n;
    }
    ptr.set(ijk_,q_-n);
    return *this;
}

//Supports the offset dereference operator ([])
c_info& container_base_2d::iterator::operator[](const difference_type& incre) const {
    c_info ci;
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=incre;
    int diff=q_+n-co_iter[ijk_];
    while(diff>=0){
        n=n-co_iter[ijk_]+q_;
        ijk_++;
        q_=0;
        diff=q_+n-co_iter[ijk_];
    }
    ci.set(ijk_,q_+n);
    return ci;
}

//swappable??

container_base_2d::iterator container_base_2d::begin() {
    return iterator(co); //first particle in the container
}
container_base_2d::iterator container_base_2d::end() {
    c_info ci;
    //find the last particle to point to
    int ijk_=nxy-1;
    while(co[ijk_]==0){
        ijk_--;
    }
    int q_=co[ijk_];  //1 over the end of the particles
    ci.set(ijk_,q_);
    return iterator(co, ci);
}

//-------------------------iterator_subset--------------------------------//
void subset_info::setup_circle(double vx,double vy,double r,bool bounds_test) {
	if(bounds_test) {mode=circle;v0=vx;v1=vy;v2=r*r; printf("circle\n");} else mode=no_check;
	ai=step_int((vx-ax-r)*xsp);
	bi=step_int((vx-ax+r)*xsp);
	aj=step_int((vy-ay-r)*ysp);
	bj=step_int((vy-ay+r)*ysp);
	setup_common();
}

void subset_info::setup_box(double xmin,double xmax,double ymin,double ymax,bool bounds_test) {
	if(bounds_test) {mode=rectangle;v0=xmin;v1=xmax;v2=ymin;v3=ymax;} else mode=no_check;
	ai=step_int((xmin-ax)*xsp);
	bi=step_int((xmax-ax)*xsp);
	aj=step_int((ymin-ay)*ysp);
	bj=step_int((ymax-ay)*ysp);
	setup_common();
}

void subset_info::setup_intbox(int ai_,int bi_,int aj_,int bj_) {
	ai=ai_;bi=bi_;aj=aj_;bj=bj_;
	mode=no_check;
	setup_common();
}

void subset_info::setup_common() {
	if(!xperiodic) {
		if(ai<0) {ai=0;if(bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if(ai>=nx) ai=nx-1;}
	}
	if(!yperiodic) {
		if(aj<0) {aj=0;if(bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if(aj>=ny) aj=ny-1;}
	}
	di=step_mod(ai,nx);
        apx=step_div(ai,nx)*sx;
	dj=step_mod(aj,ny);
        apy=step_div(aj,ny)*sy;
	inc1=nx+di-step_mod(bi,nx);

        ddi=step_mod(bi,nx);
        ddj=step_mod(bj,ny);
        aapx=step_div(bi,nx)*sx;
        aapy=step_div(bj,ny)*sy;
}

bool container_base_2d::iterator_subset::out_of_bounds_iter(int ijk_, int q_) const {
    return cl_iter->out_of_bounds_2(ijk_,q_,px,py);
}

void container_base_2d::iterator_subset::next_block_iter(int &ijk_){
    if(i<cl_iter->bi) {
        i++;
        if(ci<cl_iter->nx-1) {ci++;ijk_++;} else {ci=0;ijk_+=1-cl_iter->nx;px+=cl_iter->sx;}
    } else if(j<cl_iter->bj) {
        i=cl_iter->ai;ci=cl_iter->di;px=cl_iter->apx;j++;
        if(cj<cl_iter->ny-1) {cj++;ijk_+=cl_iter->inc1;} else {cj=0;ijk_+=cl_iter->inc1-cl_iter->nxy;py+=cl_iter->sy;}
    }
}

void container_base_2d::iterator_subset::previous_block_iter(int &ijk_){
    cl_iter->previous_block_iter2(ijk_,i,j,ci,cj,px,py);
}

void subset_info::previous_block_iter2(int &ijk_, int &i_, int &j_, int &ci_, int &cj_, double &px_, double &py_){
    if(i_>ai){
        i_--;
        if(ci_>0){ci_--; ijk_--;} else {ci_=nx-1;ijk_+=nx-1;px_-=sx;}
    } else if(j_>aj) {
        i_=bi;ci_=ddi;px_=aapx;j_--;
        if(cj_>0){cj_--; ijk_-=inc1;} else {cj_=ny-1;ijk_+=nxy-inc1;py_-=sy;}
    }
}

bool subset_info::out_of_bounds_2(int ijk_, int q_, double px_, double py_){
    double *pp=p[ijk_]+ps*q_;
    if(mode==circle) {
            double fx(*pp+px_-v0),fy(pp[1]+py_-v1);
            return fx*fx+fy*fy>v2;
    } else {
            double f(*pp+px_);if(f<v0||f>v1) return true;
            f=pp[1]+py_;return f<v2||f>v3;
    }
}

container_base_2d::iterator_subset::iterator_subset(subset_info* _si) :
            cl_iter(_si), j(cl_iter->aj), i(cl_iter->ai)
{
    //find the first particle to point to
    ci=cl_iter->step_mod(i,cl_iter->nx);
    cj=cl_iter->step_mod(j,cl_iter->ny);

    px=cl_iter->step_div(i,cl_iter->nx)*cl_iter->sx;
    py=cl_iter->step_div(j,cl_iter->ny)*cl_iter->sy;

    int ijk_=ci+cl_iter->nx*cj;
    int q_=0;

    while(cl_iter->co[ijk_]==0){
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

container_base_2d::iterator_subset::iterator_subset(subset_info* _si, c_info _ptr, int _i, int _j)
  : cl_iter(_si), ptr(_ptr),i(_i),j(_j){
    ci=cl_iter->step_mod(i,cl_iter->nx);
    cj=cl_iter->step_mod(j,cl_iter->ny);
    px=cl_iter->step_div(i,cl_iter->nx)*cl_iter->sx;
    py=cl_iter->step_div(j,cl_iter->ny)*cl_iter->sy;
}
container_base_2d::iterator_subset::iterator_subset(const iterator_subset& _ci)
  : ptr(_ci.ptr), cl_iter(_ci.cl_iter), i(_ci.i), j(_ci.j),
    ci(_ci.ci), cj(_ci.cj), px(_ci.px), py(_ci.py) {}

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

//Can be compared for equivalence using the equality/inequality operators
bool container_base_2d::iterator_subset::operator==(const iterator_subset& rhs) const
{
    if(ptr.ijk==rhs.ptr.ijk && ptr.q==rhs.ptr.q && i==rhs.i && j==rhs.j){return true;}
    else{return false;}
}
bool container_base_2d::iterator_subset::operator!=(const iterator_subset& rhs) const
{
    if(ptr.ijk==rhs.ptr.ijk && ptr.q==rhs.ptr.q && i==rhs.i && j==rhs.j){return false;}
    else{return true;}
}

//Can be dereferenced as an rvalue (if in a dereferenceable state).
c_info& container_base_2d::iterator_subset::operator*(){return ptr;}
c_info* container_base_2d::iterator_subset::operator->(){return &ptr;}

//Can be incremented (if in a dereferenceable state).
//The result is either also dereferenceable or a past-the-end iterator.
//Two iterators that compare equal, keep comparing equal after being both increased.
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator++()
{
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    while(n>0){
        q_++;
        while(q_>=cl_iter->co[ijk_]){
            q_=0;
            next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_++;
            while(q_>=cl_iter->co[ijk_]){
                q_=0;
                next_block_iter(ijk_);
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}
container_base_2d::iterator_subset container_base_2d::iterator_subset::operator++(int)
{
    iterator_subset tmp=*this;
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    while(n>0){
        q_++;
        while(q_>=cl_iter->co[ijk_]){
            q_=0;
            next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_++;
            while(q_>=cl_iter->co[ijk_]){
                q_=0;
                next_block_iter(ijk_);
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return tmp;
}
//Can be decremented (if a dereferenceable iterator value precedes it).
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator--()
{
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    while(n>0){
        q_--;
        while(q_<0){
            previous_block_iter(ijk_);
            q_=cl_iter->co[ijk_]-1;
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_--;
            while(q_<0){
                previous_block_iter(ijk_);
                q_=cl_iter->co[ijk_]-1;
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}
container_base_2d::iterator_subset container_base_2d::iterator_subset::operator--(int)
{
    iterator_subset tmp=*this;
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=1;
    while(n>0){
        q_--;
        while(q_<0){
            previous_block_iter(ijk_);
            q_=cl_iter->co[ijk_]-1;
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_--;
            while(q_<0){
                previous_block_iter(ijk_);
                q_=cl_iter->co[ijk_]-1;
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return tmp;
}
//Supports the arithmetic operators + and - between an iterator and an integer value, or subtracting an iterator from another.
container_base_2d::iterator_subset::difference_type container_base_2d::iterator_subset::operator-(const iterator_subset& rhs) const
{
    difference_type diff=0;
    if(*this==rhs){
        diff=0;
    }
    else if(*this<rhs){
        iterator_subset tmp(*this);
        while(tmp!=rhs){
            tmp++;
            diff++;
        }
        diff=-diff;
    }
    else { //*this>rhs
        iterator_subset tmp(*this);
        while(tmp!=rhs){
            tmp--;
            diff++;
        }
    }
    return diff;
}
container_base_2d::iterator_subset container_base_2d::iterator_subset::operator+(const difference_type& incre) const
{
    iterator_subset tmp=*this;
    tmp+=incre;
    return tmp;
}
container_base_2d::iterator_subset container_base_2d::iterator_subset::operator-(const difference_type& decre) const
{
    iterator_subset tmp=*this;
    tmp-=decre;
    return tmp;
}

//Can be compared with inequality relational operators (<, >, <= and >=).
bool container_base_2d::iterator_subset::operator>(const iterator_subset& rhs) const
{
    if(j==rhs.j && i==rhs.i && ptr.q>rhs.ptr.q){return true;}
    else if(j==rhs.j && i>rhs.i){return true;}
    else if(j>rhs.j){return true;}
    else {return false;}
}
bool container_base_2d::iterator_subset::operator<(const iterator_subset& rhs) const
{
    if(j==rhs.j && i==rhs.i && ptr.q<rhs.ptr.q){return true;}
    else if(j==rhs.j && i<rhs.i){return true;}
    else if(j<rhs.j){return true;}
    else {return false;}
}
bool container_base_2d::iterator_subset::operator>=(const iterator_subset& rhs) const
{
    if(j==rhs.j && i==rhs.i && ptr.q>=rhs.ptr.q){return true;}
    else if(j==rhs.j && i>rhs.i){return true;}
    else if(j>rhs.j){return true;}
    else {return false;}
}
bool container_base_2d::iterator_subset::operator<=(const iterator_subset& rhs) const
{
    if(j==rhs.j && i==rhs.i && ptr.q<=rhs.ptr.q){return true;}
    else if(j==rhs.j && i<rhs.i){return true;}
    else if(j<rhs.j){return true;}
    else {return false;}
}

//Supports compound assignment operations += and -=
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator+=(const difference_type& incre)
{
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=incre;
    while(n>0){
        q_++;
        while(q_>=cl_iter->co[ijk_]){
            q_=0;
            next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_++;
            while(q_>=cl_iter->co[ijk_]){
                q_=0;
                next_block_iter(ijk_);
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}
container_base_2d::iterator_subset& container_base_2d::iterator_subset::operator-=(const difference_type& decre)
{
    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=decre;
    while(n>0){
        q_--;
        while(q_<0){
            previous_block_iter(ijk_);
            q_=cl_iter->co[ijk_]-1;
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_--;
            while(q_<0){
                previous_block_iter(ijk_);
                q_=cl_iter->co[ijk_]-1;
            }
        }
        n--;
    }
    ptr.set(ijk_,q_);
    return *this;
}
//Supports the offset dereference operator ([])
c_info& container_base_2d::iterator_subset::operator[](const difference_type& incre) const
{
    c_info ci;
    iterator_subset tmp=*this;

    int q_=ptr.q; int ijk_=ptr.ijk;
    int n=incre;
    while(n>0){
        q_++;
        while(q_>=cl_iter->co[ijk_]){
            q_=0;
            tmp.next_block_iter(ijk_);
        }
        while(cl_iter->mode!=no_check&&out_of_bounds_iter(ijk_,q_)){
            q_++;
            while(q_>=cl_iter->co[ijk_]){
                q_=0;
                tmp.next_block_iter(ijk_);
            }
        }
        n--;
    }
    ci.set(ijk_,q_);

    return ci;
}
//swappable??

container_base_2d::iterator_subset container_base_2d::begin(subset_info& si) {
    return iterator_subset(&si); //first particle in the container
}

container_base_2d::iterator_subset container_base_2d::end(subset_info& si) {
    c_info cinfo;
    //find the last particle to point to
    int i_=si.bi; int j_=si.bj;
    int ci_=si.ddi; int cj_=si.ddj;
    int ijk_=si.ddi+si.nx*si.ddj;
    int q_=si.co[ijk_]-1;
    double px_=si.aapx; double py_=si.aapy;

    while(q_<0){
        si.previous_block_iter2(ijk_,i_,j_,ci_,cj_,px_,py_);
        q_=si.co[ijk_]-1;
    }
    while(si.mode!=no_check&&si.out_of_bounds_2(ijk_,q_,px_,py_)){
        q_--;
        while(q_<0){
            si.previous_block_iter2(ijk_,i_,j_,ci_,cj_,px_,py_);
            q_=si.co[ijk_]-1;
        }
    }
    //printf("c, %d %d, %d %d %d %d %d %d, %d %d %d\n", ijk_,q_, ai, bi, aj, bj, ak, bk,i_, j_, k_);
    cinfo.set(ijk_,q_+1); //1 over the end
    return iterator_subset(&si,cinfo,i_,j_);
}

//----------------iterator_order---------------------

//constructor: point to first particle in container
container_base_2d::iterator_order::iterator_order(particle_order& vo_) : cp_iter(vo_.o), ptr_n(0) {
    //find the first particle to point to
    int ijk_=cp_iter[0];
    int q_=cp_iter[1];
    ptr.set(ijk_,q_);
}
//copy-assignable
container_base_2d::iterator_order& container_base_2d::iterator_order::operator=(iterator_order other){
    cp_iter=other.cp_iter;
    ptr_n=other.ptr_n;
    ptr=other.ptr;
    return *this;
}

//Can be compared for equivalence using the equality/inequality operators
bool container_base_2d::iterator_order::operator==(const iterator_order& rhs) const {
    if(ptr_n==rhs.ptr_n){return true;}
    else{return false;}
}
bool container_base_2d::iterator_order::operator!=(const iterator_order& rhs) const {
    if(ptr_n==rhs.ptr_n){return false;}
    else{return true;}
}

//Can be dereferenced as an rvalue (if in a dereferenceable state).
c_info& container_base_2d::iterator_order::operator*() {
    return ptr;
}
c_info* container_base_2d::iterator_order::operator->() {
    return &ptr;
}

//Can be incremented (if in a dereferenceable state).
//The result is either also dereferenceable or a past-the-end iterator.
//Two iterators that compare equal, keep comparing equal after being both increased.
container_base_2d::iterator_order& container_base_2d::iterator_order::operator++() {
    ptr_n++;
    int ijk_=cp_iter[2*ptr_n];
    int q_=cp_iter[2*ptr_n+1];
    ptr.set(ijk_,q_);
    return *this;
}
container_base_2d::iterator_order container_base_2d::iterator_order::operator++(int) {
    iterator_order tmp(*this);
    ptr_n++;
    int ijk_=cp_iter[2*ptr_n];
    int q_=cp_iter[2*ptr_n+1];
    ptr.set(ijk_,q_);
    return tmp;
}

//Can be decremented (if a dereferenceable iterator value precedes it).
container_base_2d::iterator_order& container_base_2d::iterator_order::operator--() {
    ptr_n--;
    int ijk_=cp_iter[2*ptr_n];
    int q_=cp_iter[2*ptr_n+1];
    ptr.set(ijk_,q_);
    return *this;
}
container_base_2d::iterator_order container_base_2d::iterator_order::operator--(int) {
    iterator_order tmp(*this);
    ptr_n--;
    int ijk_=cp_iter[2*ptr_n];
    int q_=cp_iter[2*ptr_n+1];
    ptr.set(ijk_,q_);
    return tmp;
}

//Supports the arithmetic operators + and - between an iterator and an integer value, or subtracting an iterator from another.
container_base_2d::iterator_order::difference_type container_base_2d::iterator_order::operator-(const iterator_order& rhs) const {
    difference_type diff=ptr_n-rhs.ptr_n;
    return diff;
}
container_base_2d::iterator_order container_base_2d::iterator_order::operator+(const difference_type& incre) const {
    iterator_order tmp(*this);
    tmp+=incre;
    return tmp;
}
container_base_2d::iterator_order container_base_2d::iterator_order::operator-(const difference_type& decre) const {
    iterator_order tmp(*this);
    tmp-=decre;
    return tmp;
}

//Can be compared with inequality relational operators (<, >, <= and >=).
bool container_base_2d::iterator_order::operator>(const iterator_order& rhs) const {
    if(ptr_n>rhs.ptr_n){return true;}
    else{return false;}
}
bool container_base_2d::iterator_order::operator<(const iterator_order& rhs) const {
    if(ptr_n<rhs.ptr_n){return true;}
    else{return false;}
}
bool container_base_2d::iterator_order::operator>=(const iterator_order& rhs) const {
    if(ptr_n>=rhs.ptr_n){return true;}
    else{return false;}
}
bool container_base_2d::iterator_order::operator<=(const iterator_order& rhs) const {
    if(ptr_n<=rhs.ptr_n){return true;}
    else{return false;}
}

//Supports compound assignment operations += and -=
container_base_2d::iterator_order& container_base_2d::iterator_order::operator+=(const difference_type& incre) {
    ptr_n+=incre;
    int ijk_=cp_iter[2*ptr_n];
    int q_=cp_iter[2*ptr_n+1];
    ptr.set(ijk_,q_);
    return *this;
}
container_base_2d::iterator_order& container_base_2d::iterator_order::operator-=(const difference_type& decre) {
    ptr_n-=decre;
    int ijk_=cp_iter[2*ptr_n];
    int q_=cp_iter[2*ptr_n+1];
    ptr.set(ijk_,q_);
    return *this;
}

//Supports the offset dereference operator ([])
c_info& container_base_2d::iterator_order::operator[](const difference_type& incre) const {
    c_info ci;
    int ci_n=ptr_n+incre;
    int ijk_=cp_iter[2*ci_n];
    int q_=cp_iter[2*ci_n+1];
    ci.set(ijk_,q_);
    return ci;
}

//swappable??

//first particle in the container
container_base_2d::iterator_order container_base_2d::begin(particle_order &vo) {return iterator_order(vo);}
container_base_2d::iterator_order container_base_2d::end(particle_order &vo)
{    //vo, ptr, n
    int ptr_n_=0.5*(vo.op-vo.o); //1-over-the-last-particle, eg. if 0,1,2,3,4 particle, here, ptr_n=5
    c_info ci;
    int ijk_=-1; //dummy
    int q_=-1; //dummy
    ci.set(ijk_,q_);
    return iterator_order(vo, ci, ptr_n_);
}
