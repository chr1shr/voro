// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file common.cc
 * \brief Implementations of the small helper functions. */

#include "common.hh"

namespace voro {

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void voro_fatal_error(const char *p,int status) {
    fprintf(stderr,"voro++: %s\n",p);
    exit(status);
}

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_positions_2d(std::vector<double> &v,FILE *fp) {
    if(v.size()>0) {
        fprintf(fp,"(%g,%g)",v[0],v[1]);
        for(unsigned int k=2;k<v.size();k+=2)
            fprintf(fp," (%g,%g)",v[k],v[k+1]);
    }
}

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] pr the precision to print the numbers to.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_positions_2d(int pr,std::vector<double> &v,FILE *fp) {
    if(v.size()>0) {
        fprintf(fp,"(%.*g,%.*g)",pr,v[0],pr,v[1]);
        for(unsigned int k=2;k<v.size();k+=2)
            fprintf(fp," (%.*g,%.*g)",pr,v[k],pr,v[k+1]);
    }
}

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_positions_3d(std::vector<double> &v,FILE *fp) {
    if(v.size()>0) {
        fprintf(fp,"(%g,%g,%g)",v[0],v[1],v[2]);
        for(unsigned int k=3;k<v.size();k+=3)
            fprintf(fp," (%g,%g,%g)",v[k],v[k+1],v[k+2]);
    }
}

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] pr the precision to print the numbers to.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_positions_3d(int pr,std::vector<double> &v,FILE *fp) {
    if(v.size()>0) {
        fprintf(fp,"(%.*g,%.*g,%.*g)",pr,v[0],pr,v[1],pr,v[2]);
        for(unsigned int k=3;k<v.size();k+=3)
            fprintf(fp," (%.*g,%.*g,%.*g)",pr,v[k],pr,v[k+1],pr,v[k+2]);
    }
}

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation was
 * successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char *filename,const char *mode) {
    FILE *fp=fopen(filename,mode);
    if(fp==NULL) {
        fprintf(stderr,"voro++: Unable to open file '%s'\n",filename);
        exit(VOROPP_FILE_ERROR);
    }
    return fp;
}

/** \brief Prints a vector of integers.
 *
 * Prints a vector of integers.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_vector(std::vector<int> &v,FILE *fp) {
    int k=0,s=v.size();
    while(k+4<s) {
        fprintf(fp,"%d %d %d %d ",v[k],v[k+1],v[k+2],v[k+3]);
        k+=4;
    }
    if(k+3<=s) {
        if(k+4==s) fprintf(fp,"%d %d %d %d",v[k],v[k+1],v[k+2],v[k+3]);
        else fprintf(fp,"%d %d %d",v[k],v[k+1],v[k+2]);
    } else {
        if(k+2==s) fprintf(fp,"%d %d",v[k],v[k+1]);
        else fprintf(fp,"%d",v[k]);
    }
}

/** \brief Prints a vector of doubles.
 *
 * Prints a vector of doubles.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_vector(std::vector<double> &v,FILE *fp) {
    int k=0,s=v.size();
    while(k+4<s) {
        fprintf(fp,"%g %g %g %g ",v[k],v[k+1],v[k+2],v[k+3]);
        k+=4;
    }
    if(k+3<=s) {
        if(k+4==s) fprintf(fp,"%g %g %g %g",v[k],v[k+1],v[k+2],v[k+3]);
        else fprintf(fp,"%g %g %g",v[k],v[k+1],v[k+2]);
    } else {
        if(k+2==s) fprintf(fp,"%g %g",v[k],v[k+1]);
        else fprintf(fp,"%g",v[k]);
    }
}

/** \brief Prints a vector of doubles.
 *
 * Prints a vector of doubles.
 * \param[in] pr the precision to print the numbers to.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_vector(int pr,std::vector<double> &v,FILE *fp) {
    int k=0,s=v.size();
    while(k+4<s) {
        fprintf(fp,"%.*g %.*g %.*g %.*g ",pr,v[k],pr,v[k+1],pr,v[k+2],pr,v[k+3]);
        k+=4;
    }
    if(k+3<=s) {
        if(k+4==s) fprintf(fp,"%.*g %.*g %.*g %.*g",pr,v[k],pr,v[k+1],pr,v[k+2],pr,v[k+3]);
        else fprintf(fp,"%.*g %.*g %.*g",pr,v[k],pr,v[k+1],pr,v[k+2]);
    } else {
        if(k+2==s) fprintf(fp,"%.*g %.*g",pr,v[k],pr,v[k+1]);
        else fprintf(fp,"%.*g",pr,v[k]);
    }
}

/** \brief Prints a vector a face vertex information.
 *
 * Prints a vector of face vertex information. A value is read, which
 * corresponds to the number of vertices in the next face. The routine reads
 * this number of values and prints them as a bracked list. This is repeated
 * until the end of the vector is reached.
 * \param[in] v the vector to interpret and print.
 * \param[in] fp the file stream to print to. */
void voro_print_face_vertices(std::vector<int> &v,FILE *fp) {
    int j,k=0,l;
    if(v.size()>0) {
        l=v[k++];
        if(l<=1) {
            if(l==1) fprintf(fp,"(%d)",v[k++]);
            else fputs("()",fp);
        } else {
            j=k+l;
            fprintf(fp,"(%d",v[k++]);
            while(k<j) fprintf(fp,",%d",v[k++]);
            fputs(")",fp);
        }
        while((unsigned int) k<v.size()) {
            l=v[k++];
            if(l<=1) {
                if(l==1) fprintf(fp," (%d)",v[k++]);
                else fputs(" ()",fp);
            } else {
                j=k+l;
                fprintf(fp," (%d",v[k++]);
                while(k<j) fprintf(fp,",%d",v[k++]);
                fputs(")",fp);
            }
        }
    }
}

/** \brief Checks whether neighbor information is in a format string.
 *
 * Checks to see whether "%n" appears in a format sequence to determine
 * whether neighbor information is required or not.
 * \param[in] format the format string to check.
 * \return True if a "%n" is found, false otherwise. */
bool voro_contains_neighbor(const char *format) {
    char *fmp=(const_cast<char*>(format));

    // Check to see if "%n" appears in the format sequence
    while(*fmp!=0) {
        if(*fmp=='%') {
            fmp++;
            if(*fmp=='n') return true;
            else if(*fmp==0) return false;
        }
        fmp++;
    }

    return false;
}

/** \brief Reads the precision from a custom output format string.
 *
 * The 2D and 3D Voronoi cell classes can output custom information
 * about a Voronoi cell to varying precision. This common function
 * reads the precision from an output string.
 * \param[in] fp the file handle to write to.
 * \param[in,out] fmp a pointer to the position in the output string before the
 *                    numerical precision value; on exit, this is updated to
 *                    to the position after the numerical value.
 * \param[out] pr the precision that was read.
 * \return True if the precision was successfully read, false otherwise. */
bool voro_read_precision(FILE *fp,char *&fmp,int &pr) {
    fmp++;

    // Check for first digit of the precision string
    if(*fmp==0) {fputs("%.",fp);fmp--;return false;}
    else if(*fmp>='0'&&*fmp<='9') {
        pr=static_cast<int>(*fmp-'0');
    } else {fputs("%.",fp);putc(*fmp,fp);return false;}

    // Check for subsequent digits of the precision string
    fmp++;
    while(*fmp>='0'&&*fmp<='9') {
        pr=10*pr+static_cast<int>(*fmp-'0');
        fmp++;
    }
    return true;
}

}
