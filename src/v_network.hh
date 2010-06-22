class voronoi_network {
	public:	
		fpoint **pts;
		int **idmem;
		int *ptsc;
		int *ptsmem;
		int **ed;
		int **ne;
		double **raded;
		unsigned int **pered;
		int edc,edmem;
		int *nu;
		int *numem;
		int *reg;
		int *regp;
		int *vmap;
		unsigned int *perio;
		int netmem;
		void print_network(ostream &os=cout,bool reverse_remove=false);
		inline void print_network(const char* filename) {
			ofstream os;os.open(filename,ofstream::out|ofstream::trunc);
			print_network(os);os.close();
		}
		void draw_network(ostream &os=cout);
		inline void draw_network(const char* filename) {
			ofstream os;os.open(filename,ofstream::out|ofstream::trunc);
			draw_network(os);os.close();
		}
		template<class n_option>
		void add_to_network_slanted(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z,int idn);	
		template<class n_option>
		void add_to_network(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z,int idn);	
	private:
		void add_particular_vertex_memory(int l);
		void add_edge_network_memory();
		void add_network_memory(int l);
		inline unsigned int pack_periodicity(int i,int j,int k);
		inline void unpack_periodicity(unsigned int pa,int &i,int &j,int &k);
		template<class n_option>
		void add_edges_to_network(voronoicell_base<n_option> &c);
		bool not_already_there(int k,int j,unsigned int cper);
		bool search_previous_slanted(fpoint gx,fpoint gy,fpoint x,fpoint y,fpoint z,int &ijk,int &q,unsigned int &cper);
		bool safe_search_previous(fpoint x,fpoint y,fpoint z,int &ijk,int &q,unsigned int &cper);
		bool search_previous(fpoint x,fpoint y,fpoint z,int &ijk,int &q,unsigned int &cper);		
};
