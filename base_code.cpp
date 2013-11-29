#include "functionality.h"
parms p;

int main(int argc, char* argv[]) {
	
	///All entities are in mm dimensions. (the scaling is done while writing the results in vtk and stl files)
	p.alpha = 0.0;
	p.beta = 0.0;

	p.a = 0.0;
	p.b = 1.0;
	p.c = 0.0;
	p.d = 0.0;
	
	p.rb = 1.0; // bottom radius
	p.rb0= 1.1;
	p.rt = 0.8; // top daughter
	p.h1 = -1.0; // bottom height parent
	p.h2 = 0; //
	p.h3 = 1.5; //10 top height daughter 2 worked
	p.ra = p.rt;

	p.xshift = 0.0;
	p.c1 = 1.0;
	p.c2 = 2;
	p.b_bnw = p.rb;
	p.angle = pi / 1.8; // 1.80 angle between p and d
	p.ss0 = 3;
	p.ss1 = 3;
	
	p.xc = 0.8;

	p.r = 0.95; // radius of intersecting cylinder where 0.636 <= r <=0.9886

	p.s0 = 0.6;
	p.s1 = 0.8;

	p.A1 = p.c1;
	p.B1 = 0;
	p.C1 = 3 * (p.c2 * cos(p.alfa) - p.c1) - (p.ss1 * sin(p.alfa));
	p.D1 = 2 * (p.c2 * cos(p.alfa) - p.c1) + (p.ss1 * sin(p.alfa));

	p.A2 = 0;
	p.B2 = 3;
	p.C2 = (3 * p.c2 * sin(p.alfa)) - (2 * p.ss0) - (p.ss1 * cos(p.alfa));
	p.D2 = (p.ss1 * cos(p.alfa)) + (p.ss0) - (2 * p.c2 * sin(p.alfa));

	p.v1 = acos(p.xc / p.c1);
	p.v2 = acos(-(p.xc - p.r * sin(p.alfa)) / (p.c2 * cos(p.alfa)));

	p.m = 25; // number of axial points >=3
	p.n = 11; //  number of circumferential points <=m && >=3
	int parent_grid_point_skip = 0 ,daughter_grid_point_skip = 0;
	int downstream_offset,upstream_offset;

	p.nx = p.m + 2;
	p.ny = p.n + 2;
	int num_SMCs=0,num_ECs=0;
	double tol = 1e-7;
	int iflag = 4;
	int itcg = 50;
	int mm, nn, lw;



	if (iflag == 2) {
		mm = 3 * p.m;
		nn = 7 * p.n;
		lw = (int) (fmax(mm, nn)) + 2 * (p.n + p.m) + 19; //for iflag = 2
	} else if (iflag == 4) {
		mm = 3 * p.m;
		nn = 4 * p.n;
		lw = (int) (fmax(mm, nn)) + 4 * p.n + 2 * p.m
				+ 0.5 * ((p.n + 1) * (p.n + 1)) + 19;	//for iflag = 4
	}
	lw = 10000;
	double w[lw];
	int idf = p.nx;

	double **fx, **fy, **fz;
	double *bda, *bdb, *bdc, *bdd;
	
	double ****f;
	f = (double****)malloc(2*sizeof(double***));
	f[0] = (double***)malloc(3*sizeof(double***));
	f[1] = (double***)malloc(3*sizeof(double***));
	for (int i=0;i<3;i++){
		f[0][i]=allocate_dirichlet(p.nx,p.ny);
		f[1][i]=allocate_dirichlet(p.nx,p.ny);
	}
	
	fx = allocate_dirichlet(p.nx, p.ny);
	fy = allocate_dirichlet(p.nx, p.ny);
	fz = allocate_dirichlet(p.nx, p.ny);
	bda = allocate_neuman(p.n);
	bdb = allocate_neuman(p.n);
	bdc = allocate_neuman(p.m);
	bdd = allocate_neuman(p.m);

	initialize_f(p, fx, fy, fz);
	initialize_bd(p.n, bda);
	initialize_bd(p.n, bdb);
	initialize_bd(p.m, bdc);
	initialize_bd(p.m, bdd);
	double****storage = (double****) malloc(2 * sizeof(double***));
	for (int i = 0; i < 2; i++)
		storage[i] = (double***) malloc(3 * sizeof(double**));
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			storage[i][j] = (double**) malloc(p.nx * sizeof(double*));
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k=0; k<p.nx; k++)
				storage[i][j][k] = (double*) malloc(p.ny * sizeof(double));
		}
	}
	///Setting up conditions for left half of the geometry
	p.alfa = pi - p.angle;
	p.c = 3 * pi / 2;
	p.d = 2 * pi + pi / 2;

	/***** Solving for parent segment ******/
	dirichlet_boundary_parent_segment(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	downstream_offset = parent_grid_point_skip;upstream_offset =0;
	//num_SMCs+=smc_mesh(p, fx, fy, fz,"smc","parent_left",downstream_offset,upstream_offset);
	//num_ECs+=ec_mesh(p, fx, fy, fz,"ec","parent_left",downstream_offset,upstream_offset);
	store_arrays(p,fx,fy,fz,storage,0);
	format_vtk_unstructuredGrid(p, fx, fy, fz, "parent_left",downstream_offset,upstream_offset);
    
	/***** Solving for half end caps on parent segment *****/
	dirichlet_boundary_end_cap_parent(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "end_cap_parent_left",foldername);
	downstream_offset =0; upstream_offset =0;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "end_cap_parent_left",downstream_offset,upstream_offset);
	/***** Solving for outer wall Left ******/
	dirichlet_boundary_outer_wall(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s0* cos(p.angle);
		bdb[k] = p.s1* sin(p.angle);
	}

	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s0* sin(p.angle);
		bdb[k] = -p.s1* cos(p.angle);
	}

	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "outer_left",foldername);
	downstream_offset =0; upstream_offset =daughter_grid_point_skip;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "outer_left",downstream_offset,upstream_offset);
	//print_stdout( p,fx,fy,fz);
	//record_data(0,p,fx,fy,fz,f);
	//print_stdout( p,&f[0][0][0],&f[0][1][0],&f[0][2][0]);
	//format_vtk(p, &f[0][0][0],&f[0][1][0],&f[0][2][0], "outer_left");
	//num_SMCs+=smc_mesh(p, fx, fy, fz,"smc","outer_left",downstream_offset,upstream_offset);
	//num_ECs+=ec_mesh(p, fx, fy, fz,"ec","outer_left",downstream_offset,upstream_offset);
	
	/***** Solving for half end cap on left daughter segment *****/
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "end_cap_left_daughter_1",foldername);
	downstream_offset = 0; upstream_offset =0;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "end_cap_left_daughter_1",downstream_offset,upstream_offset);
	/***** Solving for inner wall Left ******/
	p.c = pi / 2;
	p.d = 3 * pi / 2;
	dirichlet_boundary_inner_wall(p, fx, fy, fz);
	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s0 * sin(p.angle);
		bdb[k] = p.s1* cos(p.angle);
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s0 * cos(p.angle);
		bdb[k] = p.s1* sin(p.angle);
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "inner_left",foldername);
	downstream_offset =0; upstream_offset =daughter_grid_point_skip;
	//num_SMCs+=smc_mesh(p, fx, fy, fz,"smc","inner_left",downstream_offset,upstream_offset);
	//num_ECs+=ec_mesh(p, fx, fy, fz,"ec","inner_left",downstream_offset,upstream_offset);
	format_vtk_unstructuredGrid(p, fx, fy, fz, "inner_left",downstream_offset,upstream_offset);
	
	/***** Solving for half end cap on left daughter segment *****/
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "end_cap_left_daughter_2",foldername);
	downstream_offset =0; upstream_offset =0;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "end_cap_left_daughter_2",downstream_offset,upstream_offset);

	/**** SOLVING FOR SYMMETRIC RIGHT HALF OF THE GEOMETRY ****/

	p.alfa = pi + p.angle;
	p.c = pi/2;
	p.d = 3*pi/2;

	/***** Solving for parent segment ******/
	dirichlet_boundary_parent_segment(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = p.h1;
		bdb[k] = -p.h1;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "parent_right",foldername);
	downstream_offset =parent_grid_point_skip; upstream_offset =0;
	store_arrays(p,fx,fy,fz,storage,1);
	format_vtk_unstructuredGrid(p, fx, fy, fz, "parent_right",downstream_offset,upstream_offset);
	//num_SMCs+=smc_mesh(p, fx, fy, fz,"smc","parent_right",downstream_offset,upstream_offset);
	//num_ECs+=ec_mesh(p, fx, fy, fz,"ec","parent_right",downstream_offset,upstream_offset);

//	format_primitive(p,store_arrays,"parent", foldername);

	/***** Solving for half end caps on parent segment *****/
	dirichlet_boundary_end_cap_parent(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "end_cap_parent_right", foldername);
	downstream_offset =0; upstream_offset =0;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "end_cap_parent_right",downstream_offset,upstream_offset);
	/***** Solving for outer wall Right******/
	dirichlet_boundary_outer_wall(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = p.s0 * cos(p.angle);
		bdb[k] = -p.s1* sin(p.angle);
	}

	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s0 * sin(p.angle);
		bdb[k] = p.s1* cos(p.angle);
	}

	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "outer_right",foldername);
	downstream_offset =0;upstream_offset =daughter_grid_point_skip;
	//num_SMCs+=smc_mesh(p, fx, fy, fz,"smc","outer_right",downstream_offset,upstream_offset);
	//num_ECs+=ec_mesh(p, fx, fy, fz,"ec","outer_right",downstream_offset,upstream_offset);
	format_vtk_unstructuredGrid(p, fx, fy, fz, "outer_right",downstream_offset,upstream_offset);
	/***** Solving for half end cap on left daughter segment *****/
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "end_cap_right_daughter_1",foldername);
	downstream_offset =0; upstream_offset = 0;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "end_cap_right_daughter_1",downstream_offset,upstream_offset);

	/***** Solving for inner wall Right ******/
	p.c = 3*pi / 2;
	p.d = 2*pi +  pi / 2;
	dirichlet_boundary_inner_wall(p, fx, fy, fz);
	for (int k = 0; k < p.n; k++) {
		bda[k] = p.s0 * sin(p.angle);
		bdb[k] = -p.s1* cos(p.angle);
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = -p.s0 * cos(p.angle);
		bdb[k] = p.s1* sin(p.angle);
	}

		fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);
	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "inner_right",foldername);
	downstream_offset =0; upstream_offset =daughter_grid_point_skip;
	//num_SMCs+=smc_mesh(p, fx, fy, fz,"smc","inner_right",downstream_offset,upstream_offset);
	//num_ECs+=ec_mesh(p, fx, fy, fz,"ec","inner_right",downstream_offset,upstream_offset);
	format_vtk_unstructuredGrid(p, fx, fy, fz, "inner_right",downstream_offset,upstream_offset);
	/***** Solving for half end cap on left daughter segment *****/
	dirichlet_boundary_end_cap_daughter(p, fx, fy, fz);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0.0;
		bdb[k] = 0.0;
	}
	fx = solve(p, bda, bdb, bdc, bdd, fx, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fy = solve(p, bda, bdb, bdc, bdd, fy, idf, iflag, tol, itcg, w, lw);

	for (int k = 0; k < p.n; k++) {
		bda[k] = 0;
		bdb[k] = 0;
	}
	fz = solve(p, bda, bdb, bdc, bdd, fz, idf, iflag, tol, itcg, w, lw);

	//format_vtk(p, fx, fy, fz, "end_cap_right_daughter_2",foldername);
	downstream_offset=0; upstream_offset = 0;
	format_vtk_unstructuredGrid(p, fx, fy, fz, "end_cap_right_daughter_2",downstream_offset,upstream_offset);
/*********************************************/
	printf("TOTAL SMCs =%d\nTOTAL ECs=%d\n",num_SMCs,num_ECs);


}
