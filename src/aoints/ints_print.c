#include <mpi.h>
#include <stdio.h>

#include "aoints.h"
#include "scf.h"

static struct basis_function *basis_funcs;
static int n_of_bfns;

/* ignore zero value and elements under diagonal */
void print_ao1eints(double (*eval_func)(), char *name)
{
	int i, j, n1, n2;
	struct basis_function *bfn1, *bfn2;
	static char *aosym[] = {"s", "px", "py", "pz"};
	
	for (i = 0; i < n_of_bfns; i++)
		for (j = i; j < n_of_bfns; j++) {
			char *aos1, *aos2;
			
			bfn1 = &basis_funcs[i];
			bfn2 = &basis_funcs[j];
			aos1 = aosym[bfn1->f->L + bfn1->m];
			aos2 = aosym[bfn2->f->L + bfn2->m];
			
			printf("    %s   %-4s%-4d%-6s%-4s%-4d%-2s%15.8f\n", name,
				searchByZ(bfn1->a->Z)->sym, i, aos1,
				searchByZ(bfn2->a->Z)->sym, j, aos2,
				eval_func(bfn1, bfn2));
		}
}

void print_ints(struct basis_function *funcs, int N)
{
	int mpi_rank;	
	
	n_of_bfns = N;
	basis_funcs = funcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	if (scf_options.print_1eov) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0) {
			printf("\n\n    ==============================================\n");
			printf("           Overlap One-Electron AO Integrals \n");
			printf("    ==============================================\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		print_ao1eints(aoint_overlap, "1eov");
	}
	
	if (scf_options.print_1eke) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0) {
			printf("\n\n    ==============================================\n");
			printf("       Kinetic Energy One-Electron AO Integrals \n");
			printf("    ==============================================\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		print_ao1eints(aoint_kinetic, "1eke");
	}
	
	if (scf_options.print_1epe) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0) {
			printf("\n\n    ===============================================\n");
			printf("      Nuclear-Attraction One-Electron AO Integrals \n");
			printf("    ===============================================\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		print_ao1eints(aoint_potential, "1epe");
	}
	
	// recalculation of ERI's on master-node specially for printing out
	if (scf_options.print_2eri && mpi_rank == 0) {
		int i, j, k, l, nzero = 0;
		int M = n_of_bfns;
		double eri;
		
		printf("\n\n    ===============================================\n");
		printf("            Electron-Repulsion AO Integrals        \n");
		printf("    ===============================================\n");
		for (i = 0; i < M; i++)
			for (j = i; j < M; j++)
				for (k = i; k < M; k++)
					for (l = (k == i) ? j : k; l < M; l++) {
						struct basis_function *fi = &basis_funcs[i],
											  *fj = &basis_funcs[j],
											  *fk = &basis_funcs[k],
											  *fl = &basis_funcs[l];
						eri = aoint_eri(fi, fj, fk, fl);
						if (eri < 1e-8)
							nzero++;
						printf("            %d   %d   %d   %d   %13.8f\n", i+1, j+1, k+1, l+1, eri);
					}
		printf("    ===============================================\n");
		printf("          Number of unique ERIs = %d\n", (M*M*M*M+2*M*M*M+3*M*M+2*M)/8);
		printf("              Estimated (N^4/8) = %d\n", M*M*M*M/8);
		printf("                 Less than 1e-8 = %d\n", nzero);
		printf("    ===============================================\n");
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
}

