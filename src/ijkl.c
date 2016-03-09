#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i, j, k, l;
    int N, eri = 0;

    printf("N = ");
    if (scanf("%d", &N) != 1 || N <= 0) {
        printf("Error!\n");
        exit(1);
    }
    
    printf("< i j | k l >\n");
    printf("-------------\n");
    for (i = 0; i < N; i++)
        for (j = i; j < N; j++)
            for (k = i; k < N; k++)
                for (l = (k == i) ? j : k; l < N; l++) {
                    printf("< %d %d | %d %d >\n", i+1, j+1, k+1, l+1);
                    eri++;
                }

    printf("\n------------------\n");
    printf("       N = %d\n", N);
    printf("  N(ERI) = %d\n", eri);
    printf("   Exact = %d\n", (N*N*N*N + 2*N*N*N + 3*N*N + 2*N)/8);
    printf("   N^4/8 = %d\n", (int)((double)N*N*N*N / 8.0));
}


