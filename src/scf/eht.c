/*    eht.c
 * 
 *  Расширенный метод Хюккеля (EHT).
 *  Служит для генерации начального приближения к молекулярным
 *  орбиталям.
 * 
 *  Идея метода состоит в том, чтобы приблизительно оценить матричные
 *  элементы оператора Фока в базисе АО, а затем диагонализацией F найти
 *  начальное разложение МО по АО, т.е. решить вариационную задачу
 *       FC = ESC
 *  Матрица S при этом рассчитывается "по-честному". Полученные векторы
 *  C используются для генерации начальной матрицы плотности метода ХФ.
 *  
 *  Будем использовать для оценки элементов F значения потенциалов
 *  ионизации АО, полученные на уровне теории HF/cc-pVTZ, то есть,
 *  фактически, энергии на хартри-фоковском пределе (сам не знаю,
 *  зачем мне базис такого качества!). Кроме того, при их расчете
 *  я не переусложнял задачу и пользовался UHF (разумеется,
 *  CASSCF был бы для этого гораздо лучше!). Однако учитывая, сколько в
 *  методе Хюккеля глубокого смысла, такое огрубление не должно сильно
 *  испортить результаты нашего расчета.
 *  Будет время - пересчитаю в CASSCF.
 *  
 *  Главная проблема, которая нас подстерегает, состоит в том, чтобы
 *  определить главное квантовое число орбитали гауссова типа, исходя
 *  только из ее показателя экспоненты.
 * 
 *  В дальнейшем планируется для генерации EHT-начального приближения не
 *  брать табличные значения ПИ, а рассчитывать их "на лету" для данного
 *  базисного набора.
 */

#include <stdio.h>
#include <stdlib.h>

#include "scf.h"
#include "../input/chem.h"
#include "../input/basis.h"
#include "../util/util.h"

/* Ionization potentials
 * Two-dimensional array of format:
 * [Z][] = {1s 2s 2p 3s 3p 3d ...}
 */
double IP[][5] = {
/* H  */  { -0.50,  0.00,  0.00, 0.0, 0.0},
/* He */  { -0.91,  0.00,  0.00, 0.0, 0.0},
/* Li */  { -2.48, -0.20,  0.00, 0.0, 0.0},
/* Be */  { -4.73, -0.31,  0.00, 0.0, 0.0},
/* B  */  { -7.69, -0.50, -0.32, 0.0, 0.0},
/* C  */  {-11.32, -0.71, -0.44, 0.0, 0.0},
/* N  */  {-15.62, -0.94, -0.57, 0.0, 0.0},
/* O  */  {-20.67, -1.24, -0.67, 0.0, 0.0},
/* F  */  {-26.38, -1.57, -0.76, 0.0, 0.0},
/* Ne */  {-32.77, -1.93, -0.85, 0.0, 0.0},
};

/*    Fii = IP[i]
 *    Fij = 1/2*K*(Fii + Fjj)*Sij, K = 1.75
 */
void guess_F_eht(double *F, double *S, struct basis_function *bfn, int N)
{
	int i, j;
	
	for (i = 0; i < N; i++) {
		int n = bfn[i].f->n;
		int l = bfn[i].f->L;
		int z = bfn[i].a->Z;
		F[i*N+i] = IP[z-1][n-1+l];
	}
	for (i = 0; i < N; i++)
		for (j = i+1; j < N; j++) {
			F[i*N+j] = 0.875 * (F[i*N+i] + F[j*N+j]) * S[i*N+j];
			F[j*N+i] = F[i*N+j];
		}
}








